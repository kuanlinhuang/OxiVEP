use axum::extract::State;
use axum::http::header;
use axum::response::IntoResponse;
use axum::Json;
use serde::Deserialize;
use std::path::PathBuf;
use std::sync::{Arc, RwLock};
use std::time::Instant;

use crate::context::AnnotationContext;
use crate::errors::AppError;

/// Shared application state.
pub struct SharedState {
    pub ctx: RwLock<AnnotationContext>,
    pub data_dir: Option<PathBuf>,
    pub sa_dir: Option<PathBuf>,
}

pub type AppState = Arc<SharedState>;

const INDEX_HTML: &str = include_str!("../../../web/index.html");

pub async fn index_html() -> impl IntoResponse {
    (
        [(header::CONTENT_TYPE, "text/html; charset=utf-8")],
        INDEX_HTML,
    )
}

/// Enhanced status: reports transcript/genome/SA state so the SPA can
/// decide whether to use server-side annotation or upload example GFF3.
pub async fn status(State(state): State<AppState>) -> Json<serde_json::Value> {
    let (transcripts, gff3_source, has_fasta, sa_sources) = {
        let guard = state.ctx.read().unwrap();
        (
            guard.transcript_count(),
            guard.gff3_source.clone(),
            guard.seq_provider.is_some(),
            guard.sa_source_names(),
        )
    };
    Json(serde_json::json!({
        "status": "ok",
        "backend": true,
        "transcripts": transcripts,
        "gff3_source": gff3_source,
        "has_fasta": has_fasta,
        "sa_sources": sa_sources,
    }))
}

/// List available genome GFF3 files from the data directory.
/// Scans for .gff3, .gff3.gz, and .fastvep.cache files.
pub async fn list_genomes(State(state): State<AppState>) -> Json<serde_json::Value> {
    let Some(ref data_dir) = state.data_dir else {
        return Json(serde_json::json!({ "genomes": [] }));
    };

    let mut genomes = Vec::new();

    // Scan subdirectories: each subdir is a genome with GFF3 + optional FASTA
    if let Ok(entries) = std::fs::read_dir(data_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                let name = path
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_default();
                let gff3 = find_file_with_ext(&path, &["gff3", "gff3.gz", "fastvep.cache"]);
                let fasta = find_file_with_ext(&path, &["fa", "fasta", "fa.gz", "fasta.gz"]);
                let has_sa = path.join("sa").is_dir()
                    && std::fs::read_dir(path.join("sa"))
                        .map(|rd| {
                            rd.flatten().any(|e| {
                                let n = e.file_name().to_string_lossy().to_string();
                                n.ends_with(".osa") || n.ends_with(".osa2")
                            })
                        })
                        .unwrap_or(false);
                if gff3.is_some() {
                    genomes.push(serde_json::json!({
                        "name": name,
                        "has_fasta": fasta.is_some(),
                        "has_sa": has_sa,
                    }));
                }
            }
        }
    }

    // Also scan top-level for loose GFF3 files
    if let Ok(entries) = std::fs::read_dir(data_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_file() {
                let fname = path
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_default();
                if fname.ends_with(".gff3")
                    || fname.ends_with(".gff3.gz")
                    || fname.ends_with(".fastvep.cache")
                {
                    let stem = fname
                        .trim_end_matches(".fastvep.cache")
                        .trim_end_matches(".gz")
                        .trim_end_matches(".gff3")
                        .to_string();
                    genomes.push(serde_json::json!({
                        "name": stem,
                        "has_fasta": false,
                    }));
                }
            }
        }
    }

    genomes.sort_by(|a, b| {
        a["name"]
            .as_str()
            .unwrap_or("")
            .cmp(b["name"].as_str().unwrap_or(""))
    });

    Json(serde_json::json!({ "genomes": genomes }))
}

/// Load a genome from the data directory by name.
#[derive(Deserialize)]
pub struct LoadGenomeRequest {
    name: String,
}

pub async fn load_genome(
    State(state): State<AppState>,
    Json(req): Json<LoadGenomeRequest>,
) -> Result<Json<serde_json::Value>, AppError> {
    let Some(ref data_dir) = state.data_dir else {
        return Err(AppError::BadRequest("No data directory configured".into()));
    };

    let (gff3_path, fasta_path) = resolve_genome_paths(data_dir, &req.name)?;

    // Check for per-genome SA directory (data_dir/<name>/sa/)
    let genome_sa_dir = data_dir.join(&req.name).join("sa");
    let sa_dir = if genome_sa_dir.is_dir() {
        Some(genome_sa_dir)
    } else {
        // Fall back to global --sa-dir
        state.sa_dir.clone()
    };

    let name = req.name.clone();
    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let (transcripts, sa_sources) = tokio::task::spawn_blocking(move || {
        let mut guard = ctx
            .ctx
            .write()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        guard.load_genome(
            gff3_path.to_str().unwrap(),
            fasta_path.as_ref().map(|p| p.to_str().unwrap()),
            sa_dir.as_ref().map(|p| p.to_str().unwrap()),
        )?;
        let tr_count = guard.transcript_count();
        let sa_names = guard.sa_source_names();
        Ok::<_, anyhow::Error>((tr_count, sa_names))
    })
    .await??;

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "name": name,
        "transcripts": transcripts,
        "sa_sources": sa_sources,
        "time_ms": time_ms,
    })))
}

#[derive(Deserialize)]
pub struct AnnotateRequest {
    vcf: Option<String>,
    pick: Option<bool>,
}

pub async fn annotate(
    State(state): State<AppState>,
    Json(req): Json<AnnotateRequest>,
) -> Result<Json<serde_json::Value>, AppError> {
    let vcf_text = req.vcf.unwrap_or_default();
    if vcf_text.is_empty() {
        return Err(AppError::BadRequest("No VCF data provided".into()));
    }

    let pick = req.pick.unwrap_or(false);
    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let results = tokio::task::spawn_blocking(move || {
        let guard = ctx
            .ctx
            .read()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        guard.annotate_vcf_text(&vcf_text, pick)
    })
    .await??;

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "results": results,
        "count": results.len(),
        "time_ms": time_ms,
    })))
}

pub async fn upload_gff3(
    State(state): State<AppState>,
    body: String,
) -> Result<Json<serde_json::Value>, AppError> {
    if body.is_empty() {
        return Err(AppError::BadRequest("No GFF3 data provided".into()));
    }

    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let (genes, transcripts) = tokio::task::spawn_blocking(move || {
        let mut guard = ctx
            .ctx
            .write()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        guard.update_gff3_text(&body)
    })
    .await??;

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "genes": genes,
        "transcripts": transcripts,
        "time_ms": time_ms,
    })))
}

// --- helpers ---

fn find_file_with_ext(dir: &std::path::Path, extensions: &[&str]) -> Option<PathBuf> {
    let entries = std::fs::read_dir(dir).ok()?;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_file() {
            let fname = path.file_name()?.to_string_lossy();
            for ext in extensions {
                if fname.ends_with(ext) {
                    return Some(path);
                }
            }
        }
    }
    None
}

fn resolve_genome_paths(
    data_dir: &std::path::Path,
    name: &str,
) -> Result<(PathBuf, Option<PathBuf>), AppError> {
    // Check if it's a subdirectory
    let subdir = data_dir.join(name);
    if subdir.is_dir() {
        let gff3 = find_file_with_ext(&subdir, &["gff3", "gff3.gz", "fastvep.cache"])
            .ok_or_else(|| AppError::BadRequest(format!("No GFF3 found in genome '{}'", name)))?;
        let fasta = find_file_with_ext(&subdir, &["fa", "fasta", "fa.gz", "fasta.gz"]);
        return Ok((gff3, fasta));
    }

    // Check top-level files
    for ext in &["gff3", "gff3.gz", "fastvep.cache"] {
        let path = data_dir.join(format!("{}.{}", name, ext));
        if path.exists() {
            return Ok((path, None));
        }
    }

    Err(AppError::BadRequest(format!(
        "Genome '{}' not found in data directory",
        name
    )))
}
