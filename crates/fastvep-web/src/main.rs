mod context;
mod errors;
mod handlers;

use axum::extract::DefaultBodyLimit;
use axum::routing::{get, post};
use axum::Router;
use clap::Parser;
use std::path::PathBuf;
use std::sync::{Arc, RwLock};
use tower::limit::ConcurrencyLimitLayer;
use tower::ServiceBuilder;
use tower_http::cors::{Any, CorsLayer};
use tower_http::trace::TraceLayer;

use crate::context::AnnotationContext;
use crate::handlers::{AppState, SharedState};

#[derive(Parser)]
#[command(name = "fastvep-web")]
#[command(about = "fastVEP production web server for variant annotation")]
struct Cli {
    /// Port to listen on
    #[arg(long, default_value_t = 8080, env = "FASTVEP_PORT")]
    port: u16,

    /// Bind address
    #[arg(long, default_value = "0.0.0.0", env = "FASTVEP_BIND")]
    bind: String,

    /// GFF3 annotation file to load at startup
    #[arg(long, env = "FASTVEP_GFF3")]
    gff3: Option<String>,

    /// FASTA reference file (requires .fai index for mmap)
    #[arg(long, env = "FASTVEP_FASTA")]
    fasta: Option<String>,

    /// Supplementary annotation directory (.osa/.osa2 files)
    #[arg(long, env = "FASTVEP_SA_DIR")]
    sa_dir: Option<String>,

    /// Directory containing genome data (subdirs with GFF3 + FASTA per organism)
    #[arg(long, env = "FASTVEP_DATA_DIR")]
    data_dir: Option<String>,

    /// Upstream/downstream distance in bp
    #[arg(long, default_value_t = 5000)]
    distance: u64,

    /// Maximum request body size in bytes
    #[arg(long, default_value_t = 10_485_760)]
    max_body_size: usize,

    /// Maximum concurrent annotation requests
    #[arg(long, default_value_t = 64)]
    max_concurrent: usize,
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "fastvep_web=info".into()),
        )
        .init();

    let cli = Cli::parse();

    tracing::info!("Loading annotation data...");
    let ctx = AnnotationContext::new(
        cli.gff3.as_deref(),
        cli.fasta.as_deref(),
        cli.sa_dir.as_deref(),
        cli.distance,
    )?;

    let data_dir = cli.data_dir.map(PathBuf::from);
    let sa_dir = cli.sa_dir.as_ref().map(PathBuf::from);
    if let Some(ref dir) = data_dir {
        tracing::info!("Data directory: {}", dir.display());
    }
    if let Some(ref dir) = sa_dir {
        tracing::info!("SA directory: {}", dir.display());
    }

    let state: AppState = Arc::new(SharedState {
        ctx: RwLock::new(ctx),
        data_dir,
        sa_dir,
    });

    let app = Router::new()
        .route("/", get(handlers::index_html))
        .route("/index.html", get(handlers::index_html))
        .route("/api/status", get(handlers::status))
        .route("/api/genomes", get(handlers::list_genomes))
        .route("/api/load-genome", post(handlers::load_genome))
        .route("/api/annotate", post(handlers::annotate))
        .route("/api/upload-gff3", post(handlers::upload_gff3))
        .with_state(state)
        .layer(DefaultBodyLimit::max(cli.max_body_size))
        .layer(
            ServiceBuilder::new()
                .layer(TraceLayer::new_for_http())
                .layer(
                    CorsLayer::new()
                        .allow_origin(Any)
                        .allow_methods(Any)
                        .allow_headers(Any),
                )
                .layer(ConcurrencyLimitLayer::new(cli.max_concurrent)),
        );

    let addr = format!("{}:{}", cli.bind, cli.port);
    let listener = tokio::net::TcpListener::bind(&addr).await?;
    tracing::info!("fastVEP web server listening on http://{}", addr);

    axum::serve(listener, app)
        .with_graceful_shutdown(shutdown_signal())
        .await?;

    Ok(())
}

async fn shutdown_signal() {
    tokio::signal::ctrl_c()
        .await
        .expect("Failed to listen for ctrl-c");
    tracing::info!("Shutdown signal received, draining connections...");
}
