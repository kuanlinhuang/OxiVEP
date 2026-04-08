use anyhow::Result;
use clap::{Parser, Subcommand};

mod pipeline;
mod webserver;

#[derive(Parser)]
#[command(name = "oxivep")]
#[command(about = "OxiVEP - A high-performance variant effect predictor")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Annotate variants with predicted consequences
    Annotate {
        /// Input file (VCF format). Use "-" for stdin.
        #[arg(short, long)]
        input: String,

        /// Output file. Use "-" for stdout.
        #[arg(short, long, default_value = "-")]
        output: String,

        /// GFF3 annotation file for transcript models
        #[arg(long)]
        gff3: Option<String>,

        /// Path to FASTA reference file
        #[arg(long)]
        fasta: Option<String>,

        /// Output format (vcf, tab, json)
        #[arg(long, default_value = "vcf")]
        output_format: String,

        /// Turn on all common annotation flags
        #[arg(long)]
        everything: bool,

        /// Number of variants to buffer
        #[arg(long, default_value_t = 5000)]
        buffer_size: usize,

        /// Pick one consequence per variant (most severe)
        #[arg(long)]
        pick: bool,

        /// Include gene symbol in output
        #[arg(long)]
        symbol: bool,

        /// Include HGVS notations
        #[arg(long)]
        hgvs: bool,

        /// Include canonical transcript flag
        #[arg(long)]
        canonical: bool,

        /// Upstream/downstream distance (bp)
        #[arg(long, default_value_t = 5000)]
        distance: u64,

        /// Path to VEP cache directory for known variant annotation
        #[arg(long)]
        cache_dir: Option<String>,

        /// Path to binary transcript cache file (auto-generated if not specified)
        #[arg(long)]
        transcript_cache: Option<String>,

        /// Directory containing supplementary annotation files (.osa, .osi, .oga)
        #[arg(long)]
        sa_dir: Option<String>,
    },

    /// Launch the web interface for interactive variant annotation
    Web {
        /// Port to listen on
        #[arg(long, default_value_t = 8080)]
        port: u16,

        /// GFF3 annotation file
        #[arg(long)]
        gff3: Option<String>,

        /// Path to FASTA reference file
        #[arg(long)]
        fasta: Option<String>,
    },

    /// Build a binary transcript cache for fast startup
    Cache {
        /// GFF3 annotation file
        #[arg(long)]
        gff3: String,

        /// Path to FASTA reference file (for pre-building sequences)
        #[arg(long)]
        fasta: Option<String>,

        /// Output cache file path
        #[arg(short, long)]
        output: String,
    },

    /// Build a supplementary annotation database (.osa) from a source file
    SaBuild {
        /// Source type: clinvar, gnomad, dbsnp
        #[arg(long)]
        source: String,

        /// Input file (VCF, typically gzipped)
        #[arg(short, long)]
        input: String,

        /// Output base path (will create .osa and .osa.idx)
        #[arg(short, long)]
        output: String,

        /// Genome assembly (e.g., GRCh38)
        #[arg(long, default_value = "GRCh38")]
        assembly: String,
    },

    /// Filter annotated VEP output
    Filter {
        /// Input file (VEP-annotated VCF)
        #[arg(short, long)]
        input: String,

        /// Output file
        #[arg(short, long, default_value = "-")]
        output: String,

        /// Filter expression
        #[arg(long)]
        filter: String,
    },
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Annotate {
            input,
            output,
            gff3,
            fasta,
            output_format,
            everything: _,
            buffer_size: _,
            pick,
            symbol: _,
            hgvs,
            canonical: _,
            distance,
            cache_dir,
            transcript_cache,
            sa_dir,
        } => {
            pipeline::run_annotate(pipeline::AnnotateConfig {
                input,
                output,
                gff3,
                fasta,
                output_format,
                pick,
                hgvs,
                distance,
                cache_dir,
                transcript_cache,
                sa_dir,
            })?;
        }
        Commands::Cache { gff3, fasta, output } => {
            pipeline::run_cache_build(&gff3, fasta.as_deref(), &output)?;
        }
        Commands::Web { port, gff3, fasta } => {
            webserver::run_server(port, gff3, fasta)?;
        }
        Commands::SaBuild {
            source,
            input,
            output,
            assembly,
        } => {
            pipeline::run_sa_build(&source, &input, &output, &assembly)?;
        }
        Commands::Filter { input, filter, .. } => {
            log::info!("OxiVEP filter: input={}, filter={}", input, filter);
            eprintln!("OxiVEP filter engine is under construction.");
        }
    }

    Ok(())
}
