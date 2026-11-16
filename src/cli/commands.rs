use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "reality",
    about = "Standard Model particle physics simulator",
    version
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// List Standard Model particles and their properties
    Particles {
        /// Filter: "quarks", "leptons", "bosons", or "all"
        #[arg(default_value = "all")]
        filter: String,
    },

    /// Compute cross section for a process at given center-of-mass energy
    #[command(name = "cross-section")]
    CrossSection {
        /// Process: "ee_mumu", "ee_mumu_ew", "ee_qqbar"
        process: String,

        /// Center-of-mass energy in GeV
        #[arg(long)]
        sqrt_s: f64,

        /// Number of integration points
        #[arg(long, default_value = "64")]
        n_points: usize,
    },

    /// Compute decay table for a particle
    Decay {
        /// Particle: "w", "z", "higgs", "top"
        particle: String,
    },

    /// Generate Monte Carlo events
    Generate {
        /// Process: "ee_mumu", "ee_mumu_ew"
        process: String,

        /// Center-of-mass energy in GeV
        #[arg(long)]
        sqrt_s: f64,

        /// Number of events to generate
        #[arg(long, default_value = "100")]
        n_events: usize,

        /// Output format: "json" or "table"
        #[arg(long, default_value = "json")]
        format: String,

        /// Random seed (optional)
        #[arg(long)]
        seed: Option<u64>,
    },

    /// Show the value of a coupling constant at a given energy scale
    Coupling {
        /// Coupling: "alpha_em", "alpha_s", "sin2w", "gf"
        name: String,

        /// Energy scale Q in GeV (for running couplings)
        #[arg(long)]
        scale: Option<f64>,
    },

    /// Run the three gauge couplings to high energy and test grand unification
    Unification {
        /// Superpartner mass scale in GeV; above it the MSSM beta functions apply
        #[arg(long)]
        susy_scale: Option<f64>,
    },

    /// Evaluate a model amendment file: consistency gates and held-out predictions
    Amend {
        /// Path to a model file, e.g. models/mssm.toml
        file: std::path::PathBuf,

        /// Loop order of the running: 1 or 2
        #[arg(long, default_value = "2", value_parser = clap::value_parser!(u8).range(1..=2))]
        loops: u8,
    },
}
