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

    /// Scan an amendment over the mass of its new fields and find the best fit
    Scan {
        /// Path to a model file, e.g. models/mssm.toml
        file: std::path::PathBuf,

        /// Lowest mass of the lightest new field (GeV)
        #[arg(long, default_value = "100")]
        from: f64,

        /// Highest mass of the lightest new field (GeV)
        #[arg(long, default_value = "1e7")]
        to: f64,

        /// Number of masses, evenly spaced in log
        #[arg(long, default_value = "31")]
        points: usize,

        /// Loop order of the running: 1 or 2
        #[arg(long, default_value = "2", value_parser = clap::value_parser!(u8).range(1..=2))]
        loops: u8,

        /// Output format: "table" or "csv"
        #[arg(long, default_value = "table")]
        format: String,
    },

    /// Search for non-supersymmetric additions to the SM that unify and survive
    Search {
        /// Largest number of new multiplets per candidate
        #[arg(long, default_value = "3")]
        max_fields: usize,

        /// Comma-separated mass grid in GeV (default 1e3,1e4,1e5,1e6,1e8,1e10,1e12)
        #[arg(long)]
        masses: Option<String>,

        /// Loop order for the final evaluation: 1 or 2
        #[arg(long, default_value = "2", value_parser = clap::value_parser!(u8).range(1..=2))]
        loops: u8,

        /// Accept candidates with |pull| at or below this
        #[arg(long, default_value = "1.0")]
        accept: f64,

        /// Candidates to show for each number of fields
        #[arg(long, default_value = "5")]
        top: usize,

        /// Write the shown candidates as model files into this directory
        #[arg(long)]
        write: Option<std::path::PathBuf>,

        /// Rank accepted models by robustness, using this many random draws each (0 = off)
        #[arg(long, default_value = "0")]
        robust: usize,
    },

    /// Measure how much of the nearby model space still works for a model
    Robust {
        /// Path to a model file
        file: std::path::PathBuf,

        /// Number of random draws
        #[arg(long, default_value = "200")]
        samples: usize,

        /// Each added mass moves log-uniformly within this many decades either way
        #[arg(long, default_value = "1")]
        spread: f64,

        /// Draws with any added mass below this (GeV) fail as excluded by colliders
        #[arg(long, default_value = "1000")]
        min_mass: f64,

        /// Generic SU(5) GUT-field ratios range over [1/R, R] if the model has no [[heavy]] fields (1 = off)
        #[arg(long, default_value = "10")]
        thresholds: f64,

        /// Random seed
        #[arg(long, default_value = "1")]
        seed: u64,
    },
}
