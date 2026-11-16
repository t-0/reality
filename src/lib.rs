pub mod couplings;
pub mod cross_section;
pub mod decay;
pub mod dirac;
pub mod kinematics;
pub mod matrix_elements;
pub mod model;
pub mod monte_carlo;
pub mod particles;
pub mod phase_space;
pub mod vertices;

pub mod cli;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum RealityError {
    #[error("Unknown particle with PDG ID {0}")]
    UnknownParticle(i32),

    #[error("Process '{0}' is not implemented")]
    UnknownProcess(String),

    #[error("Kinematically forbidden: sqrt(s) = {sqrt_s} GeV is below threshold {threshold} GeV")]
    BelowThreshold { sqrt_s: f64, threshold: f64 },

    #[error("Phase space generation failed after {0} attempts")]
    PhaseSpaceFailure(usize),

    #[error("Integration did not converge: chi^2/dof = {0:.2}")]
    IntegrationNotConverged(f64),

    #[error("Numerical error: {0}")]
    Numerical(String),

    #[error("Invalid model: {0}")]
    InvalidModel(String),
}

pub type Result<T> = std::result::Result<T, RealityError>;
