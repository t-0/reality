use crate::particles::PdgId;

/// Trait for computing spin-averaged/summed squared matrix elements |M|^2.
pub trait MatrixElement {
    /// Compute the spin-averaged/summed |M|^2 for given Mandelstam variables.
    /// For decay processes, s = m_parent^2 and t, u are decay kinematics.
    fn amplitude_squared(&self, s: f64, t: f64) -> f64;

    /// Initial state particle PDG IDs.
    fn initial_state(&self) -> Vec<PdgId>;

    /// Final state particle PDG IDs.
    fn final_state(&self) -> Vec<PdgId>;

    /// Number of initial-state spin/color degrees of freedom for averaging.
    /// For e+e-: 2*2 = 4 (spins). For qq: 2*2*3*3 = 36 (spin x color).
    fn initial_state_dof(&self) -> f64;

    /// Symmetry factor for identical final-state particles. Default = 1.
    fn symmetry_factor(&self) -> f64 {
        1.0
    }

    /// Final state masses.
    fn final_masses(&self) -> (f64, f64);

    /// Description of this process.
    fn name(&self) -> &str;
}
