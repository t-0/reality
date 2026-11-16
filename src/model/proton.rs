//! Proton lifetime from gauge-boson exchange at the unification scale.
//!
//! The superheavy X and Y gauge bosons of SU(5)/SO(10) turn quarks into leptons,
//! giving p -> e+ pi0 through a dimension-6 operator with rate
//!   Gamma ~ alpha_GUT^2 m_p^5 / M_X^4.
//! This is the naive dimensional estimate. Hadronic matrix elements, the
//! difference between M_X and the crossing scale, and threshold corrections
//! make it uncertain by about a factor of 10 either way.
//! Supersymmetric models also have dimension-5 decays (p -> K+ nu-bar) that
//! depend on the superpartner spectrum and are not estimated here.

use crate::couplings::constants::M_PROTON;

/// hbar in GeV * s.
const HBAR_GEV_S: f64 = 6.582_119_569e-25;

/// Seconds per year.
const SECONDS_PER_YEAR: f64 = 3.155_76e7;

/// Factor by which the naive estimate is uncertain in each direction.
pub const UNCERTAINTY_FACTOR: f64 = 10.0;

/// Super-Kamiokande lower limit on tau(p -> e+ pi0), 90% CL, in years
/// (Phys. Rev. D 102, 112011, 2020).
pub const SUPER_K_LIMIT_YEARS: f64 = 2.4e34;

/// Approximate Hyper-Kamiokande reach for p -> e+ pi0 after about 10 years
/// of running (90% CL), in years. Hyper-K is due to start around 2027-28.
pub const HYPER_K_REACH_YEARS: f64 = 1.0e35;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// Even the longest lifetime in the uncertainty band is below the Super-K limit.
    Excluded,
    /// The uncertainty band straddles the Super-K limit.
    Borderline,
    /// Allowed today, and Hyper-K can see it or rule it out.
    TestableAtHyperK,
    /// Allowed, and beyond Hyper-K's reach.
    OutOfReach,
}

#[derive(Debug, Clone, Copy)]
pub struct ProtonDecay {
    /// Mass of the X/Y gauge bosons, taken as the unification scale (GeV).
    pub m_x: f64,
    pub alpha_gut: f64,
    /// Central estimate of tau(p -> e+ pi0) in years.
    pub lifetime_years: f64,
    pub status: Status,
}

impl ProtonDecay {
    pub fn estimate(m_gut: f64, inverse_alpha_gut: f64) -> Self {
        let alpha_gut = 1.0 / inverse_alpha_gut;
        let width = alpha_gut * alpha_gut * M_PROTON.powi(5) / m_gut.powi(4);
        let lifetime_years = HBAR_GEV_S / width / SECONDS_PER_YEAR;

        let low = lifetime_years / UNCERTAINTY_FACTOR;
        let high = lifetime_years * UNCERTAINTY_FACTOR;
        let status = if high < SUPER_K_LIMIT_YEARS {
            Status::Excluded
        } else if low < SUPER_K_LIMIT_YEARS {
            Status::Borderline
        } else if low < HYPER_K_REACH_YEARS {
            Status::TestableAtHyperK
        } else {
            Status::OutOfReach
        };

        Self {
            m_x: m_gut,
            alpha_gut,
            lifetime_years,
            status,
        }
    }

    pub fn range_years(&self) -> (f64, f64) {
        (
            self.lifetime_years / UNCERTAINTY_FACTOR,
            self.lifetime_years * UNCERTAINTY_FACTOR,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_textbook_value() {
        // M_X = 1e16 GeV, alpha_GUT = 1/25 gives the standard ~1e35 year estimate.
        let p = ProtonDecay::estimate(1e16, 25.0);
        assert!(p.lifetime_years > 1e35 && p.lifetime_years < 3e35, "{:e}", p.lifetime_years);
    }

    #[test]
    fn test_scales_as_fourth_power() {
        let a = ProtonDecay::estimate(1e15, 25.0);
        let b = ProtonDecay::estimate(1e16, 25.0);
        assert!((b.lifetime_years / a.lifetime_years / 1e4 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_low_scale_is_excluded() {
        // The SM alpha_1 = alpha_2 crossing near 1e13 GeV: the fate of minimal SU(5).
        assert_eq!(ProtonDecay::estimate(1e13, 42.0).status, Status::Excluded);
        assert_eq!(ProtonDecay::estimate(1e17, 25.0).status, Status::OutOfReach);
    }
}
