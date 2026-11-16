use serde::Deserialize;

use crate::couplings::constants::M_Z;
use crate::{RealityError, Result};

/// How a field contributes to the beta functions and anomalies.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Kind {
    /// One left-handed Weyl fermion in the given representation (chiral).
    Weyl,
    /// A vectorlike pair: Weyl fermions in R and its conjugate (anomaly-free).
    Dirac,
    /// A complex scalar.
    ComplexScalar,
    /// A real scalar (only sensible in real representations with Y = 0).
    RealScalar,
}

impl Kind {
    /// Weight multiplying the Dynkin index in the one-loop beta coefficient.
    pub fn beta_weight(self) -> f64 {
        match self {
            Kind::Weyl => 2.0 / 3.0,
            Kind::Dirac => 4.0 / 3.0,
            Kind::ComplexScalar => 1.0 / 3.0,
            Kind::RealScalar => 1.0 / 6.0,
        }
    }
}

/// SU(3) colour representations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Su3 {
    Singlet,
    Triplet,
    AntiTriplet,
    Sextet,
    AntiSextet,
    Octet,
}

impl Su3 {
    pub fn parse(s: &str) -> Result<Self> {
        match s.trim() {
            "1" => Ok(Su3::Singlet),
            "3" => Ok(Su3::Triplet),
            "3bar" | "3b" | "3̄" => Ok(Su3::AntiTriplet),
            "6" => Ok(Su3::Sextet),
            "6bar" | "6b" => Ok(Su3::AntiSextet),
            "8" => Ok(Su3::Octet),
            other => Err(RealityError::InvalidModel(format!(
                "unknown SU(3) representation '{}' (use 1, 3, 3bar, 6, 6bar, 8)",
                other
            ))),
        }
    }

    pub fn dim(self) -> f64 {
        match self {
            Su3::Singlet => 1.0,
            Su3::Triplet | Su3::AntiTriplet => 3.0,
            Su3::Sextet | Su3::AntiSextet => 6.0,
            Su3::Octet => 8.0,
        }
    }

    /// Dynkin index T(R), normalized so T(fundamental) = 1/2.
    pub fn index(self) -> f64 {
        match self {
            Su3::Singlet => 0.0,
            Su3::Triplet | Su3::AntiTriplet => 0.5,
            Su3::Sextet | Su3::AntiSextet => 2.5,
            Su3::Octet => 3.0,
        }
    }

    /// Cubic anomaly coefficient A(R), normalized so A(3) = 1.
    pub fn cubic_anomaly(self) -> f64 {
        match self {
            Su3::Singlet | Su3::Octet => 0.0,
            Su3::Triplet => 1.0,
            Su3::AntiTriplet => -1.0,
            Su3::Sextet => 7.0,
            Su3::AntiSextet => -7.0,
        }
    }

    pub fn conjugate(self) -> Self {
        match self {
            Su3::Triplet => Su3::AntiTriplet,
            Su3::AntiTriplet => Su3::Triplet,
            Su3::Sextet => Su3::AntiSextet,
            Su3::AntiSextet => Su3::Sextet,
            other => other,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Su3::Singlet => "1",
            Su3::Triplet => "3",
            Su3::AntiTriplet => "3bar",
            Su3::Sextet => "6",
            Su3::AntiSextet => "6bar",
            Su3::Octet => "8",
        }
    }
}

/// SU(2) representation of dimension n = 2j + 1.
pub fn su2_index(n: u32) -> f64 {
    let n = n as f64;
    (n - 1.0) * n * (n + 1.0) / 12.0
}

/// Parse a hypercharge written as "1/6", "-2/3", "0" or a decimal.
pub fn parse_rational(s: &str) -> Result<f64> {
    let s = s.trim();
    let bad = || RealityError::InvalidModel(format!("cannot read hypercharge '{}'", s));
    match s.split_once('/') {
        Some((n, d)) => {
            let n: f64 = n.trim().parse().map_err(|_| bad())?;
            let d: f64 = d.trim().parse().map_err(|_| bad())?;
            if d == 0.0 {
                return Err(bad());
            }
            Ok(n / d)
        }
        None => s.parse().map_err(|_| bad()),
    }
}

/// A gauge multiplet with its quantum numbers and mass threshold.
/// Hypercharge uses the convention Q = T3 + Y.
#[derive(Debug, Clone)]
pub struct Field {
    pub name: String,
    pub kind: Kind,
    pub su3: Su3,
    pub su2: u32,
    pub y: f64,
    pub copies: u32,
    /// Mass in GeV; the field contributes to running above this scale.
    pub mass: f64,
}

impl Field {
    pub fn new(name: &str, kind: Kind, su3: Su3, su2: u32, y: f64, copies: u32) -> Self {
        Self {
            name: name.to_string(),
            kind,
            su3,
            su2,
            y,
            copies,
            mass: M_Z,
        }
    }

    pub fn with_mass(mut self, mass: f64) -> Self {
        self.mass = mass;
        self
    }

    /// Contribution (b1, b2, b3) to the one-loop beta coefficients,
    /// with U(1)_Y in GUT normalization (factor 3/5).
    pub fn beta_contribution(&self) -> [f64; 3] {
        let w = self.kind.beta_weight() * self.copies as f64;
        let d3 = self.su3.dim();
        let d2 = self.su2 as f64;
        [
            w * 0.6 * self.y * self.y * d3 * d2,
            w * su2_index(self.su2) * d3,
            w * self.su3.index() * d2,
        ]
    }

    /// The field as a list of left-handed Weyl multiplets, for anomaly counting.
    /// Scalars return nothing; a Dirac field returns R and its conjugate.
    pub fn weyl_components(&self) -> Vec<(Su3, u32, f64, u32)> {
        match self.kind {
            Kind::Weyl => vec![(self.su3, self.su2, self.y, self.copies)],
            Kind::Dirac => vec![
                (self.su3, self.su2, self.y, self.copies),
                (self.su3.conjugate(), self.su2, -self.y, self.copies),
            ],
            Kind::ComplexScalar | Kind::RealScalar => vec![],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_rational() {
        assert!((parse_rational("1/6").unwrap() - 1.0 / 6.0).abs() < 1e-12);
        assert!((parse_rational("-2/3").unwrap() + 2.0 / 3.0).abs() < 1e-12);
        assert_eq!(parse_rational("0").unwrap(), 0.0);
        assert!(parse_rational("1/0").is_err());
        assert!(parse_rational("x").is_err());
    }

    #[test]
    fn test_su2_index() {
        assert_eq!(su2_index(1), 0.0);
        assert_eq!(su2_index(2), 0.5);
        assert_eq!(su2_index(3), 2.0);
        assert_eq!(su2_index(4), 5.0);
    }
}
