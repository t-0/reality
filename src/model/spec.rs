use std::path::Path;

use serde::Deserialize;

use super::field::{parse_rational, Field, Kind, Su3};
use crate::couplings::constants::M_Z;
use crate::{RealityError, Result};

/// Names of the three gauge couplings, in order: U(1)_Y, SU(2)_L, SU(3)_c.
pub const COUPLING_NAMES: [&str; 3] = ["alpha_1", "alpha_2", "alpha_3"];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Base {
    /// Start from the Standard Model field content.
    #[default]
    Sm,
    /// Start from gauge bosons only; list every field explicitly.
    None,
}

#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum Number {
    Text(String),
    Value(f64),
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct FieldSpec {
    name: String,
    kind: Kind,
    #[serde(default = "default_su3")]
    su3: String,
    #[serde(default = "default_su2")]
    su2: u32,
    #[serde(default = "default_y")]
    y: Number,
    #[serde(default = "default_copies")]
    copies: u32,
    #[serde(default = "default_mass")]
    mass: f64,
}

fn default_su3() -> String {
    "1".to_string()
}
fn default_su2() -> u32 {
    1
}
fn default_y() -> Number {
    Number::Value(0.0)
}
fn default_copies() -> u32 {
    1
}
fn default_mass() -> f64 {
    M_Z
}
fn default_inputs() -> Vec<String> {
    vec!["alpha_1".to_string(), "alpha_2".to_string()]
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct SpecFile {
    name: String,
    #[serde(default)]
    description: Option<String>,
    #[serde(default)]
    base: Base,
    #[serde(default = "default_inputs")]
    inputs: Vec<String>,
    #[serde(default)]
    field: Vec<FieldSpec>,
}

/// A model amendment: new fields on top of a base, plus the declared inputs.
/// The coupling not listed in `inputs` is held out and predicted.
#[derive(Debug, Clone)]
pub struct Amendment {
    pub name: String,
    pub description: Option<String>,
    pub base: Base,
    /// Indices of the two couplings taken from measurement.
    pub inputs: (usize, usize),
    /// Index of the held-out coupling.
    pub predicted: usize,
    /// Fields added on top of the base.
    pub added: Vec<Field>,
}

impl Amendment {
    pub fn from_file(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|e| {
            RealityError::InvalidModel(format!("cannot read {}: {}", path.display(), e))
        })?;
        Self::from_toml(&text)
    }

    pub fn from_toml(text: &str) -> Result<Self> {
        let spec: SpecFile =
            toml::from_str(text).map_err(|e| RealityError::InvalidModel(e.to_string()))?;

        let idx = |name: &str| {
            COUPLING_NAMES.iter().position(|n| *n == name).ok_or_else(|| {
                RealityError::InvalidModel(format!(
                    "unknown input '{}' (use alpha_1, alpha_2, alpha_3)",
                    name
                ))
            })
        };
        if spec.inputs.len() != 2 {
            return Err(RealityError::InvalidModel(
                "inputs must name exactly two couplings; the third is predicted".to_string(),
            ));
        }
        let (i, j) = (idx(&spec.inputs[0])?, idx(&spec.inputs[1])?);
        if i == j {
            return Err(RealityError::InvalidModel(
                "the two inputs must be different couplings".to_string(),
            ));
        }
        let predicted = 3 - i - j;

        let mut added = Vec::new();
        for f in spec.field {
            let y = match &f.y {
                Number::Text(s) => parse_rational(s)?,
                Number::Value(v) => *v,
            };
            if f.su2 == 0 {
                return Err(RealityError::InvalidModel(format!(
                    "field '{}': su2 is the dimension of the representation and must be >= 1",
                    f.name
                )));
            }
            if f.mass <= 0.0 {
                return Err(RealityError::InvalidModel(format!(
                    "field '{}': mass must be positive",
                    f.name
                )));
            }
            added.push(Field {
                name: f.name,
                kind: f.kind,
                su3: Su3::parse(&f.su3)?,
                su2: f.su2,
                y,
                copies: f.copies,
                mass: f.mass,
            });
        }

        Ok(Self {
            name: spec.name,
            description: spec.description,
            base: spec.base,
            inputs: (i.min(j), i.max(j)),
            predicted,
            added,
        })
    }

    /// Base fields followed by the added fields.
    pub fn all_fields(&self) -> Vec<Field> {
        let mut fields = match self.base {
            Base::Sm => super::sm::standard_model(),
            Base::None => Vec::new(),
        };
        fields.extend(self.added.iter().cloned());
        fields
    }

    /// Distinct mass scales above M_Z introduced by the added fields.
    /// Each one is a parameter set by hand.
    pub fn mass_scales(&self) -> Vec<f64> {
        let mut scales: Vec<f64> = self
            .added
            .iter()
            .map(|f| f.mass)
            .filter(|&m| m > M_Z)
            .collect();
        scales.sort_by(|a, b| a.partial_cmp(b).unwrap());
        scales.dedup_by(|a, b| (*a / *b - 1.0).abs() < 1e-9);
        scales
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_minimal() {
        let a = Amendment::from_toml("name = \"SM\"").unwrap();
        assert_eq!(a.inputs, (0, 1));
        assert_eq!(a.predicted, 2);
        assert!(a.added.is_empty());
    }

    #[test]
    fn test_parse_field() {
        let a = Amendment::from_toml(
            r#"
            name = "test"
            inputs = ["alpha_3", "alpha_1"]
            [[field]]
            name = "vlq"
            kind = "dirac"
            su3 = "3"
            su2 = 2
            y = "1/6"
            mass = 2000
            "#,
        )
        .unwrap();
        assert_eq!(a.inputs, (0, 2));
        assert_eq!(a.predicted, 1);
        assert_eq!(a.added[0].su3, Su3::Triplet);
        assert!((a.added[0].y - 1.0 / 6.0).abs() < 1e-12);
        assert_eq!(a.mass_scales(), vec![2000.0]);
    }

    #[test]
    fn test_rejects_bad_input() {
        assert!(Amendment::from_toml("name = \"x\"\ninputs = [\"alpha_1\"]").is_err());
        assert!(Amendment::from_toml("name = \"x\"\ninputs = [\"alpha_1\", \"alpha_1\"]").is_err());
        assert!(Amendment::from_toml("name = \"x\"\ninputs = [\"alpha_1\", \"g\"]").is_err());
        assert!(Amendment::from_toml("name = \"x\"\nbogus = 1").is_err());
    }
}
