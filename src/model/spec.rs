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

/// Which top Yukawa running to include at two loops.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum TopYukawa {
    /// MSSM coefficients in supersymmetric energy ranges, SM coefficients elsewhere.
    #[default]
    Auto,
    /// Leave the Yukawa contribution out.
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
#[serde(untagged)]
enum RatioSpec {
    One(f64),
    Range([f64; 2]),
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct HeavySpec {
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
    /// Mass relative to the X boson mass M_X, either fixed or a range [low, high].
    ratio: RatioSpec,
    /// Heavy fields with the same group share one mass ratio.
    #[serde(default)]
    group: Option<String>,
}

/// A field at the unification scale, entering through threshold matching.
#[derive(Debug, Clone)]
pub struct HeavyField {
    pub field: Field,
    /// Range of M / M_X; equal ends mean a fixed ratio.
    pub ratio: (f64, f64),
    /// Fields in the same group share one ratio.
    pub group: Option<String>,
}

impl HeavyField {
    /// Geometric middle of the ratio range.
    pub fn central_ratio(&self) -> f64 {
        (self.ratio.0 * self.ratio.1).sqrt()
    }

    pub fn is_varied(&self) -> bool {
        self.ratio.1 > self.ratio.0 * (1.0 + 1e-12)
    }
}

/// Largest number of varied heavy fields; corners grow as 2^n.
pub const MAX_VARIED_HEAVY: usize = 8;

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
    top_yukawa: TopYukawa,
    /// Only used for dimension-5 proton decay in supersymmetric models.
    #[serde(default)]
    tan_beta: Option<f64>,
    #[serde(default)]
    field: Vec<FieldSpec>,
    #[serde(default)]
    heavy: Vec<HeavySpec>,
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
    /// Top Yukawa scheme for two-loop running.
    pub top_yukawa: TopYukawa,
    /// tan(beta) for dimension-5 proton decay; None means the default.
    pub tan_beta: Option<f64>,
    /// Fields added on top of the base.
    pub added: Vec<Field>,
    /// Fields near the unification scale, matched at M_X instead of run.
    pub heavy: Vec<HeavyField>,
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

        let read_y = |y: &Number| match y {
            Number::Text(s) => parse_rational(s),
            Number::Value(v) => Ok(*v),
        };

        let mut added = Vec::new();
        for f in spec.field {
            if f.kind.is_heavy_only() {
                return Err(RealityError::InvalidModel(format!(
                    "field '{}': kind {:?} is only allowed in [[heavy]]",
                    f.name, f.kind
                )));
            }
            let y = read_y(&f.y)?;
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

        let mut heavy = Vec::new();
        for h in spec.heavy {
            let ratio = match h.ratio {
                RatioSpec::One(r) => (r, r),
                RatioSpec::Range([a, b]) => (a, b),
            };
            if !(ratio.0 > 0.0 && ratio.1 >= ratio.0) {
                return Err(RealityError::InvalidModel(format!(
                    "heavy field '{}': ratio must be positive, as r or [low, high] with low <= high",
                    h.name
                )));
            }
            if h.su2 == 0 {
                return Err(RealityError::InvalidModel(format!(
                    "heavy field '{}': su2 must be >= 1",
                    h.name
                )));
            }
            heavy.push(HeavyField {
                field: Field {
                    name: h.name,
                    kind: h.kind,
                    su3: Su3::parse(&h.su3)?,
                    su2: h.su2,
                    y: read_y(&h.y)?,
                    copies: h.copies,
                    mass: 0.0,
                },
                ratio,
                group: h.group,
            });
        }
        for a in &heavy {
            for b in &heavy {
                if a.group.is_some() && a.group == b.group && a.ratio != b.ratio {
                    return Err(RealityError::InvalidModel(format!(
                        "heavy fields '{}' and '{}' share group '{}' but have different ratios",
                        a.field.name,
                        b.field.name,
                        a.group.as_deref().unwrap_or("")
                    )));
                }
            }
        }
        if spec.tan_beta.map_or(false, |t| !(t > 0.0)) {
            return Err(RealityError::InvalidModel("tan_beta must be positive".to_string()));
        }
        if Self::varied_units_of(&heavy).len() > MAX_VARIED_HEAVY {
            return Err(RealityError::InvalidModel(format!(
                "at most {} independent heavy-field ratios may be varied",
                MAX_VARIED_HEAVY
            )));
        }

        Ok(Self {
            name: spec.name,
            description: spec.description,
            base: spec.base,
            inputs: (i.min(j), i.max(j)),
            predicted,
            top_yukawa: spec.top_yukawa,
            tan_beta: spec.tan_beta,
            added,
            heavy,
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

    /// One-loop threshold matching at M_X: heavy field h with M_h = r_h M_X shifts
    /// each 1/alpha_i(M_X) by b_i^h / (2 pi) * ln(r_h), so unification reads
    /// 1/alpha_i(M_X) + shift_i = 1/alpha_GUT.
    pub fn threshold_shift(&self, ratios: &[f64]) -> [f64; 3] {
        let mut shift = [0.0; 3];
        for (h, &r) in self.heavy.iter().zip(ratios) {
            let b = h.field.beta_contribution();
            for i in 0..3 {
                shift[i] += b[i] / (2.0 * std::f64::consts::PI) * r.ln();
            }
        }
        shift
    }

    pub fn central_ratios(&self) -> Vec<f64> {
        self.heavy.iter().map(|h| h.central_ratio()).collect()
    }

    /// Independently varied ratios: each entry lists the heavy fields that share it.
    fn varied_units_of(heavy: &[HeavyField]) -> Vec<Vec<usize>> {
        let mut units: Vec<Vec<usize>> = Vec::new();
        for (k, h) in heavy.iter().enumerate().filter(|(_, h)| h.is_varied()) {
            let existing = h.group.as_ref().and_then(|g| {
                units.iter().position(|u| heavy[u[0]].group.as_ref() == Some(g))
            });
            match existing {
                Some(u) => units[u].push(k),
                None => units.push(vec![k]),
            }
        }
        units
    }

    /// Number of independent mass ratios varied over a range.
    pub fn varied_ratio_count(&self) -> usize {
        Self::varied_units_of(&self.heavy).len()
    }

    /// Every combination of range ends for the varied ratios (fixed ones keep
    /// their value, grouped fields move together). Empty if nothing is varied.
    pub fn ratio_corners(&self) -> Vec<Vec<f64>> {
        let units = Self::varied_units_of(&self.heavy);
        if units.is_empty() {
            return Vec::new();
        }
        (0..1usize << units.len())
            .map(|mask| {
                let mut r = self.central_ratios();
                for (bit, unit) in units.iter().enumerate() {
                    for &k in unit {
                        let (lo, hi) = self.heavy[k].ratio;
                        r[k] = if mask & (1 << bit) != 0 { hi } else { lo };
                    }
                }
                r
            })
            .collect()
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
    fn test_heavy_fields() {
        let a = Amendment::from_toml(
            r#"
            name = "su5"
            [[heavy]]
            name = "XY"
            kind = "vector_superfield"
            su3 = "3"
            su2 = 2
            y = "-5/6"
            ratio = 1.0
            [[heavy]]
            name = "H_C"
            kind = "chiral"
            su3 = "3"
            y = "-1/3"
            copies = 2
            ratio = [0.1, 10]
            "#,
        )
        .unwrap();
        // X, Y vector superfield reproduces the SU(5) values (-10, -6, -4).
        let b = a.heavy[0].field.beta_contribution();
        assert!((b[0] + 10.0).abs() < 1e-12 && (b[1] + 6.0).abs() < 1e-12 && (b[2] + 4.0).abs() < 1e-12);
        assert_eq!(a.ratio_corners().len(), 2);
        assert!((a.central_ratios()[1] - 1.0).abs() < 1e-12);
        // Ratio 1 everywhere gives no shift.
        assert_eq!(a.threshold_shift(&[1.0, 1.0]), [0.0; 3]);
    }

    #[test]
    fn test_grouped_heavy_fields_move_together() {
        let a = Amendment::from_toml(
            r#"
            name = "g"
            [[heavy]]
            name = "s8"
            kind = "chiral"
            su3 = "8"
            ratio = [0.1, 10]
            group = "sigma"
            [[heavy]]
            name = "s3"
            kind = "chiral"
            su2 = 3
            ratio = [0.1, 10]
            group = "sigma"
            "#,
        )
        .unwrap();
        let corners = a.ratio_corners();
        assert_eq!(corners.len(), 2);
        assert!(corners.iter().all(|r| r[0] == r[1]));
        assert!(Amendment::from_toml(
            "name = \"g\"\n[[heavy]]\nname = \"a\"\nkind = \"chiral\"\nratio = 1.0\ngroup = \"x\"\n[[heavy]]\nname = \"b\"\nkind = \"chiral\"\nratio = 2.0\ngroup = \"x\""
        )
        .is_err());
    }

    #[test]
    fn test_heavy_kind_rejected_in_light_fields() {
        assert!(Amendment::from_toml("name = \"x\"\n[[field]]\nname = \"v\"\nkind = \"vector\"").is_err());
    }

    #[test]
    fn test_rejects_bad_input() {
        assert!(Amendment::from_toml("name = \"x\"\ninputs = [\"alpha_1\"]").is_err());
        assert!(Amendment::from_toml("name = \"x\"\ninputs = [\"alpha_1\", \"alpha_1\"]").is_err());
        assert!(Amendment::from_toml("name = \"x\"\ninputs = [\"alpha_1\", \"g\"]").is_err());
        assert!(Amendment::from_toml("name = \"x\"\nbogus = 1").is_err());
    }
}
