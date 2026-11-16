use super::field::{Field, Kind, Su3};

/// Standard Model field content as left-handed Weyl multiplets plus the Higgs.
/// Right-handed fields appear as their conjugates (u^c, d^c, e^c).
/// All fields are treated as active from M_Z, as in the usual one-loop analysis.
pub fn standard_model() -> Vec<Field> {
    vec![
        Field::new("Q", Kind::Weyl, Su3::Triplet, 2, 1.0 / 6.0, 3),
        Field::new("u^c", Kind::Weyl, Su3::AntiTriplet, 1, -2.0 / 3.0, 3),
        Field::new("d^c", Kind::Weyl, Su3::AntiTriplet, 1, 1.0 / 3.0, 3),
        Field::new("L", Kind::Weyl, Su3::Singlet, 2, -0.5, 3),
        Field::new("e^c", Kind::Weyl, Su3::Singlet, 1, 1.0, 3),
        Field::new("H", Kind::ComplexScalar, Su3::Singlet, 2, 0.5, 1),
    ]
}
