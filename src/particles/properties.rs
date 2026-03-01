use super::PdgId;

/// Spin of a particle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Spin {
    Zero,
    Half,
    One,
}

impl Spin {
    /// Number of spin states (2s+1), accounting for massless spin-1 having 2 polarizations.
    pub fn degrees_of_freedom(self, is_massless: bool) -> u32 {
        match (self, is_massless) {
            (Spin::Zero, _) => 1,
            (Spin::Half, _) => 2,
            (Spin::One, true) => 2,
            (Spin::One, false) => 3,
        }
    }
}

impl std::fmt::Display for Spin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Spin::Zero => write!(f, "0"),
            Spin::Half => write!(f, "1/2"),
            Spin::One => write!(f, "1"),
        }
    }
}

/// Color representation of a particle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ColorRep {
    Singlet,
    Triplet,
    AntiTriplet,
    Octet,
}

impl ColorRep {
    /// Number of color degrees of freedom.
    pub fn dimension(self) -> u32 {
        match self {
            ColorRep::Singlet => 1,
            ColorRep::Triplet => 3,
            ColorRep::AntiTriplet => 3,
            ColorRep::Octet => 8,
        }
    }

    /// Color representation of the antiparticle.
    pub fn conjugate(self) -> Self {
        match self {
            ColorRep::Triplet => ColorRep::AntiTriplet,
            ColorRep::AntiTriplet => ColorRep::Triplet,
            other => other,
        }
    }
}

impl std::fmt::Display for ColorRep {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ColorRep::Singlet => write!(f, "1"),
            ColorRep::Triplet => write!(f, "3"),
            ColorRep::AntiTriplet => write!(f, "3*"),
            ColorRep::Octet => write!(f, "8"),
        }
    }
}

/// Physical properties of a Standard Model particle.
#[derive(Debug, Clone)]
pub struct ParticleData {
    pub pdg_id: PdgId,
    pub name: &'static str,
    pub anti_name: Option<&'static str>,
    pub mass: f64,
    pub width: Option<f64>,
    pub charge: f64,
    pub spin: Spin,
    pub color: ColorRep,
    pub weak_isospin_3: f64,
    pub hypercharge: f64,
}

impl ParticleData {
    /// Create the antiparticle data from this particle.
    pub fn conjugate(&self) -> Self {
        if self.pdg_id.is_self_conjugate() {
            return self.clone();
        }
        Self {
            pdg_id: self.pdg_id.antiparticle(),
            name: self.anti_name.unwrap_or(self.name),
            anti_name: Some(self.name),
            mass: self.mass,
            width: self.width,
            charge: -self.charge,
            spin: self.spin,
            color: self.color.conjugate(),
            weak_isospin_3: -self.weak_isospin_3,
            hypercharge: -self.hypercharge,
        }
    }

    pub fn is_stable(&self) -> bool {
        self.width.is_none()
    }

    pub fn is_massless(&self) -> bool {
        self.mass == 0.0
    }
}

impl std::fmt::Display for ParticleData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} (PDG: {}, mass: {:.6} GeV)", self.name, self.pdg_id, self.mass)
    }
}
