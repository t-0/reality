/// PDG Monte Carlo particle numbering scheme.
/// Positive values are particles, negative are antiparticles.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PdgId(pub i32);

impl PdgId {
    // Quarks
    pub const DOWN: Self = Self(1);
    pub const UP: Self = Self(2);
    pub const STRANGE: Self = Self(3);
    pub const CHARM: Self = Self(4);
    pub const BOTTOM: Self = Self(5);
    pub const TOP: Self = Self(6);

    // Leptons
    pub const ELECTRON: Self = Self(11);
    pub const ELECTRON_NEUTRINO: Self = Self(12);
    pub const MUON: Self = Self(13);
    pub const MUON_NEUTRINO: Self = Self(14);
    pub const TAU: Self = Self(15);
    pub const TAU_NEUTRINO: Self = Self(16);

    // Gauge bosons
    pub const GLUON: Self = Self(21);
    pub const PHOTON: Self = Self(22);
    pub const Z_BOSON: Self = Self(23);
    pub const W_BOSON: Self = Self(24);

    // Higgs
    pub const HIGGS: Self = Self(25);

    // Antiparticles (for convenience)
    pub const ANTI_DOWN: Self = Self(-1);
    pub const ANTI_UP: Self = Self(-2);
    pub const ANTI_STRANGE: Self = Self(-3);
    pub const ANTI_CHARM: Self = Self(-4);
    pub const ANTI_BOTTOM: Self = Self(-5);
    pub const ANTI_TOP: Self = Self(-6);

    pub const POSITRON: Self = Self(-11);
    pub const ANTI_ELECTRON_NEUTRINO: Self = Self(-12);
    pub const ANTI_MUON: Self = Self(-13);
    pub const ANTI_MUON_NEUTRINO: Self = Self(-14);
    pub const ANTI_TAU: Self = Self(-15);
    pub const ANTI_TAU_NEUTRINO: Self = Self(-16);

    pub const W_MINUS: Self = Self(-24);

    pub fn antiparticle(self) -> Self {
        if self.is_self_conjugate() {
            self
        } else {
            Self(-self.0)
        }
    }

    pub fn is_self_conjugate(self) -> bool {
        matches!(
            self,
            Self::PHOTON | Self::Z_BOSON | Self::GLUON | Self::HIGGS
        )
    }

    pub fn is_quark(self) -> bool {
        let id = self.0.abs();
        (1..=6).contains(&id)
    }

    pub fn is_lepton(self) -> bool {
        let id = self.0.abs();
        (11..=16).contains(&id)
    }

    pub fn is_neutrino(self) -> bool {
        let id = self.0.abs();
        matches!(id, 12 | 14 | 16)
    }

    pub fn is_boson(self) -> bool {
        let id = self.0.abs();
        (21..=25).contains(&id)
    }

    pub fn is_anti(self) -> bool {
        self.0 < 0
    }

    /// Generation (1, 2, 3) for quarks and leptons, None for bosons.
    pub fn generation(self) -> Option<u8> {
        let id = self.0.abs();
        match id {
            1 | 2 | 11 | 12 => Some(1),
            3 | 4 | 13 | 14 => Some(2),
            5 | 6 | 15 | 16 => Some(3),
            _ => None,
        }
    }
}

impl std::fmt::Display for PdgId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_antiparticle() {
        assert_eq!(PdgId::ELECTRON.antiparticle(), PdgId::POSITRON);
        assert_eq!(PdgId::POSITRON.antiparticle(), PdgId::ELECTRON);
        assert_eq!(PdgId::PHOTON.antiparticle(), PdgId::PHOTON);
    }

    #[test]
    fn test_classification() {
        assert!(PdgId::UP.is_quark());
        assert!(PdgId::ANTI_TOP.is_quark());
        assert!(PdgId::ELECTRON.is_lepton());
        assert!(PdgId::MUON_NEUTRINO.is_neutrino());
        assert!(PdgId::PHOTON.is_boson());
        assert!(PdgId::HIGGS.is_boson());
    }

    #[test]
    fn test_generation() {
        assert_eq!(PdgId::UP.generation(), Some(1));
        assert_eq!(PdgId::CHARM.generation(), Some(2));
        assert_eq!(PdgId::TOP.generation(), Some(3));
        assert_eq!(PdgId::PHOTON.generation(), None);
    }
}
