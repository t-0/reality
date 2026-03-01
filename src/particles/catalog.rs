use std::collections::HashMap;
use std::sync::LazyLock;

use crate::couplings::constants::*;

use super::pdg::PdgId;
use super::properties::{ColorRep, ParticleData, Spin};

/// Global catalog of all Standard Model particles, keyed by PDG ID.
pub static CATALOG: LazyLock<HashMap<PdgId, ParticleData>> = LazyLock::new(build_catalog);

fn build_catalog() -> HashMap<PdgId, ParticleData> {
    let particles = vec![
        // === Quarks ===
        ParticleData {
            pdg_id: PdgId::UP,
            name: "u",
            anti_name: Some("u~"),
            mass: M_UP,
            width: None,
            charge: 2.0 / 3.0,
            spin: Spin::Half,
            color: ColorRep::Triplet,
            weak_isospin_3: 0.5,
            hypercharge: 1.0 / 3.0,
        },
        ParticleData {
            pdg_id: PdgId::DOWN,
            name: "d",
            anti_name: Some("d~"),
            mass: M_DOWN,
            width: None,
            charge: -1.0 / 3.0,
            spin: Spin::Half,
            color: ColorRep::Triplet,
            weak_isospin_3: -0.5,
            hypercharge: 1.0 / 3.0,
        },
        ParticleData {
            pdg_id: PdgId::CHARM,
            name: "c",
            anti_name: Some("c~"),
            mass: M_CHARM,
            width: None,
            charge: 2.0 / 3.0,
            spin: Spin::Half,
            color: ColorRep::Triplet,
            weak_isospin_3: 0.5,
            hypercharge: 1.0 / 3.0,
        },
        ParticleData {
            pdg_id: PdgId::STRANGE,
            name: "s",
            anti_name: Some("s~"),
            mass: M_STRANGE,
            width: None,
            charge: -1.0 / 3.0,
            spin: Spin::Half,
            color: ColorRep::Triplet,
            weak_isospin_3: -0.5,
            hypercharge: 1.0 / 3.0,
        },
        ParticleData {
            pdg_id: PdgId::TOP,
            name: "t",
            anti_name: Some("t~"),
            mass: M_TOP,
            width: Some(GAMMA_TOP),
            charge: 2.0 / 3.0,
            spin: Spin::Half,
            color: ColorRep::Triplet,
            weak_isospin_3: 0.5,
            hypercharge: 1.0 / 3.0,
        },
        ParticleData {
            pdg_id: PdgId::BOTTOM,
            name: "b",
            anti_name: Some("b~"),
            mass: M_BOTTOM,
            width: None,
            charge: -1.0 / 3.0,
            spin: Spin::Half,
            color: ColorRep::Triplet,
            weak_isospin_3: -0.5,
            hypercharge: 1.0 / 3.0,
        },
        // === Leptons ===
        ParticleData {
            pdg_id: PdgId::ELECTRON,
            name: "e-",
            anti_name: Some("e+"),
            mass: M_ELECTRON,
            width: None,
            charge: -1.0,
            spin: Spin::Half,
            color: ColorRep::Singlet,
            weak_isospin_3: -0.5,
            hypercharge: -1.0,
        },
        ParticleData {
            pdg_id: PdgId::ELECTRON_NEUTRINO,
            name: "nu_e",
            anti_name: Some("nu_e~"),
            mass: 0.0,
            width: None,
            charge: 0.0,
            spin: Spin::Half,
            color: ColorRep::Singlet,
            weak_isospin_3: 0.5,
            hypercharge: -1.0,
        },
        ParticleData {
            pdg_id: PdgId::MUON,
            name: "mu-",
            anti_name: Some("mu+"),
            mass: M_MUON,
            width: None,
            charge: -1.0,
            spin: Spin::Half,
            color: ColorRep::Singlet,
            weak_isospin_3: -0.5,
            hypercharge: -1.0,
        },
        ParticleData {
            pdg_id: PdgId::MUON_NEUTRINO,
            name: "nu_mu",
            anti_name: Some("nu_mu~"),
            mass: 0.0,
            width: None,
            charge: 0.0,
            spin: Spin::Half,
            color: ColorRep::Singlet,
            weak_isospin_3: 0.5,
            hypercharge: -1.0,
        },
        ParticleData {
            pdg_id: PdgId::TAU,
            name: "tau-",
            anti_name: Some("tau+"),
            mass: M_TAU,
            width: None,
            charge: -1.0,
            spin: Spin::Half,
            color: ColorRep::Singlet,
            weak_isospin_3: -0.5,
            hypercharge: -1.0,
        },
        ParticleData {
            pdg_id: PdgId::TAU_NEUTRINO,
            name: "nu_tau",
            anti_name: Some("nu_tau~"),
            mass: 0.0,
            width: None,
            charge: 0.0,
            spin: Spin::Half,
            color: ColorRep::Singlet,
            weak_isospin_3: 0.5,
            hypercharge: -1.0,
        },
        // === Gauge Bosons ===
        ParticleData {
            pdg_id: PdgId::PHOTON,
            name: "gamma",
            anti_name: None,
            mass: 0.0,
            width: None,
            charge: 0.0,
            spin: Spin::One,
            color: ColorRep::Singlet,
            weak_isospin_3: 0.0,
            hypercharge: 0.0,
        },
        ParticleData {
            pdg_id: PdgId::GLUON,
            name: "g",
            anti_name: None,
            mass: 0.0,
            width: None,
            charge: 0.0,
            spin: Spin::One,
            color: ColorRep::Octet,
            weak_isospin_3: 0.0,
            hypercharge: 0.0,
        },
        ParticleData {
            pdg_id: PdgId::Z_BOSON,
            name: "Z",
            anti_name: None,
            mass: M_Z,
            width: Some(GAMMA_Z),
            charge: 0.0,
            spin: Spin::One,
            color: ColorRep::Singlet,
            weak_isospin_3: 0.0,
            hypercharge: 0.0,
        },
        ParticleData {
            pdg_id: PdgId::W_BOSON,
            name: "W+",
            anti_name: Some("W-"),
            mass: M_W,
            width: Some(GAMMA_W),
            charge: 1.0,
            spin: Spin::One,
            color: ColorRep::Singlet,
            weak_isospin_3: 1.0,
            hypercharge: 0.0,
        },
        // === Higgs ===
        ParticleData {
            pdg_id: PdgId::HIGGS,
            name: "H",
            anti_name: None,
            mass: M_H,
            width: Some(GAMMA_H),
            charge: 0.0,
            spin: Spin::Zero,
            color: ColorRep::Singlet,
            weak_isospin_3: -0.5,
            hypercharge: 1.0,
        },
    ];

    let mut map = HashMap::new();
    for p in &particles {
        map.insert(p.pdg_id, p.clone());
        // Insert antiparticle if not self-conjugate
        if !p.pdg_id.is_self_conjugate() {
            let anti = p.conjugate();
            map.insert(anti.pdg_id, anti);
        }
    }
    map
}

/// Look up a particle by PDG ID.
pub fn get_particle(pdg_id: PdgId) -> Option<&'static ParticleData> {
    CATALOG.get(&pdg_id)
}

/// Iterate over all particles in the catalog.
pub fn all_particles() -> impl Iterator<Item = &'static ParticleData> {
    CATALOG.values()
}

/// Iterate over quarks (particles only, not antiparticles).
pub fn quarks() -> impl Iterator<Item = &'static ParticleData> {
    CATALOG
        .values()
        .filter(|p| p.pdg_id.is_quark() && !p.pdg_id.is_anti())
}

/// Iterate over leptons (particles only).
pub fn leptons() -> impl Iterator<Item = &'static ParticleData> {
    CATALOG
        .values()
        .filter(|p| p.pdg_id.is_lepton() && !p.pdg_id.is_anti())
}

/// Iterate over bosons.
pub fn bosons() -> impl Iterator<Item = &'static ParticleData> {
    CATALOG.values().filter(|p| p.pdg_id.is_boson() && !p.pdg_id.is_anti())
}

/// Look up a particle by name (case-insensitive).
pub fn find_by_name(name: &str) -> Option<&'static ParticleData> {
    let lower = name.to_lowercase();
    CATALOG.values().find(|p| {
        p.name.to_lowercase() == lower
            || p.anti_name.map_or(false, |n| n.to_lowercase() == lower)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_catalog_size() {
        // 17 particles + antiparticles for non-self-conjugate ones
        // Self-conjugate: gamma, Z, gluon, Higgs = 4
        // Non-self-conjugate: 6 quarks + 6 leptons + W = 13, so 13 antiparticles
        // Total: 17 + 13 = 30
        assert_eq!(CATALOG.len(), 30);
    }

    #[test]
    fn test_lookup_electron() {
        let e = get_particle(PdgId::ELECTRON).unwrap();
        assert_eq!(e.name, "e-");
        assert!((e.mass - M_ELECTRON).abs() < 1e-6);
        assert!((e.charge - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_lookup_positron() {
        let e = get_particle(PdgId::POSITRON).unwrap();
        assert_eq!(e.name, "e+");
        assert!((e.charge - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lookup_photon() {
        let g = get_particle(PdgId::PHOTON).unwrap();
        assert_eq!(g.mass, 0.0);
        assert!(g.pdg_id.is_self_conjugate());
    }

    #[test]
    fn test_higgs_mass() {
        let h = get_particle(PdgId::HIGGS).unwrap();
        assert!((h.mass - M_H).abs() < 0.01);
    }

    #[test]
    fn test_top_quark() {
        let t = get_particle(PdgId::TOP).unwrap();
        assert!((t.mass - M_TOP).abs() < 0.01);
        assert!(t.width.is_some());
    }
}
