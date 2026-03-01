use std::collections::HashMap;
use std::sync::LazyLock;

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
            mass: 0.00216,
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
            mass: 0.00467,
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
            mass: 1.27,
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
            mass: 0.0934,
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
            mass: 172.57,
            width: Some(1.42),
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
            mass: 4.18,
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
            mass: 0.000511,
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
            mass: 0.10566,
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
            mass: 1.777,
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
            mass: 91.1876,
            width: Some(2.4955),
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
            mass: 80.377,
            width: Some(2.085),
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
            mass: 125.20,
            width: Some(0.00407),
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
        assert!((e.mass - 0.000511).abs() < 1e-6);
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
        assert!((h.mass - 125.20).abs() < 0.01);
    }

    #[test]
    fn test_top_quark() {
        let t = get_particle(PdgId::TOP).unwrap();
        assert!((t.mass - 172.57).abs() < 0.01);
        assert!(t.width.is_some());
    }
}
