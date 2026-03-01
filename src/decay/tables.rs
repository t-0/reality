use crate::couplings::ckm::{ckm, Generation};
use crate::couplings::constants::*;
use crate::decay::partial_width;
use crate::decay::total_width::DecayTable;
use crate::matrix_elements::higgs_decay_gammagamma::HiggsToGammaGamma;
use crate::particles::PdgId;

/// Build the W boson decay table.
pub fn w_boson_decay_table() -> DecayTable {
    let v = ckm();

    let mut channels = Vec::new();

    // Leptonic channels: W -> l nu (3 generations)
    let lepton_names = ["e+ nu_e", "mu+ nu_mu", "tau+ nu_tau"];
    for name in lepton_names {
        channels.push((vec![], name.to_string(), partial_width::w_to_lepton_nu()));
    }

    // Hadronic channels: W -> u_i d_j (kinematically allowed)
    let up_quarks = [
        (Generation::First, "u"),
        (Generation::Second, "c"),
    ];
    let down_quarks = [
        (Generation::First, "d"),
        (Generation::Second, "s"),
        (Generation::Third, "b"),
    ];

    for &(up_gen, up_name) in &up_quarks {
        for &(down_gen, down_name) in &down_quarks {
            let v_sq = v.abs_squared(up_gen, down_gen);
            let label = format!("{} {}", up_name, down_name);
            channels.push((vec![], label, partial_width::w_to_quarks(v_sq)));
        }
    }

    DecayTable::from_channels(PdgId::W_BOSON, "W+", channels)
}

/// Build the Z boson decay table.
pub fn z_boson_decay_table() -> DecayTable {
    let mut channels = Vec::new();

    // Neutrinos: Z -> nu nu_bar (3 generations)
    // T3 = +1/2, Q = 0, N_c = 1
    let nu_width = partial_width::z_to_ffbar(0.5, 0.0, 1.0);
    for name in ["nu_e nu_e~", "nu_mu nu_mu~", "nu_tau nu_tau~"] {
        channels.push((vec![], name.to_string(), nu_width));
    }

    // Charged leptons: Z -> l+ l- (3 generations)
    // T3 = -1/2, Q = -1, N_c = 1
    let lepton_width = partial_width::z_to_ffbar(-0.5, -1.0, 1.0);
    for name in ["e+ e-", "mu+ mu-", "tau+ tau-"] {
        channels.push((vec![], name.to_string(), lepton_width));
    }

    // Up-type quarks: T3 = +1/2, Q = +2/3, N_c = 3
    let up_width = partial_width::z_to_ffbar(0.5, 2.0 / 3.0, N_C);
    for name in ["u u~", "c c~"] {
        channels.push((vec![], name.to_string(), up_width));
    }

    // Down-type quarks: T3 = -1/2, Q = -1/3, N_c = 3
    let down_width = partial_width::z_to_ffbar(-0.5, -1.0 / 3.0, N_C);
    for name in ["d d~", "s s~", "b b~"] {
        channels.push((vec![], name.to_string(), down_width));
    }

    DecayTable::from_channels(PdgId::Z_BOSON, "Z", channels)
}

/// Build the Higgs boson decay table.
pub fn higgs_decay_table() -> DecayTable {
    let mut channels = Vec::new();

    // H -> bb (dominant)
    channels.push((vec![], "b b~".to_string(), partial_width::h_to_ffbar(M_BOTTOM, N_C)));

    // H -> cc
    channels.push((vec![], "c c~".to_string(), partial_width::h_to_ffbar(M_CHARM, N_C)));

    // H -> tau+ tau-
    channels.push((vec![], "tau+ tau-".to_string(), partial_width::h_to_ffbar(M_TAU, 1.0)));

    // H -> mu+ mu-
    channels.push((vec![], "mu+ mu-".to_string(), partial_width::h_to_ffbar(0.10566, 1.0)));

    // H -> gamma gamma (loop)
    channels.push((vec![], "gamma gamma".to_string(), HiggsToGammaGamma.partial_width()));

    // H -> WW* and H -> ZZ* are important but require off-shell treatment
    // For now, add approximate values from known total width
    // H -> WW*: ~21.5% of total width
    // H -> ZZ*: ~2.6% of total width
    // Using total width = 4.07 MeV
    let total_approx = 0.00407;
    channels.push((vec![], "W W*".to_string(), 0.215 * total_approx));
    channels.push((vec![], "Z Z*".to_string(), 0.026 * total_approx));

    // H -> gg (loop, via top)
    // ~8.6% of total width
    channels.push((vec![], "g g".to_string(), 0.086 * total_approx));

    DecayTable::from_channels(PdgId::HIGGS, "H", channels)
}

/// Build the top quark decay table.
pub fn top_quark_decay_table() -> DecayTable {
    let channels = vec![
        (vec![], "W+ b".to_string(), partial_width::top_to_wb()),
    ];
    DecayTable::from_channels(PdgId::TOP, "t", channels)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_w_total_width() {
        let table = w_boson_decay_table();
        // Expected: ~2.085 GeV
        assert!(
            (table.total_width - 2.085).abs() < 0.15,
            "Gamma_W = {} GeV",
            table.total_width
        );
    }

    #[test]
    fn test_w_leptonic_branching() {
        let table = w_boson_decay_table();
        // Each leptonic channel should be ~10.8%
        let br_enu = table.channels[0].branching_ratio;
        assert!(
            (br_enu - 0.108).abs() < 0.02,
            "BR(W -> e nu) = {}",
            br_enu
        );
    }

    #[test]
    fn test_z_total_width() {
        let table = z_boson_decay_table();
        // Expected: ~2.4955 GeV
        assert!(
            (table.total_width - 2.4955).abs() < 0.1,
            "Gamma_Z = {} GeV",
            table.total_width
        );
    }

    #[test]
    fn test_z_ee_branching() {
        let table = z_boson_decay_table();
        // BR(Z -> e+e-) ~ 3.36%
        let ee_channel = table.channels.iter().find(|c| c.label == "e+ e-").unwrap();
        assert!(
            (ee_channel.branching_ratio - 0.0336).abs() < 0.005,
            "BR(Z -> ee) = {}",
            ee_channel.branching_ratio
        );
    }

    #[test]
    fn test_z_invisible_width() {
        let table = z_boson_decay_table();
        // Invisible width (3 neutrinos) ~ 499 MeV
        let inv_width: f64 = table.channels.iter()
            .filter(|c| c.label.contains("nu_"))
            .map(|c| c.partial_width)
            .sum();
        assert!(
            (inv_width - 0.499).abs() < 0.02,
            "Gamma_inv = {} GeV",
            inv_width
        );
    }

    #[test]
    fn test_higgs_bb_dominant() {
        let table = higgs_decay_table();
        let bb = table.channels.iter().find(|c| c.label == "b b~").unwrap();
        // H -> bb should be the dominant channel, BR ~ 58%
        assert!(
            bb.branching_ratio > 0.40,
            "BR(H -> bb) = {} should be dominant",
            bb.branching_ratio
        );
    }

    #[test]
    fn test_top_width() {
        let table = top_quark_decay_table();
        assert!(
            (table.total_width - 1.4).abs() < 0.3,
            "Gamma_top = {} GeV",
            table.total_width
        );
        // t -> Wb should be ~100%
        assert!((table.channels[0].branching_ratio - 1.0).abs() < 1e-10);
    }
}
