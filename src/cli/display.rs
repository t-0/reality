use comfy_table::{Cell, Table};

use crate::couplings::constants::*;
use crate::couplings::running;
use crate::cross_section::total::total_cross_section_2to2_with_points;
use crate::decay::tables;
use crate::matrix_elements::ee_mumu::EpEmToMuPMuM;
use crate::matrix_elements::ee_qqbar::EpEmToQQbar;
use crate::matrix_elements::ee_ww::EpEmToWpWm;
use crate::monte_carlo::sampler::EventGenerator;
use crate::particles;

pub fn print_particles(filter: &str) {
    let mut table = Table::new();
    table.set_header(vec![
        "Name", "PDG ID", "Mass (GeV)", "Charge", "Spin", "Color", "Width (GeV)",
    ]);

    let mut items: Vec<_> = match filter {
        "quarks" => particles::quarks().collect(),
        "leptons" => particles::leptons().collect(),
        "bosons" => particles::bosons().collect(),
        _ => {
            let mut all: Vec<_> = particles::all_particles()
                .filter(|p| !p.pdg_id.is_anti())
                .collect();
            all.sort_by_key(|p| p.pdg_id.0);
            all
        }
    };

    items.sort_by_key(|p| p.pdg_id.0);

    for p in &items {
        table.add_row(vec![
            Cell::new(p.name),
            Cell::new(format!("{}", p.pdg_id)),
            Cell::new(format!("{:.6}", p.mass)),
            Cell::new(format!("{:+.2}", p.charge)),
            Cell::new(format!("{}", p.spin)),
            Cell::new(format!("{}", p.color)),
            Cell::new(
                p.width
                    .map_or("stable".to_string(), |w| format!("{:.4e}", w)),
            ),
        ]);
    }

    println!("Standard Model Particles ({})\n", filter);
    println!("{table}");
}

pub fn print_cross_section(process: &str, sqrt_s: f64, n_points: usize) {
    let sigma = match process {
        "ee_mumu" => {
            let me = EpEmToMuPMuM::qed_only();
            total_cross_section_2to2_with_points(&me, sqrt_s, n_points)
        }
        "ee_mumu_ew" => {
            let me = EpEmToMuPMuM::full_electroweak();
            total_cross_section_2to2_with_points(&me, sqrt_s, n_points)
        }
        "ee_qqbar" | "ee_uubar" => {
            let me = EpEmToQQbar::up();
            total_cross_section_2to2_with_points(&me, sqrt_s, n_points)
        }
        "ee_ddbar" => {
            let me = EpEmToQQbar::down();
            total_cross_section_2to2_with_points(&me, sqrt_s, n_points)
        }
        "ee_ccbar" => {
            let me = EpEmToQQbar::charm();
            total_cross_section_2to2_with_points(&me, sqrt_s, n_points)
        }
        "ee_bbbar" => {
            let me = EpEmToQQbar::bottom();
            total_cross_section_2to2_with_points(&me, sqrt_s, n_points)
        }
        "ee_ww" => {
            let me = EpEmToWpWm;
            total_cross_section_2to2_with_points(&me, sqrt_s, n_points)
        }
        _ => {
            eprintln!("Unknown process '{}'. Available:", process);
            eprintln!("  ee_mumu     - e+e- -> mu+mu- (QED)");
            eprintln!("  ee_mumu_ew  - e+e- -> mu+mu- (QED + Z)");
            eprintln!("  ee_uubar    - e+e- -> u ubar");
            eprintln!("  ee_ddbar    - e+e- -> d dbar");
            eprintln!("  ee_ccbar    - e+e- -> c cbar");
            eprintln!("  ee_bbbar    - e+e- -> b bbar");
            eprintln!("  ee_ww       - e+e- -> W+W-");
            return;
        }
    };

    println!("Process: {}", process);
    println!("sqrt(s) = {:.2} GeV", sqrt_s);
    println!("sigma   = {:.6e} nb", sigma);
    println!("        = {:.6e} pb", sigma * 1e3);
    println!("        = {:.6e} fb", sigma * 1e6);
}

pub fn print_decay_table(particle: &str) {
    let table = match particle.to_lowercase().as_str() {
        "w" | "w+" | "w-" => tables::w_boson_decay_table(),
        "z" | "z0" => tables::z_boson_decay_table(),
        "h" | "higgs" => tables::higgs_decay_table(),
        "t" | "top" => tables::top_quark_decay_table(),
        _ => {
            eprintln!("Unknown particle '{}'. Available: w, z, higgs, top", particle);
            return;
        }
    };

    println!("{} Decay Table", table.parent_name);
    println!("Total width: {:.4e} GeV", table.total_width);
    if table.total_width > 0.0 {
        let tau_gev = 1.0 / table.total_width;
        let tau_s = tau_gev * 6.582e-25;
        println!("Lifetime:    {:.4e} s", tau_s);
    }
    println!();

    let mut display = Table::new();
    display.set_header(vec!["Channel", "Partial Width (GeV)", "BR (%)"]);

    let mut channels = table.channels.clone();
    channels.sort_by(|a, b| b.branching_ratio.partial_cmp(&a.branching_ratio).unwrap());

    for ch in &channels {
        display.add_row(vec![
            Cell::new(&ch.label),
            Cell::new(format!("{:.4e}", ch.partial_width)),
            Cell::new(format!("{:.3}", ch.branching_ratio * 100.0)),
        ]);
    }

    println!("{display}");
}

pub fn print_generated_events(
    process: &str,
    sqrt_s: f64,
    n_events: usize,
    format: &str,
    seed: Option<u64>,
) {
    let mut rng = if let Some(s) = seed {
        use rand::SeedableRng;
        rand::rngs::StdRng::seed_from_u64(s)
    } else {
        use rand::SeedableRng;
        rand::rngs::StdRng::from_os_rng()
    };

    match process {
        "ee_mumu" => {
            let me = EpEmToMuPMuM::qed_only();
            generate_and_print(&me, sqrt_s, n_events, format, &mut rng);
        }
        "ee_mumu_ew" => {
            let me = EpEmToMuPMuM::full_electroweak();
            generate_and_print(&me, sqrt_s, n_events, format, &mut rng);
        }
        _ => {
            eprintln!(
                "Unknown process '{}'. Available for generation: ee_mumu, ee_mumu_ew",
                process
            );
        }
    }
}

fn generate_and_print(
    me: &impl crate::matrix_elements::MatrixElement,
    sqrt_s: f64,
    n_events: usize,
    format: &str,
    rng: &mut impl rand::Rng,
) {
    let mut gen = EventGenerator::new(me, sqrt_s);
    gen.warmup(rng, 10000);
    let sample = gen.generate_sample(rng, n_events);

    match format {
        "json" => {
            println!("{}", serde_json::to_string_pretty(&sample).unwrap());
        }
        "table" => {
            println!("Process: {}", sample.process);
            println!("sqrt(s) = {:.2} GeV", sample.sqrt_s);
            println!(
                "Cross section: {:.4e} +/- {:.4e} nb",
                sample.cross_section_nb, sample.cross_section_error_nb
            );
            println!("Events generated: {}\n", sample.n_events);

            let mut table = Table::new();
            table.set_header(vec!["Event", "PDG", "Status", "E", "px", "py", "pz"]);

            for event in &sample.events {
                for p in &event.particles {
                    table.add_row(vec![
                        Cell::new(format!("{}", event.event_number)),
                        Cell::new(format!("{}", p.pdg_id)),
                        Cell::new(format!("{:?}", p.status)),
                        Cell::new(format!("{:.4}", p.momentum[0])),
                        Cell::new(format!("{:.4}", p.momentum[1])),
                        Cell::new(format!("{:.4}", p.momentum[2])),
                        Cell::new(format!("{:.4}", p.momentum[3])),
                    ]);
                }
            }
            println!("{table}");
        }
        _ => {
            eprintln!("Unknown format '{}'. Available: json, table", format);
        }
    }
}

pub fn print_coupling(name: &str, scale: Option<f64>) {
    match name {
        "alpha_em" | "alpha" => {
            println!("alpha_em (q^2 = 0) = {:.10}", ALPHA_EM);
            println!("alpha_em (M_Z)     = {:.10}", ALPHA_EM_MZ);
            println!("e = sqrt(4*pi*alpha) = {:.6}", e_charge());
        }
        "alpha_s" => {
            if let Some(q) = scale {
                let a = running::alpha_s(q);
                println!("alpha_s(Q = {:.2} GeV) = {:.6}", q, a);
                println!("n_f = {} active flavors", running::n_active_flavors(q));
            } else {
                println!("alpha_s(M_Z = {:.4} GeV) = {:.6}", M_Z, ALPHA_S_MZ);
                println!("\nRunning values:");
                for q in [1.0, 2.0, 5.0, 10.0, 50.0, M_Z, 200.0, 1000.0] {
                    println!(
                        "  alpha_s({:>8.2} GeV) = {:.6}  (n_f = {})",
                        q,
                        running::alpha_s(q),
                        running::n_active_flavors(q)
                    );
                }
            }
        }
        "sin2w" | "sin2_theta_w" | "weinberg" => {
            println!("sin^2(theta_W) = {:.6}", SIN2_THETA_W);
            println!("cos^2(theta_W) = {:.6}", COS2_THETA_W);
            println!("sin(theta_W)   = {:.6}", SIN_THETA_W);
            println!("cos(theta_W)   = {:.6}", COS_THETA_W);
        }
        "gf" | "fermi" => {
            println!("G_F = {:.6e} GeV^-2", G_FERMI);
            println!("v_EW = {:.2} GeV", V_EW);
        }
        "g_weak" | "gw" => {
            println!("g_weak = e / sin(theta_W) = {:.6}", g_weak());
        }
        _ => {
            eprintln!("Unknown coupling '{}'. Available:", name);
            eprintln!("  alpha_em  - Fine structure constant");
            eprintln!("  alpha_s   - Strong coupling (use --scale for running)");
            eprintln!("  sin2w     - Weinberg angle");
            eprintln!("  gf        - Fermi constant");
            eprintln!("  gw        - Weak coupling constant");
        }
    }
}
