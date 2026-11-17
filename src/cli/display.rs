use comfy_table::{Cell, Table};

use crate::couplings::constants::*;
use crate::couplings::running;
use crate::couplings::unification;
use crate::cross_section::total::total_cross_section_2to2_with_points;
use crate::decay::tables;
use crate::matrix_elements::ee_mumu::EpEmToMuPMuM;
use crate::matrix_elements::ee_qqbar::EpEmToQQbar;
use crate::matrix_elements::ee_ww::EpEmToWpWm;
use crate::model::{self, gates::Verdict, proton, spec::COUPLING_NAMES, Amendment};
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
        let tau_s = tau_gev * HBAR_GEV_S;
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

pub fn print_unification(susy_scale: Option<f64>) {
    let names = ["1/alpha_1 (U(1)_Y)", "1/alpha_2 (SU(2)_L)", "1/alpha_3 (SU(3)_c)"];
    let short = ["alpha_1", "alpha_2", "alpha_3"];

    match susy_scale {
        Some(ms) => println!("Model: SM below {:.0} GeV, MSSM above (one-loop)", ms),
        None => println!("Model: Standard Model (one-loop)"),
    }
    println!("alpha_1 uses GUT normalization: alpha_1 = (5/3) alpha_em / cos^2(theta_W)\n");

    let mut table = Table::new();
    table.set_header(vec!["Q (GeV)", names[0], names[1], names[2]]);
    for exp in [2, 4, 6, 8, 10, 12, 14, 15, 16, 17, 18, 19] {
        let q = if exp == 2 { M_Z } else { 10f64.powi(exp) };
        let a = unification::inverse_couplings(q, susy_scale);
        table.add_row(vec![
            Cell::new(format!("{:.3e}", q)),
            Cell::new(format!("{:.2}", a[0])),
            Cell::new(format!("{:.2}", a[1])),
            Cell::new(format!("{:.2}", a[2])),
        ]);
    }
    println!("{table}");

    println!("\nPairwise crossings:");
    for c in unification::crossings(susy_scale).iter().flatten() {
        println!(
            "  {} = {}  at Q = {:.2e} GeV  (1/alpha = {:.2})",
            short[c.pair.0], short[c.pair.1], c.scale, c.inverse_alpha
        );
    }

    let (q, spread) = unification::closest_approach(susy_scale);
    println!(
        "\nClosest approach: Q = {:.2e} GeV, spread in 1/alpha = {:.2}",
        q, spread
    );
    if spread < 1.0 {
        println!("The three couplings meet within one-loop accuracy.");
    } else {
        println!("The three couplings miss each other. Try --susy-scale 1000.");
    }
}

fn fmt_b(b: [f64; 3]) -> String {
    format!("({:.3}, {:.3}, {:.3})", b[0], b[1], b[2])
}

pub fn print_amendment(path: &std::path::Path, loops: u8) {
    let amendment = match Amendment::from_file(path) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };
    let report = model::evaluate(&amendment, loops);
    let (i, j) = amendment.inputs;
    let k = amendment.predicted;

    println!("{}", amendment.name);
    if let Some(d) = &amendment.description {
        println!("{}", d);
    }
    println!(
        "\nInputs (measured): {}, {}    Held out (predicted): {}",
        COUPLING_NAMES[i], COUPLING_NAMES[j], COUPLING_NAMES[k]
    );

    if !amendment.added.is_empty() {
        let mut table = Table::new();
        table.set_header(vec!["Added field", "Kind", "SU(3)", "SU(2)", "Y", "Copies", "Mass (GeV)"]);
        for f in &amendment.added {
            table.add_row(vec![
                Cell::new(&f.name),
                Cell::new(format!("{:?}", f.kind)),
                Cell::new(f.su3.label()),
                Cell::new(f.su2),
                Cell::new(format!("{:+.4}", f.y)),
                Cell::new(f.copies),
                Cell::new(format!("{:.3e}", f.mass)),
            ]);
        }
        println!("{table}");
    }

    if !amendment.heavy.is_empty() {
        let mut table = Table::new();
        table.set_header(vec!["Heavy field (at M_X)", "Kind", "SU(3)", "SU(2)", "Y", "Copies", "M / M_X", "b_i"]);
        for h in &amendment.heavy {
            let f = &h.field;
            let ratio = if h.is_varied() {
                format!("{} - {}", h.ratio.0, h.ratio.1)
            } else {
                format!("{}", h.ratio.0)
            };
            table.add_row(vec![
                Cell::new(&f.name),
                Cell::new(format!("{:?}", f.kind)),
                Cell::new(f.su3.label()),
                Cell::new(f.su2),
                Cell::new(format!("{:+.4}", f.y)),
                Cell::new(f.copies),
                Cell::new(ratio),
                Cell::new(fmt_b(f.beta_contribution())),
            ]);
        }
        println!("{table}");
    }

    println!("\nOne-loop beta coefficients (b1, b2, b3) by energy range:");
    let fields = amendment.all_fields();
    for s in &report.segments {
        let active: Vec<_> = fields.iter().filter(|f| f.mass <= s.lo * (1.0 + 1e-9)).collect();
        let susy = model::two_loop::is_supersymmetric(&active);
        println!(
            "  {:.2e} - {:.2e} GeV   {}{}",
            s.lo,
            s.hi,
            fmt_b(s.b),
            if susy { "   supersymmetric" } else { "" }
        );
    }

    println!("\nConsistency gates:");
    for g in &report.gates {
        let mark = match g.verdict {
            Verdict::Pass => "PASS",
            Verdict::Warn => "WARN",
            Verdict::Fail => "FAIL",
        };
        println!("  [{}] {:<20} {}", mark, g.name, g.detail);
    }
    match report.landau_pole {
        Some(q) => println!("  [FAIL] {:<20} a coupling diverges at {:.2e} GeV", "perturbativity", q),
        None => println!("  [PASS] {:<20} all couplings finite up to 1e19 GeV", "perturbativity"),
    }
    if !report.gates_passed() {
        println!("\nThe model is inconsistent; its predictions are shown for reference only.");
    }

    println!(
        "\nPrediction ({}-loop running{}):",
        loops,
        if loops >= 2 {
            format!(", top Yukawa: {:?}", amendment.top_yukawa).to_lowercase()
                .replace("auto", "on")
        } else {
            String::new()
        }
    );
    match &report.prediction {
        Some(p) => {
            println!(
                "  {} = {} at M_GUT = {:.2e} GeV, 1/alpha_GUT = {:.2}",
                COUPLING_NAMES[i], COUPLING_NAMES[j], p.m_gut, p.inverse_alpha_gut
            );
            println!(
                "  {}(M_Z) predicted = {:.4}   measured = {:.4}   ({:+.1}%)",
                COUPLING_NAMES[k],
                p.alpha_predicted(),
                p.alpha_measured(),
                100.0 * p.relative_error()
            );
            println!(
                "  1/{} predicted = {:.3}   measured = {:.3}",
                COUPLING_NAMES[k], p.inverse_predicted, p.inverse_measured
            );
            println!(
                "  pull = {:+.1} sigma (experiment + {:.0}% theory), {:+.1} sigma (experiment only)",
                p.pull,
                100.0 * p.theory_rel,
                p.pull_experimental
            );
            if let Some(b) = &report.thresholds {
                println!(
                    "  GUT thresholds ({} corners): {}(M_Z) = {:.4} - {:.4}, M_X = {:.2e} - {:.2e} GeV",
                    b.corners, COUPLING_NAMES[k], b.alpha_min, b.alpha_max, b.m_gut_min, b.m_gut_max
                );
                println!(
                    "  Measured value is {} the threshold band.",
                    if b.contains(p.alpha_measured()) { "inside" } else { "outside" }
                );
            }
        }
        None => println!(
            "  {} and {} never meet below 1e19 GeV (or only beyond a Landau pole): no prediction.",
            COUPLING_NAMES[i], COUPLING_NAMES[j]
        ),
    }

    if let Some(pd) = &report.proton {
        let (lo, hi) = pd.range_years();
        println!("\nProton decay (p -> e+ pi0 via X/Y gauge bosons at M_GUT):");
        println!(
            "  tau = {:.1e} years (range {:.0e} - {:.0e})",
            pd.lifetime_years, lo, hi
        );
        println!(
            "  Super-K limit: > {:.1e} years    Hyper-K reach: ~{:.0e} years",
            proton::SUPER_K_LIMIT_YEARS,
            proton::HYPER_K_REACH_YEARS
        );
        let verdict = match pd.status {
            proton::Status::Excluded => "EXCLUDED: protons would already have been seen decaying.".to_string(),
            proton::Status::Borderline => "BORDERLINE: the estimate straddles the Super-K limit.".to_string(),
            proton::Status::TestableAtHyperK => format!(
                "TESTABLE. Kill condition: Hyper-K sees no p -> e+ pi0 down to ~{:.0e} years.",
                proton::HYPER_K_REACH_YEARS
            ),
            proton::Status::OutOfReach => {
                "ALLOWED, beyond Hyper-K's reach. No near-term kill condition from this channel.".to_string()
            }
        };
        println!("  {}", verdict);
    }

    if let Some(d) = &report.dimension5 {
        let (lo, hi) = d.range_years();
        println!("\nProton decay (p -> K+ nu-bar via colour-triplet higgsinos, supersymmetric models):");
        println!(
            "  M_HC = {:.2e} GeV, M_S = {:.2e} GeV, mu_H = {:.2e} GeV, tan(beta) = {}",
            d.m_hc, d.m_s, d.mu_h, d.tan_beta
        );
        println!("  tau = {:.1e} years (range {:.0e} - {:.0e})", d.lifetime_years, lo, hi);
        if let Some((a, b)) = report.thresholds.and_then(|t| t.dimension5_years) {
            println!("  Over the GUT-threshold ranges: {:.1e} - {:.1e} years", a, b);
        }
        println!(
            "  Super-K limit: > {:.1e} years    Hyper-K reach: ~{:.0e} years",
            proton::SUPER_K_LIMIT_K_NU_YEARS,
            proton::HYPER_K_REACH_K_NU_YEARS
        );
        let verdict = match d.status {
            proton::Status::Excluded => "EXCLUDED: p -> K+ nu-bar would already have been seen.".to_string(),
            proton::Status::Borderline => "BORDERLINE: the estimate straddles the Super-K limit.".to_string(),
            proton::Status::TestableAtHyperK => format!(
                "TESTABLE. Kill condition: Hyper-K sees no p -> K+ nu-bar down to ~{:.0e} years.",
                proton::HYPER_K_REACH_K_NU_YEARS
            ),
            proton::Status::OutOfReach => "ALLOWED, beyond Hyper-K's reach in this channel.".to_string(),
        };
        println!("  {}", verdict);
        println!("  Estimate: Hisano, Kuwahara & Nagata, arXiv:1304.3651 eq. (19), higgsino dressing.");
    }

    let scales = amendment.mass_scales();
    println!("\nLedger:");
    println!(
        "  Parameters set by hand: {} mass scale{}{}",
        scales.len(),
        if scales.len() == 1 { "" } else { "s" },
        if scales.is_empty() {
            String::new()
        } else {
            format!(
                " ({})",
                scales.iter().map(|m| format!("{:.2e} GeV", m)).collect::<Vec<_>>().join(", ")
            )
        }
    );
    if report.dimension5.is_some() {
        println!(
            "  tan(beta): {}{} (used only for dimension-5 proton decay)",
            amendment.tan_beta.unwrap_or(model::predict::DEFAULT_TAN_BETA),
            if amendment.tan_beta.is_none() { ", default" } else { "" }
        );
    }
    let ratios = amendment.varied_ratio_count();
    if !amendment.heavy.is_empty() {
        println!(
            "  GUT-scale mass ratios: {} varied over a range (reported as a band, not fitted)",
            ratios
        );
    }
    println!("  Parameters fitted to the prediction: 0");
    let mut predicted = Vec::new();
    if report.prediction.is_some() {
        predicted.push(format!("{}(M_Z)", COUPLING_NAMES[k]));
    }
    if report.proton.is_some() {
        predicted.push("proton lifetime".to_string());
    }
    if report.dimension5.is_some() {
        predicted.push("p -> K+ nu-bar lifetime".to_string());
    }
    println!(
        "  Predictions: {}{}",
        predicted.len(),
        if predicted.is_empty() { String::new() } else { format!(" ({})", predicted.join(", ")) }
    );
    if loops >= 2 && amendment.heavy.is_empty() {
        println!("\nNote: sharp thresholds; no GUT-scale fields declared, so GUT threshold corrections are not included.");
    } else if loops >= 2 {
        println!("\nNote: sharp superpartner thresholds; GUT thresholds matched at one loop at M_X.");
    } else {
        println!("\nNote: one-loop running. Two-loop and threshold corrections shift alpha_s by a few percent.");
    }
}

pub fn print_scan(path: &std::path::Path, from: f64, to: f64, n: usize, loops: u8, format: &str) {
    use crate::model::scan;

    if !(from > 0.0 && to > from && n >= 2) {
        eprintln!("Need 0 < --from < --to and --points >= 2.");
        std::process::exit(1);
    }
    let amendment = match Amendment::from_file(path) {
        Ok(a) if !a.added.is_empty() => a,
        Ok(_) => {
            eprintln!("{} adds no fields, so there is no mass to scan.", path.display());
            std::process::exit(1);
        }
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };
    let k = amendment.predicted;
    let points = scan::scan(&amendment, loops, &scan::log_grid(from, to, n));

    if format == "csv" {
        println!("mass_gev,{0}_predicted,{0}_measured,pull,m_gut_gev,inverse_alpha_gut,proton_lifetime_years,{0}_band_min,{0}_band_max", COUPLING_NAMES[k]);
        for p in &points {
            match (&p.prediction, &p.proton) {
                (Some(pr), Some(pd)) => {
                    let (lo, hi) = p
                        .thresholds
                        .map(|b| (format!("{:.6}", b.alpha_min), format!("{:.6}", b.alpha_max)))
                        .unwrap_or_default();
                    println!(
                        "{:.6e},{:.6},{:.6},{:.4},{:.6e},{:.4},{:.6e},{},{}",
                        p.mass,
                        pr.alpha_predicted(),
                        pr.alpha_measured(),
                        pr.pull,
                        pr.m_gut,
                        pr.inverse_alpha_gut,
                        pd.lifetime_years,
                        lo,
                        hi
                    )
                }
                _ => println!("{:.6e},,,,,,,,", p.mass),
            }
        }
        return;
    }

    println!("{}  ({}-loop)", amendment.name, loops);
    println!("Scanning the lightest new mass; heavier fields keep their ratio to it.\n");
    let mut table = Table::new();
    table.set_header(vec![
        "Mass (GeV)".to_string(),
        format!("{} predicted", COUPLING_NAMES[k]),
        "Pull".to_string(),
        "M_GUT (GeV)".to_string(),
        "Proton lifetime (yr)".to_string(),
        "GUT threshold band".to_string(),
    ]);
    for p in &points {
        match (&p.prediction, &p.proton) {
            (Some(pr), Some(pd)) => table.add_row(vec![
                Cell::new(format!("{:.2e}", p.mass)),
                Cell::new(format!("{:.4}", pr.alpha_predicted())),
                Cell::new(format!("{:+.2}", pr.pull)),
                Cell::new(format!("{:.2e}", pr.m_gut)),
                Cell::new(format!("{:.1e}", pd.lifetime_years)),
                Cell::new(
                    p.thresholds
                        .map(|b| format!("{:.4} - {:.4}", b.alpha_min, b.alpha_max))
                        .unwrap_or_default(),
                ),
            ]),
            _ => table.add_row(vec![
                Cell::new(format!("{:.2e}", p.mass)),
                Cell::new("no unification"),
                Cell::new(""),
                Cell::new(""),
                Cell::new(""),
                Cell::new(""),
            ]),
        };
    }
    println!("{table}");

    if let Some(m) = points.iter().find_map(|p| p.prediction.as_ref()).map(|p| p.alpha_measured()) {
        println!("\nMeasured {}(M_Z) = {:.4}", COUPLING_NAMES[k], m);
    }
    match scan::best_fit(&points) {
        Some((m, pull)) if pull == 0.0 => println!("Best fit: {:.2e} GeV (pull crosses zero)", m),
        Some((m, pull)) => println!(
            "Closest: {:.2e} GeV with pull {:+.2} (never reaches zero in this range)",
            m, pull
        ),
        None => println!("No mass in this range gives unification."),
    }
    for n_sigma in [1.0, 2.0] {
        let ranges = scan::allowed_ranges(&points, n_sigma);
        let text = if ranges.is_empty() {
            "none".to_string()
        } else {
            ranges
                .iter()
                .map(|(a, b)| format!("{:.2e} - {:.2e} GeV", a, b))
                .collect::<Vec<_>>()
                .join(", ")
        };
        println!("Within {} sigma: {}", n_sigma, text);
    }
    if points.iter().any(|p| p.thresholds.is_some()) {
        let ranges = scan::band_ranges(&points);
        println!(
            "Measured value inside the GUT threshold band: {}",
            if ranges.is_empty() {
                "nowhere in this range".to_string()
            } else {
                ranges
                    .iter()
                    .map(|(a, b)| format!("{:.2e} - {:.2e} GeV (grid points)", a, b))
                    .collect::<Vec<_>>()
                    .join(", ")
            }
        );
    }
}

pub fn print_search(
    max_fields: usize,
    masses: Option<&str>,
    loops: u8,
    accept: f64,
    top: usize,
    write: Option<&std::path::Path>,
    robust_samples: usize,
) {
    use crate::model::robust::{robustness, RobustConfig};
    use crate::model::search::{self, SearchConfig};

    let mut cfg = SearchConfig {
        max_fields,
        loops,
        accept_sigma: accept,
        ..SearchConfig::default()
    };
    if let Some(list) = masses {
        let parsed: Result<Vec<f64>, _> = list.split(',').map(|m| m.trim().parse::<f64>()).collect();
        match parsed {
            Ok(mut m) if !m.is_empty() && m.iter().all(|&x| x > 0.0) => {
                m.sort_by(|a, b| a.partial_cmp(b).unwrap());
                m.dedup();
                cfg.masses = m;
            }
            _ => {
                eprintln!("--masses must be a comma-separated list of positive numbers, e.g. 1e3,1e6,1e9");
                std::process::exit(1);
            }
        }
    }
    if max_fields == 0 || max_fields > 4 {
        eprintln!("--max-fields must be between 1 and 4.");
        std::process::exit(1);
    }

    let menu = search::menu();
    println!(
        "Searching SM + up to {} of {} multiplets at {} masses ({:.0e} - {:.0e} GeV)...",
        max_fields,
        menu.len(),
        cfg.masses.len(),
        cfg.masses[0],
        cfg.masses[cfg.masses.len() - 1]
    );
    let mut result = search::search(&cfg);
    // Robustness for every accepted model, then sort each group by it.
    let robust: Vec<f64> = if robust_samples > 0 {
        let rcfg = RobustConfig { samples: robust_samples, loops, accept_sigma: accept, ..RobustConfig::default() };
        eprintln!("Measuring robustness of {} models ({} draws each)...", result.hits.len(), robust_samples);
        let mut scored: Vec<(f64, search::Hit)> = result
            .hits
            .drain(..)
            .map(|h| (robustness(&h.amendment, &rcfg).pass, h))
            .collect();
        scored.sort_by(|a, b| {
            a.1.fields()
                .cmp(&b.1.fields())
                .then(b.0.partial_cmp(&a.0).unwrap())
                .then(a.1.mass_scales.cmp(&b.1.mass_scales))
        });
        let (r, h): (Vec<f64>, Vec<search::Hit>) = scored.into_iter().unzip();
        result.hits = h;
        r
    } else {
        Vec::new()
    };
    println!(
        "Evaluated {} candidates at one loop; {} passed the screen; {} accepted at {} loop{}.",
        result.evaluated,
        result.screened,
        result.hits.len(),
        loops,
        if loops == 1 { "" } else { "s" }
    );
    println!(
        "Acceptance: |pull| <= {}, proton lifetime above Super-K, M_GUT <= 1e18 GeV, perturbative.\n",
        accept
    );
    if let Some(b) = &result.baseline {
        println!(
            "Baseline SM: alpha_s(M_Z) = {:.4}, pull {:+.1} sigma.\n",
            b.alpha_predicted(),
            b.pull
        );
    }

    let mut shown = Vec::new();
    for n in 1..=max_fields {
        let group: Vec<(usize, &search::Hit)> =
            result.hits.iter().enumerate().filter(|(_, h)| h.fields() == n).collect();
        println!(
            "{} new multiplet{}: {} accepted",
            n,
            if n == 1 { "" } else { "s" },
            group.len()
        );
        if group.is_empty() {
            println!();
            continue;
        }
        let mut table = Table::new();
        let mut header = vec![
            "Added (mass in GeV)",
            "alpha_s",
            "Pull",
            "M_GUT (GeV)",
            "tau_p (yr)",
            "Scales",
            "States",
            "Test",
        ];
        if !robust.is_empty() {
            header.push("Robust");
        }
        table.set_header(header);
        for &(idx, h) in group.iter().take(top) {
            let fields = h
                .amendment
                .added
                .iter()
                .map(|f| {
                    format!(
                        "{}{} @ {:.0e}",
                        if f.copies > 1 { format!("{} x ", f.copies) } else { String::new() },
                        f.name,
                        f.mass
                    )
                })
                .collect::<Vec<_>>()
                .join("\n");
            let mut tests = Vec::new();
            if h.proton.lifetime_years < proton::HYPER_K_REACH_YEARS * proton::UNCERTAINTY_FACTOR {
                tests.push("Hyper-K");
            }
            if h.lightest_coloured.map_or(false, |m| m <= 1e4) {
                tests.push("collider");
            }
            let mut row = vec![
                Cell::new(fields),
                Cell::new(format!("{:.4}", h.prediction.alpha_predicted())),
                Cell::new(format!("{:+.2}", h.prediction.pull)),
                Cell::new(format!("{:.1e}", h.prediction.m_gut)),
                Cell::new(format!("{:.0e}", h.proton.lifetime_years)),
                Cell::new(h.mass_scales),
                Cell::new(h.new_states),
                Cell::new(if tests.is_empty() { "-".to_string() } else { tests.join(", ") }),
            ];
            if let Some(r) = robust.get(idx) {
                row.push(Cell::new(format!("{:.0}%", 100.0 * r)));
            }
            table.add_row(row);
            shown.push(h);
        }
        println!("{table}\n");
    }
    if !result.hits.is_empty() {
        let mut counts = vec![0usize; menu.len()];
        for h in &result.hits {
            let mut seen: Vec<usize> = h.picks.iter().map(|p| p.0).collect();
            seen.dedup();
            for i in seen {
                counts[i] += 1;
            }
        }
        let mut order: Vec<usize> = (0..menu.len()).collect();
        order.sort_by(|&a, &b| counts[b].cmp(&counts[a]).then(a.cmp(&b)));
        let mut table = Table::new();
        table.set_header(vec!["Multiplet", "In accepted models", "Share", "Known as"]);
        for i in order {
            table.add_row(vec![
                Cell::new(menu[i].label),
                Cell::new(counts[i]),
                Cell::new(format!("{:.0}%", 100.0 * counts[i] as f64 / result.hits.len() as f64)),
                Cell::new(menu[i].note),
            ]);
        }
        println!("How often each multiplet appears among the {} accepted models:", result.hits.len());
        println!("{table}\n");
    }
    if robust.is_empty() {
        println!("Sorted by: fewest multiplets, then fewest mass scales, then fewest new states, then |pull|.");
    } else {
        println!("Sorted by: fewest multiplets, then robustness (share of {} random draws that still pass:", robust_samples);
        println!("masses within 1 decade, any below 1 TeV fails; SU(5) GUT fields within 0.1-10 x M_X).");
        println!(
            "Statistical error on robustness with {} draws: up to +-{:.0}%.",
            robust_samples,
            100.0 * 0.5 / (robust_samples as f64).sqrt()
        );
    }
    println!("Test: Hyper-K = proton decay within ~10x of Hyper-K reach; collider = coloured state at or below 10 TeV.");

    if let Some(dir) = write {
        if let Err(e) = std::fs::create_dir_all(dir) {
            eprintln!("cannot create {}: {}", dir.display(), e);
            std::process::exit(1);
        }
        let mut counts = std::collections::HashMap::new();
        for h in &shown {
            let rank = counts.entry(h.fields()).or_insert(0);
            *rank += 1;
            let path = dir.join(format!("search_{}f_{}.toml", h.fields(), rank));
            if let Err(e) = std::fs::write(&path, search::to_toml(h)) {
                eprintln!("cannot write {}: {}", path.display(), e);
                std::process::exit(1);
            }
        }
        println!("\nWrote {} model files to {}.", shown.len(), dir.display());
    }
}

pub fn print_robust(path: &std::path::Path, samples: usize, spread: f64, min_mass: f64, thresholds: f64, seed: u64) {
    use crate::model::robust::{robustness, RobustConfig};

    let amendment = match Amendment::from_file(path) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };
    if samples == 0 || spread < 0.0 || min_mass <= 0.0 || thresholds < 1.0 {
        eprintln!("Need --samples >= 1, --spread >= 0, --min-mass > 0 and --thresholds >= 1.");
        std::process::exit(1);
    }
    let cfg = RobustConfig {
        samples,
        spread_decades: spread,
        min_mass,
        threshold_range: thresholds,
        seed,
        ..RobustConfig::default()
    };
    let r = robustness(&amendment, &cfg);
    println!("{}", amendment.name);
    println!(
        "\nDraws: {}. Each added mass moves log-uniformly within {} decade{}; draws below {:.0e} GeV fail.",
        r.samples,
        spread,
        if spread == 1.0 { "" } else { "s" },
        min_mass
    );
    if r.generic_thresholds {
        println!(
            "GUT fields: generic {}SU(5) (colour-triplet Higgs, adjoint octet and triplet) at {}-{} x M_X.",
            if r.supersymmetric { "supersymmetric " } else { "non-SUSY " },
            1.0 / thresholds,
            thresholds
        );
    } else if !amendment.heavy.is_empty() {
        println!("GUT fields: the model's own [[heavy]] fields, drawn within their ranges.");
    } else {
        println!("GUT fields: none.");
    }
    println!(
        "\nRobustness: {:.0}% +- {:.0}% of draws pass everything.",
        100.0 * r.pass,
        100.0 * r.pass_error
    );
    println!("Failures (a draw can fail more than one):");
    println!("  new particle below {:.0e} GeV: {:.0}%", min_mass, 100.0 * r.fail_collider);
    println!("  alpha_s outside 1 sigma:     {:.0}%", 100.0 * r.fail_alpha);
    println!("  proton below Super-K limit:  {:.0}%", 100.0 * r.fail_proton);
    println!("  no unification below 1e18:   {:.0}%", 100.0 * r.fail_unify);
    if let Some([lo, mid, hi]) = r.alpha_quantiles {
        println!("  alpha_s over draws: {:.4} (median), 90% between {:.4} and {:.4}", mid, lo, hi);
    }
}
