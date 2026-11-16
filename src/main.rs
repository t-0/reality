use clap::Parser;
use reality::cli::commands::{Cli, Commands};
use reality::cli::display;

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Particles { filter } => {
            display::print_particles(&filter);
        }
        Commands::CrossSection {
            process,
            sqrt_s,
            n_points,
        } => {
            display::print_cross_section(&process, sqrt_s, n_points);
        }
        Commands::Decay { particle } => {
            display::print_decay_table(&particle);
        }
        Commands::Generate {
            process,
            sqrt_s,
            n_events,
            format,
            seed,
        } => {
            display::print_generated_events(&process, sqrt_s, n_events, &format, seed);
        }
        Commands::Coupling { name, scale } => {
            display::print_coupling(&name, scale);
        }
        Commands::Unification { susy_scale } => {
            display::print_unification(susy_scale);
        }
        Commands::Amend { file, loops } => {
            display::print_amendment(&file, loops);
        }
    }
}
