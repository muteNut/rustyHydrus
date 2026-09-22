use hydrus_core::{Simulation, StepStatus};
use std::path::PathBuf;

fn usage() -> ! {
    eprintln!("hydrus1d-cli <project-dir | project.h1dr> [-o outdir] [--quiet]\n  Runs a HYDRUS-1D project (legacy directory or native .h1dr JSON) and writes CSV results.");
    std::process::exit(2)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let mut input: Option<PathBuf> = None;
    let mut out: Option<PathBuf> = None;
    let mut quiet = false;
    while let Some(a) = args.next() {
        match a.as_str() {
            "-o" => out = args.next().map(PathBuf::from),
            "--quiet" | "-q" => quiet = true,
            "-h" | "--help" => usage(),
            _ => input = Some(PathBuf::from(a)),
        }
    }
    let input = input.unwrap_or_else(|| usage());
    let prj = if input.is_dir() {
        hydrus_io::read_legacy_project(&input)
    } else {
        hydrus_io::load_project(&input)
    };
    let prj = match prj {
        Ok(p) => p,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1)
        }
    };
    let mut sim = match Simulation::new(prj) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1)
        }
    };
    let mut last = -1.0;
    let status = sim.run(|s| {
        if !quiet {
            let pct = 100.0 * (s.t - s.t_init) / (s.t_max - s.t_init);
            if pct - last >= 5.0 {
                eprint!("\r{:5.1}%  t = {:.5}  dt = {:.3e}", pct, s.t, s.dt);
                last = pct;
            }
        }
        true
    });
    if !quiet {
        eprintln!();
    }
    for m in &sim.res.messages {
        eprintln!("{}", m);
    }
    let outdir = out.unwrap_or_else(|| if input.is_dir() { input.join("rust_out") } else { input.with_extension("results") });
    if let Err(e) = hydrus_io::writer::write_all(&sim.res, &outdir) {
        eprintln!("{}", e);
        std::process::exit(1)
    }
    eprintln!("Results written to {}", outdir.display());
    if status == StepStatus::Failed {
        std::process::exit(3)
    }
}
