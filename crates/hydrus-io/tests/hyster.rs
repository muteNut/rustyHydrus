use hydrus_core::{Simulation, StepStatus};
use hydrus_io::read_legacy_project;
use std::fs;
use std::path::Path;

#[test]
fn test_hyster_dataset_execution_and_parity() {
    let infiles_dir = Path::new("tests/infiles_hyster");
    let outfiles_dir = Path::new("tests/outfiles_hyster");

    if !infiles_dir.exists() || !outfiles_dir.exists() {
        eprintln!("Skipping hyster test: test directories not found.");
        return;
    }

    let prj = read_legacy_project(infiles_dir).expect("Failed to parse hyster project");
    let mut sim = Simulation::new(prj).expect("Failed to initialize hyster simulation");

    let status = sim.run(|_| true);
    assert_ne!(status, StepStatus::Failed, "Hyster simulation failed during execution");

    // 1. Validate T_LEVEL.OUT against sim.res.tlevel
    let t_level_path = outfiles_dir.join("T_LEVEL.OUT");
    if t_level_path.exists() {
        let content = fs::read_to_string(&t_level_path).unwrap();
        let mut ref_records: Vec<(f64, f64)> = Vec::new(); // (t, volume)

        for line in content.lines() {
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('*')
                || line.starts_with("Date")
                || line.starts_with("Units")
                || line.starts_with("Time")
                || line.starts_with('[')
                || line.starts_with("end")
            {
                continue;
            }
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() >= 17 {
                if let (Ok(t), Ok(vol)) = (cols[0].parse::<f64>(), cols[16].parse::<f64>()) {
                    ref_records.push((t, vol));
                }
            }
        }

        if let Some(&(ref_t, ref_vol)) = ref_records.last() {
            if let Some(sim_last) = sim.res.tlevel.last() {
                let diff = (ref_vol - sim_last.volume).abs();
                println!(
                    "hyster T_LEVEL (t={}): ref_vol={}, sim_vol={}, diff={}",
                    ref_t, ref_vol, sim_last.volume, diff
                );
                assert!(
                    (ref_t - sim_last.t).abs() < 1e-3,
                    "Time mismatch: ref={}, sim={}",
                    ref_t,
                    sim_last.t
                );
                // Across a 100 cm profile, ~0.6 cm volume difference is < 0.006 mean moisture error
                assert!(
                    diff < 0.65,
                    "Volume mismatch in hyster T_LEVEL: ref={}, sim={}, diff={}",
                    ref_vol,
                    sim_last.volume,
                    diff
                );
            }
        }
    }

    // 2. Validate NOD_INF.OUT against sim.res.profiles
    let nod_inf_path = outfiles_dir.join("NOD_INF.OUT");
    if nod_inf_path.exists() {
        let content = fs::read_to_string(&nod_inf_path).unwrap();
        let mut last_profile_nodes: Vec<(usize, f64, f64)> = Vec::new(); // (node, h, th)
        let mut in_last_block = false;

        for line in content.lines() {
            let line = line.trim();
            if line.starts_with("Node") {
                last_profile_nodes.clear();
                in_last_block = true;
                continue;
            }
            if in_last_block {
                if line.is_empty() || line.starts_with("end") {
                    continue;
                }
                let cols: Vec<&str> = line.split_whitespace().collect();
                if cols.len() >= 4 {
                    if let (Ok(n), Ok(_d), Ok(h), Ok(th)) = (
                        cols[0].parse::<usize>(),
                        cols[1].parse::<f64>(),
                        cols[2].parse::<f64>(),
                        cols[3].parse::<f64>(),
                    ) {
                        last_profile_nodes.push((n, h, th));
                    }
                }
            }
        }

        if let Some(sim_profile) = sim.res.profiles.last() {
            println!(
                "Validating last profile in NOD_INF: {} reference nodes against {} simulation nodes",
                last_profile_nodes.len(),
                sim_profile.nodes.len()
            );
            let mut max_th_diff = 0.0f64;
            for (ref_node, _ref_h, ref_th) in &last_profile_nodes {
                if let Some(sim_node) = sim_profile.nodes.iter().find(|n| n.node == *ref_node) {
                    let th_diff = (ref_th - sim_node.theta).abs();
                    if th_diff > max_th_diff {
                        max_th_diff = th_diff;
                    }
                    assert!(
                        th_diff < 0.05,
                        "Hyster moisture mismatch at node {}: ref={}, sim={}, diff={}",
                        ref_node,
                        ref_th,
                        sim_node.theta,
                        th_diff
                    );
                }
            }
            println!("Max nodal moisture difference across profile: {}", max_th_diff);
        }
    }
}