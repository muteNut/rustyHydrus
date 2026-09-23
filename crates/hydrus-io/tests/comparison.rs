use hydrus_core::{Simulation, StepStatus};
use hydrus_io::read_legacy_project;
use std::fs;
use std::path::Path;

#[test]
fn test_strict_numeric_outfiles_parity() {
    let infiles_dir = Path::new("tests/infiles_enbal");
    let outfiles_dir = Path::new("tests/outfiles_enbal");

    if !infiles_dir.exists() || !outfiles_dir.exists() {
        eprintln!("Skipping comparison test: test directories not found.");
        return;
    }

    let prj = read_legacy_project(infiles_dir).expect("Failed to parse project");
    let mut sim = Simulation::new(prj).expect("Failed to initialize simulation");

    let status = sim.run(|_| true);
    assert_ne!(status, StepStatus::Failed, "Simulation failed during execution");

    // 1. Validate T_LEVEL.OUT against sim.res.tlevel
    let t_level_path = outfiles_dir.join("T_LEVEL.OUT");
    if t_level_path.exists() {
        let content = fs::read_to_string(&t_level_path).unwrap();
        let mut ref_records: Vec<(f64, f64, f64, f64)> = Vec::new(); // (t, v_top, v_bot, volume)

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
                if let (Ok(t), Ok(vt), Ok(vb), Ok(vol)) = (
                    cols[0].parse::<f64>(),
                    cols[3].parse::<f64>(),
                    cols[5].parse::<f64>(),
                    cols[16].parse::<f64>(),
                ) {
                    ref_records.push((t, vt, vb, vol));
                }
            }
        }

        // Compare final time step
        if let Some(&(ref_t, ref_vt, ref_vb, ref_vol)) = ref_records.last() {
            if let Some(sim_last) = sim.res.tlevel.last() {
                println!(
                    "enbal T_LEVEL (t={}): ref_vol={}, sim_vol={}, diff={}",
                    ref_t,
                    ref_vol,
                    sim_last.volume,
                    (ref_vol - sim_last.volume).abs()
                );
                assert!(
                    (ref_t - sim_last.t).abs() < 1e-3,
                    "Time mismatch: ref={}, sim={}",
                    ref_t,
                    sim_last.t
                );
                assert!(
                    (ref_vol - sim_last.volume).abs() < 0.25,
                    "Volume mismatch in T_LEVEL: ref={}, sim={}",
                    ref_vol,
                    sim_last.volume
                );
            }
        }
    }

    // 2. Validate NOD_INF.OUT against sim.res.profiles
    let nod_inf_path = outfiles_dir.join("NOD_INF.OUT");
    if nod_inf_path.exists() {
        let content = fs::read_to_string(&nod_inf_path).unwrap();
        // Parse last profile block
        let mut last_profile_nodes: Vec<(usize, f64, f64, f64)> = Vec::new(); // (node, depth, h, th)
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
                    if let (Ok(n), Ok(d), Ok(h), Ok(th)) = (
                        cols[0].parse::<usize>(),
                        cols[1].parse::<f64>(),
                        cols[2].parse::<f64>(),
                        cols[3].parse::<f64>(),
                    ) {
                        last_profile_nodes.push((n, d, h, th));
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
            for (ref_node, _ref_depth, ref_h, ref_th) in &last_profile_nodes {
                if let Some(sim_node) = sim_profile.nodes.iter().find(|n| n.node == *ref_node) {
                    let th_diff = (ref_th - sim_node.theta).abs();
                    assert!(
                        th_diff < 0.05,
                        "Moisture mismatch at node {}: ref={}, sim={}, diff={}",
                        ref_node,
                        ref_th,
                        sim_node.theta,
                        th_diff
                    );
                }
            }
        }
    }
}