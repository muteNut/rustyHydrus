use hydrus_core::{Simulation, StepStatus};
use hydrus_io::{read_legacy_project, read_nod_inf, sanitize_fortran_line};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

// Strict benchmark tolerances
const REL_TOL: f64 = 0.01; // 1% relative tolerance
const ABS_TOL: f64 = 0.01; // 0.01 cm / flux unit absolute tolerance
const TIME_MATCH_WINDOW: f64 = 0.005; // Tight time-matching window for print records

#[derive(Debug, Default)]
struct TableData {
    headers: Vec<String>,
    records: Vec<Vec<f64>>,
}

fn parse_hydrus_table(path: &Path) -> Option<TableData> {
    let file = File::open(path).ok()?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines().filter_map(|l| l.ok());

    let mut headers = Vec::new();
    let mut records = Vec::new();

    while let Some(line) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('*') {
            continue;
        }

        let upper = trimmed.to_ascii_uppercase();
        if upper.contains("TIME")
            || upper.contains("NODE")
            || upper.contains("VTOP")
            || upper.contains("VBOT")
            || upper.contains("HEAD")
        {
            let cols: Vec<String> = trimmed
                .split_whitespace()
                .map(|s| s.trim_matches(|c: char| !c.is_alphanumeric() && c != '_').to_string())
                .filter(|s| !s.is_empty())
                .collect();

            if !cols.is_empty() && headers.is_empty() {
                headers = cols;
                continue;
            }
        }

        if trimmed.starts_with('[') || trimmed.contains("[L]") || trimmed.contains("[T]") {
            continue;
        }

        let sanitized = sanitize_fortran_line(trimmed);
        let row_vals: Vec<f64> = sanitized
            .split_whitespace()
            .filter_map(|s| s.replace(['d', 'D'], "e").parse::<f64>().ok())
            .collect();

        if !row_vals.is_empty() && !headers.is_empty() && row_vals.len() >= 2 {
            records.push(row_vals);
        }
    }

    if headers.is_empty() || records.is_empty() {
        None
    } else {
        Some(TableData { headers, records })
    }
}

fn compare_series(
    sim_times: &[f64],
    sim_data: &HashMap<String, Vec<f64>>,
    ref_table: &TableData,
    file_name: &str,
    folder_name: &str,
) -> Vec<String> {
    let mut mismatches = Vec::new();
    let time_col_idx = match ref_table
        .headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case("time") || h.eq_ignore_ascii_case("t"))
    {
        Some(idx) => idx,
        None => return mismatches,
    };

    for (col_idx, header) in ref_table.headers.iter().enumerate() {
        if col_idx == time_col_idx {
            continue;
        }
        let norm_header = header.to_ascii_lowercase();

        if let Some(sim_col) = sim_data.get(&norm_header) {
            let mut failed_count = 0;
            let mut checked_count = 0;
            let mut max_diff = 0.0;
            let mut worst_sim = 0.0;
            let mut worst_ref = 0.0;

            for ref_row in &ref_table.records {
                if col_idx >= ref_row.len() {
                    continue;
                }
                let ref_t = ref_row[time_col_idx];
                let ref_val = ref_row[col_idx];

                let mut best_idx = None;
                let mut min_dt = 1e30;
                for (s_idx, &st) in sim_times.iter().enumerate() {
                    let dt = (st - ref_t).abs();
                    if dt < min_dt {
                        min_dt = dt;
                        best_idx = Some(s_idx);
                    }
                }

                if min_dt < TIME_MATCH_WINDOW {
                    if let Some(s_idx) = best_idx {
                        let sim_val = sim_col[s_idx];
                        let diff = (sim_val - ref_val).abs();
                        let denom = ref_val.abs().max(sim_val.abs());
                        let rel_diff = if denom > 1e-12 { diff / denom } else { diff };

                        checked_count += 1;
                        if diff > ABS_TOL && rel_diff > REL_TOL {
                            failed_count += 1;
                            if diff > max_diff {
                                max_diff = diff;
                                worst_sim = sim_val;
                                worst_ref = ref_val;
                            }
                        }
                    }
                }
            }

            if failed_count > 0 {
                mismatches.push(format!(
                    "[{}/{}] Col '{}': {}/{} steps failed (max diff: {:.5}, sim: {:.5}, ref: {:.5})",
                    folder_name, file_name, header, failed_count, checked_count, max_diff, worst_sim, worst_ref
                ));
            }
        }
    }

    mismatches
}

fn run_and_verify_example_strict(folder: &Path) -> Result<(), Vec<String>> {
    let folder_name = folder.file_name().unwrap().to_string_lossy();
    let prj = match read_legacy_project(folder) {
        Ok(p) => p,
        Err(e) => return Err(vec![format!("[{}] Project parse error: {}", folder_name, e)]),
    };

    let mut sim = match Simulation::new(prj) {
        Ok(s) => s,
        Err(e) => return Err(vec![format!("[{}] Init error: {}", folder_name, e)]),
    };

    while sim.step() == StepStatus::Running {}

    let res = &sim.res;
    let mut errors = Vec::new();

    // 1. T_Level.out (all columns)
    let t_level_path = folder.join("T_Level.out");
    if let Some(ref_tlevel) = parse_hydrus_table(&t_level_path) {
        let sim_times: Vec<f64> = res.tlevel.iter().map(|r| r.t).collect();
        let mut sim_tlevel: HashMap<String, Vec<f64>> = HashMap::new();
        sim_tlevel.insert("rtop".into(), res.tlevel.iter().map(|r| r.r_top).collect());
        sim_tlevel.insert("rroot".into(), res.tlevel.iter().map(|r| r.r_root).collect());
        sim_tlevel.insert("vtop".into(), res.tlevel.iter().map(|r| r.v_top).collect());
        sim_tlevel.insert("vroot".into(), res.tlevel.iter().map(|r| r.v_root).collect());
        sim_tlevel.insert("vbot".into(), res.tlevel.iter().map(|r| r.v_bot).collect());
        sim_tlevel.insert("htop".into(), res.tlevel.iter().map(|r| r.h_top).collect());
        sim_tlevel.insert("hroot".into(), res.tlevel.iter().map(|r| r.h_root).collect());
        sim_tlevel.insert("hbot".into(), res.tlevel.iter().map(|r| r.h_bot).collect());
        sim_tlevel.insert("volume".into(), res.tlevel.iter().map(|r| r.volume).collect());
        sim_tlevel.insert("sum_infil".into(), res.tlevel.iter().map(|r| r.cum_infil).collect());
        sim_tlevel.insert("sum_evap".into(), res.tlevel.iter().map(|r| r.cum_evap).collect());

        errors.extend(compare_series(&sim_times, &sim_tlevel, &ref_tlevel, "T_Level.out", &folder_name));
    }

    // 2. Nod_Inf.out (all print times, all columns: h, theta, temp)
    let nod_inf_path = folder.join("Nod_Inf.out");
    if nod_inf_path.exists() {
        if let Ok(ref_blocks) = read_nod_inf(&nod_inf_path) {
            for ref_block in &ref_blocks {
                let matched_profile = res.profiles.iter().find(|p| {
                    (p.t - ref_block.time).abs() < TIME_MATCH_WINDOW
                });

                if let Some(sim_profile) = matched_profile {
                    for ref_node in &ref_block.nodes {
                        if let Some(sim_node) = sim_profile.nodes.iter().find(|n| n.node == ref_node.node) {
                            // Check Head
                            let diff_h = (sim_node.h - ref_node.head).abs();
                            let rel_h = diff_h / ref_node.head.abs().max(sim_node.h.abs()).max(1e-12);
                            if diff_h > ABS_TOL && rel_h > REL_TOL {
                                errors.push(format!(
                                    "[{}/Nod_Inf @ t={:.2}] Node {} 'Head' mismatch: sim={:.4}, ref={:.4}",
                                    folder_name, ref_block.time, ref_node.node, sim_node.h, ref_node.head
                                ));
                            }

                            // Check Theta (tolerance: 0.005 cm3/cm3)
                            let diff_th = (sim_node.theta - ref_node.moisture).abs();
                            if diff_th > 0.005 {
                                errors.push(format!(
                                    "[{}/Nod_Inf @ t={:.2}] Node {} 'Theta' mismatch: sim={:.4}, ref={:.4}",
                                    folder_name, ref_block.time, ref_node.node, sim_node.theta, ref_node.moisture
                                ));
                            }

                            // Check Temperature if heat flow is active
                            if sim.prj.processes.heat {
                                let diff_temp = (sim_node.temp - ref_node.temp).abs();
                                if diff_temp > 0.1 {
                                    errors.push(format!(
                                        "[{}/Nod_Inf @ t={:.2}] Node {} 'Temp' mismatch: sim={:.4}, ref={:.4}",
                                        folder_name, ref_block.time, ref_node.node, sim_node.temp, ref_node.temp
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if errors.is_empty() {
        Ok(())
    } else {
        Err(errors)
    }
}

#[test]
fn test_all_direct_examples_comparison() {
    let base_path = PathBuf::from(r"C:\Users\Public\Documents\PC-Progress\Hydrus-1D 4.xx\Examples\Direct");
    if !base_path.exists() {
        println!("Skipping benchmark tests: Directory {:?} not found", base_path);
        return;
    }

    let entries = std::fs::read_dir(&base_path).expect("Must read Direct examples directory");
    let mut total_run = 0;
    let mut passed = 0;
    let mut failure_reports = Vec::new();

    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() && path.join("SELECTOR.IN").exists() {
            total_run += 1;
            let folder_name = path.file_name().unwrap().to_string_lossy().to_string();
            print!("Testing {:<12} ... ", folder_name);

            match run_and_verify_example_strict(&path) {
                Ok(_) => {
                    println!("PASSED");
                    passed += 1;
                }
                Err(errs) => {
                    println!("FAILED");
                    failure_reports.extend(errs);
                }
            }
        }
    }

    println!("\n==========================================");
    println!("Benchmark Results: {}/{} passed", passed, total_run);
    println!("==========================================");

    if !failure_reports.is_empty() {
        for report in &failure_reports {
            eprintln!("{}", report);
        }
        panic!(
            "{} benchmark examples failed strict verification against reference .out files",
            failure_reports.len()
        );
    }
}