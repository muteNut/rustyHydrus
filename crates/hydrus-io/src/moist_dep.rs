//! Parser and serializer for HYDRUS-1D MoistDep.in files.

use hydrus_core::model::{MoistDepTable, Project};
use hydrus_core::HydrusError;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

pub fn read_moist_dep<P: AsRef<Path>>(path: P, prj: &mut Project) -> Result<(), HydrusError> {
    let file = File::open(path).map_err(|e| HydrusError::Io(e.to_string()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines().filter_map(|l| l.ok());

    // Skip main header line
    lines.next();

    let n_mat = prj.water.materials.len();
    let n_sol = prj.n_solutes();
    let mut tables = Vec::with_capacity(n_sol);

    for _ in 0..n_mat {
        lines.next(); // Skip "Material M"
        for _ in 0..n_sol {
            lines.next(); // Skip "Solute j"
            let mut tab = MoistDepTable::default();
            for r in 0..9 {
                if let Some(line) = lines.next() {
                    let vals: Vec<f64> = line
                        .split_whitespace()
                        .filter_map(|s| s.parse().ok())
                        .collect();
                    if vals.len() >= 6 {
                        if r == 0 && tab.theta.is_empty() {
                            tab.theta = vec![vals[1], vals[2], vals[3], vals[4]];
                        }
                        tab.factors[r] = vec![vals[0], 1.0, 1.0, vals[5]];
                    }
                }
            }
            // Skip trailing reaction rows (10..13) if present
            for _ in 9..13 {
                lines.next();
            }
            tables.push(tab);
        }
    }

    prj.solute.moist_tables = tables;
    Ok(())
}

pub fn write_moist_dep<P: AsRef<Path>>(path: P, prj: &Project) -> Result<(), HydrusError> {
    let mut file = File::create(path).map_err(|e| HydrusError::Io(e.to_string()))?;
    writeln!(file, "Heading: Reaction rate dependence on water content")
        .map_err(|e| HydrusError::Io(e.to_string()))?;

    for (m_idx, _) in prj.water.materials.iter().enumerate() {
        writeln!(file, "Material {}", m_idx + 1).map_err(|e| HydrusError::Io(e.to_string()))?;
        for (j_idx, _) in prj.solute.species.iter().enumerate() {
            writeln!(file, "Solute {}", j_idx + 1).map_err(|e| HydrusError::Io(e.to_string()))?;
            let tab = prj.solute.moist_tables.get(j_idx);

            for r in 0..13 {
                if let Some(t) = tab {
                    if r < 9 && t.factors[r].len() >= 4 && t.theta.len() >= 4 {
                        writeln!(
                            file,
                            "{:10.4} {:10.4} {:10.4} {:10.4} {:10.4} {:10.4}",
                            t.factors[r][0], t.theta[0], t.theta[1], t.theta[2], t.theta[3], t.factors[r][3]
                        ).map_err(|e| HydrusError::Io(e.to_string()))?;
                        continue;
                    }
                }
                writeln!(file, "    1.0000     0.0000     0.0000     0.0000     0.0000     1.0000")
                    .map_err(|e| HydrusError::Io(e.to_string()))?;
            }
        }
    }
    Ok(())
}