//! Result export (CSV) — one file per HYDRUS output group.

use hydrus_core::{HydrusError, Results};
use std::fmt::Write as _;
use std::path::Path;

fn w(dir: &Path, name: &str, s: String) -> Result<(), HydrusError> {
    std::fs::write(dir.join(name), s).map_err(|e| HydrusError::Io(format!("{}: {}", name, e)))
}

pub fn tlevel_csv(r: &Results) -> String {
    let mut s = String::new();
    s.push_str("time,dt,tlevel,iter_w,it_cum,kod_top,kod_bot,rTop,rRoot,vTop,vRoot,vBot,sum_rTop,sum_rRoot,sum_vTop,sum_vRoot,sum_vBot,hTop,hRoot,hBot,runoff,sum_runoff,volume,sum_infil,sum_evap,precip");
    for j in 0..r.n_solutes {
        let _ = write!(s, ",cvTop{0},cvBot{0},sum_cvTop{0},sum_cvBot{0},sum_cvCh0_{0},sum_cvCh1_{0},cTop{0},cRoot{0},cBot{0},cvRoot{0},sum_cvRoot{0},sum_cvNEql{0}", j + 1);
    }
    s.push_str(",temp_top,temp_bot,peclet,courant\n");
    for t in &r.tlevel {
        let _ = write!(
            s,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            t.t, t.dt, t.tlevel, t.iter_w, t.it_cum, t.kod_top, t.kod_bot, t.r_top, t.r_root, t.v_top, t.v_root, t.v_bot, t.cum_r_top, t.cum_r_root, t.cum_v_top, t.cum_v_root, t.cum_v_bot,
            t.h_top, t.h_root, t.h_bot, t.run_off, t.cum_run_off, t.volume, t.cum_infil, t.cum_evap, t.precip
        );
        for c in &t.solutes {
            let _ = write!(
                s,
                ",{},{},{},{},{},{},{},{},{},{},{},{}",
                c.cv_top, c.cv_bot, c.cum_top, c.cum_bot, c.cum_ch0, c.cum_ch1, c.c_top, c.c_root, c.c_bot, c.cv_root, c.cum_root, c.cum_neq
            );
        }
        let _ = writeln!(s, ",{},{},{},{}", t.temp_top, t.temp_bot, t.peclet, t.courant);
    }
    s
}

pub fn profiles_csv(r: &Results) -> String {
    let mut s = String::from("time,node,depth,head,theta,K,C,flux,sink,kappa,v_over_Ks,temp,h_matrix,th_matrix");
    for j in 0..r.n_solutes {
        let _ = write!(s, ",conc{0},sorb{0},conc_matrix{0},flux_conc{0}", j + 1);
    }
    s.push('\n');
    for p in &r.profiles {
        for n in &p.nodes {
            let h_m_str = n.h_matrix.map(|v| v.to_string()).unwrap_or_default();
            let th_m_str = n.th_matrix.map(|v| v.to_string()).unwrap_or_default();
            let _ = write!(
                s,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                p.t, n.node, n.depth, n.h, n.theta, n.k, n.c, n.flux, n.sink, n.kappa, n.v_over_ks, n.temp, h_m_str, th_m_str
            );
            for j in 0..r.n_solutes {
                let c_val = n.conc.get(j).copied().unwrap_or(0.0);
                let s_val = n.sorb.get(j).copied().unwrap_or(0.0);
                let cm_val = n.conc_matrix.as_ref().and_then(|m| m.get(j)).map(|v| v.to_string()).unwrap_or_default();
                let fc_val = n.flux_conc.as_ref().and_then(|f| f.get(j)).map(|v| v.to_string()).unwrap_or_default();
                let _ = write!(s, ",{},{},{},{}", c_val, s_val, cm_val, fc_val);
            }
            s.push('\n');
        }
    }
    s
}

pub fn obs_csv(r: &Results) -> String {
    let mut s = String::from("time,node,head,theta,temp,flux,h_matrix,th_matrix");
    for j in 0..r.n_solutes {
        let _ = write!(s, ",conc{0},conc_matrix{0},flux_conc{0}", j + 1);
    }
    s.push('\n');
    for o in &r.obs {
        for p in &o.points {
            let h_m_str = p.h_matrix.map(|v| v.to_string()).unwrap_or_default();
            let th_m_str = p.th_matrix.map(|v| v.to_string()).unwrap_or_default();
            let _ = write!(s, "{},{},{},{},{},{},{},{}", o.t, p.node, p.h, p.theta, p.temp, p.flux, h_m_str, th_m_str);
            for j in 0..r.n_solutes {
                let c_val = p.conc.get(j).copied().unwrap_or(0.0);
                let cm_val = p.conc_matrix.as_ref().and_then(|m| m.get(j)).map(|v| v.to_string()).unwrap_or_default();
                let fc_val = p.flux_conc.as_ref().and_then(|f| f.get(j)).map(|v| v.to_string()).unwrap_or_default();
                let _ = write!(s, ",{},{},{}", c_val, cm_val, fc_val);
            }
            s.push('\n');
        }
    }
    s
}

pub fn balance_csv(r: &Results) -> String {
    let mut s = String::from("time,area,volume,inflow,h_mean,top_flux,bot_flux,WatBalT,WatBalR_percent\n");
    for b in &r.balance {
        let _ = writeln!(s, "{},{},{},{},{},{},{},{},{}", b.t, b.total.area, b.total.volume, b.total.change, b.total.h_mean, b.top_flux, b.bot_flux, b.wat_bal_t, b.wat_bal_r);
    }
    s
}

pub fn alevel_csv(r: &Results) -> String {
    let mut s = String::from("time,sum_rTop,sum_rRoot,sum_vTop,sum_vRoot,sum_vBot,hTop,hRoot,hBot\n");
    for a in &r.alevel {
        let _ = writeln!(s, "{},{},{},{},{},{},{},{},{}", a.t, a.cum[0], a.cum[1], a.cum[2], a.cum[3], a.cum[4], a.h_top, a.h_root, a.h_bot);
    }
    s
}

/// Write all result tables as CSV files into `dir`.
pub fn write_all(r: &Results, dir: &Path) -> Result<(), HydrusError> {
    std::fs::create_dir_all(dir).map_err(|e| HydrusError::Io(e.to_string()))?;
    w(dir, "T_LEVEL.csv", tlevel_csv(r))?;
    w(dir, "NOD_INF.csv", profiles_csv(r))?;
    w(dir, "OBS_NODE.csv", obs_csv(r))?;
    w(dir, "BALANCE.csv", balance_csv(r))?;
    w(dir, "A_LEVEL.csv", alevel_csv(r))?;
    Ok(())
}
