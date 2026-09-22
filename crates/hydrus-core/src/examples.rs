//! Built-in example projects.

use crate::model::*;

fn uniform_profile(depth: f64, n: usize, h: f64, mat_of: impl Fn(f64) -> usize) -> Profile {
    let dz = depth / (n as f64 - 1.0);
    let nodes = (0..n)
        .map(|i| {
            let z = i as f64 * dz;
            let m = mat_of(z);
            let mut nd = ProfileNode::new(z, h, m);
            nd.layer = m;
            nd
        })
        .collect();
    Profile { nodes, observation_nodes: vec![] }
}

/// Infiltration and evaporation in a two-layer soil with root water uptake.
pub fn infiltration_evaporation() -> Project {
    let mut p = Project::default();
    p.title = "Infiltration and evaporation in a layered soil with root water uptake".into();
    p.processes.root_water_uptake = true;
    p.water.max_iter = 20;
    p.water.materials = vec![
        SoilMaterial { name: "Loam".into(), qr: 0.078, qs: 0.43, alpha: 0.036, n: 1.56, ks: 24.96, l: 0.5, ..Default::default() },
        SoilMaterial { name: "Sand".into(), qr: 0.045, qs: 0.43, alpha: 0.145, n: 2.68, ks: 712.8, l: 0.5, ..Default::default() },
    ];
    p.water.bc = WaterBc { atmospheric: true, top_time_variable: true, surface_layer: false, kod_top: -1, free_drainage: true, kod_bot: -1, ..Default::default() };
    p.profile = uniform_profile(100.0, 101, -200.0, |z| if z < 50.0 { 1 } else { 2 });
    for nd in p.profile.nodes.iter_mut() {
        nd.beta = if nd.depth < 30.0 { 1.0 } else { 0.0 };
        if nd.mat == 2 {
            nd.h = -80.0;
        }
    }
    p.profile.observation_nodes = vec![21, 51, 101];
    p.time = TimeSettings { dt: 0.001, dt_min: 1e-5, dt_max: 5.0, t_max: 30.0, print_times: vec![5.0, 10.0, 15.0, 20.0, 25.0, 30.0], ..Default::default() };
    let rec = |t, prec, evap, tr| AtmRecord { t, prec, evap, transp: tr, h_crit_a: 10000.0, ..Default::default() };
    p.atmosphere.records = vec![rec(2.0, 3.0, 0.3, 0.4), rec(4.0, 0.0, 0.5, 0.5), rec(8.0, 0.0, 0.5, 0.5), rec(12.0, 8.0, 0.3, 0.4), rec(20.0, 0.0, 0.6, 0.6), rec(30.0, 0.0, 0.6, 0.6)];
    p
}

/// Tracer pulse transported through a uniform loam column under steady infiltration,
/// with linear sorption and first-order decay.
pub fn solute_pulse() -> Project {
    let mut p = Project::default();
    p.title = "Tracer pulse with sorption and decay under steady infiltration".into();
    p.processes.solute = true;
    p.water.materials = vec![SoilMaterial { name: "Loam".into(), qr: 0.078, qs: 0.43, alpha: 0.036, n: 1.56, ks: 24.96, l: 0.5, ..Default::default() }];
    p.water.bc = WaterBc { kod_top: -1, r_top: -2.0, free_drainage: true, kod_bot: -1, ..Default::default() };
    p.profile = uniform_profile(100.0, 201, -40.0, |_| 1);
    p.profile.observation_nodes = vec![41, 101, 201];
    p.time = TimeSettings { dt: 0.001, dt_min: 1e-6, dt_max: 0.5, t_max: 60.0, print_times: vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0], ..Default::default() };
    p.solute = SoluteSettings {
        epsi: 0.5,
        artificial_dispersion: true,
        tortuosity: true,
        max_iter: 20,
        tol_abs: 0.001,
        tol_rel: 0.01,
        t_pulse: 10.0,
        k_top: -1,
        k_bot: 0,
        materials: vec![SoluteMaterial { bulk_density: 1.3, disp_l: 2.0, frac: 1.0, th_immobile: 0.0 }],
        species: vec![Species {
            name: "Tracer".into(),
            diff_w: 0.5,
            diff_g: 0.0,
            per_material: vec![SpeciesMaterial { ks: 0.2, beta: 1.0, mu_w: 0.01, mu_s: 0.01, ..Default::default() }],
            c_top: 10.0,
            c_bot: 0.0,
        }],
        ..Default::default()
    };
    // a wet, steady profile
    for nd in p.profile.nodes.iter_mut() {
        nd.h = -40.0;
    }
    p
}

pub fn all() -> Vec<(&'static str, fn() -> Project)> {
    vec![("Infiltration & evaporation with roots", infiltration_evaporation), ("Solute pulse with sorption and decay", solute_pulse)]
}
