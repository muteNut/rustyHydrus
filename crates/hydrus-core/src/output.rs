//! Result containers produced by a simulation.

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SoluteTLevel {
    pub cv_top: f64,
    pub cv_bot: f64,
    pub cum_top: f64,
    pub cum_bot: f64,
    pub cum_ch0: f64,
    pub cum_ch1: f64,
    pub c_top: f64,
    pub c_root: f64,
    pub c_bot: f64,
    pub cv_root: f64,
    pub cum_root: f64,
    pub cum_neq: f64,
}

/// One row of T_LEVEL.OUT (water flow + boundary information)
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct TLevel {
    pub tlevel: usize,
    pub t: f64,
    pub dt: f64,
    pub iter_w: usize,
    pub iter_c: usize,
    pub it_cum: usize,
    pub kod_top: i32,
    pub kod_bot: i32,
    pub converged: bool,
    pub r_top: f64,
    pub r_root: f64,
    pub v_top: f64,
    pub v_root: f64,
    pub v_bot: f64,
    pub cum_r_top: f64,
    pub cum_r_root: f64,
    pub cum_v_top: f64,
    pub cum_v_root: f64,
    pub cum_v_bot: f64,
    pub h_top: f64,
    pub h_root: f64,
    pub h_bot: f64,
    pub run_off: f64,
    pub cum_run_off: f64,
    pub volume: f64,
    pub cum_infil: f64,
    pub cum_evap: f64,
    pub precip: f64,
    pub peclet: f64,
    pub courant: f64,
    pub temp_top: f64,
    pub temp_bot: f64,
    pub solutes: Vec<SoluteTLevel>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct NodeOut {
    pub node: usize,
    /// depth below the surface (>= 0)
    pub depth: f64,
    pub h: f64,
    pub theta: f64,
    pub k: f64,
    pub c: f64,
    pub flux: f64,
    pub sink: f64,
    pub kappa: i32,
    pub v_over_ks: f64,
    pub temp: f64,
    pub conc: Vec<f64>,
    pub sorb: Vec<f64>,
    pub sorb2: Vec<f64>,
    pub h_matrix: Option<f64>,
    pub th_matrix: Option<f64>,
    pub conc_matrix: Option<Vec<f64>>,
    pub sorb_matrix: Option<Vec<f64>>,
    pub flux_conc: Option<Vec<f64>>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ProfileOut {
    pub t: f64,
    pub nodes: Vec<NodeOut>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ObsPoint {
    pub node: usize,
    pub h: f64,
    pub theta: f64,
    pub temp: f64,
    pub flux: f64,
    pub conc: Vec<f64>,
    pub flux_conc: Option<Vec<f64>>,
    pub h_matrix: Option<f64>,
    pub th_matrix: Option<f64>,
    pub conc_matrix: Option<Vec<f64>>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ObsOut {
    pub t: f64,
    pub points: Vec<ObsPoint>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ALevel {
    pub t: f64,
    pub cum: [f64; 5],
    pub h_top: f64,
    pub h_root: f64,
    pub h_bot: f64,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct SubRegion {
    pub area: f64,
    pub volume: f64,
    pub change: f64,
    pub h_mean: f64,
    pub c_vol: Vec<f64>,
    pub c_mean: Vec<f64>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct BalanceOut {
    pub t: f64,
    pub total: SubRegion,
    pub sub: Vec<SubRegion>,
    pub top_flux: f64,
    pub bot_flux: f64,
    pub wat_bal_t: f64,
    pub wat_bal_r: f64,
    pub sol_bal_t: Vec<f64>,
    pub sol_bal_r: Vec<f64>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Results {
    pub tlevel: Vec<TLevel>,
    pub profiles: Vec<ProfileOut>,
    pub obs: Vec<ObsOut>,
    pub alevel: Vec<ALevel>,
    pub balance: Vec<BalanceOut>,
    pub messages: Vec<String>,
    pub n_solutes: usize,
    pub obs_nodes: Vec<usize>,
    pub finished: bool,
    pub failed: bool,
}
