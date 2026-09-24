//! Simulation driver: state container, initialisation and the main time loop
//! (port of HYDRUS.FOR, INPUT.FOR initialisation and TIME.FOR).

use crate::error::HydrusError;
use crate::heat::HeatState;
use crate::material::*;
use crate::model::*;
use crate::output::*;
use crate::solute::SoluteState;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StepStatus {
    Running,
    Finished,
    Failed,
}

pub struct Simulation {
    pub prj: Project,
    pub res: Results,
    pub n: usize,
    pub n_mat: usize,
    pub model: SoilModel,
    pub i_hyst: i32,
    pub x_conv: f64,
    pub t_conv: f64,
    // process switches
    pub l_wat: bool,
    pub l_chem: bool,
    pub l_temp: bool,
    pub sink_f: bool,
    pub l_root: bool,
    pub short_o: bool,
    pub top_inf: bool,
    pub bot_inf: bool,
    pub atm_bc: bool,
    pub w_layer: bool,
    pub free_d: bool,
    pub seep_f: bool,
    pub gwl_f: bool,
    pub q_drain: bool,
    pub l_var_bc: bool,
    pub max_it: usize,
    pub tol_th: f64,
    pub tol_h: f64,
    // materials
    pub par_d: Vec<Par>,
    pub par_w: Vec<Par>,
    pub ah_w: Vec<f64>,
    pub ath_w: Vec<f64>,
    pub ak_w: Vec<f64>,
    pub tabs: Vec<MatTable>,
    pub l_table: bool,
    pub h_tab_first: f64,
    pub h_tab_last: f64,
    pub h_sat: Vec<f64>,
    pub con_sat: Vec<f64>,
    pub thr: Vec<f64>,
    pub ths: Vec<f64>,
    pub con_s_max: f64,
    // mesh & node properties (index 0 = bottom)
    pub x: Vec<f64>,
    pub mat: Vec<usize>,
    pub lay: Vec<usize>,
    pub beta: Vec<f64>,
    pub ah: Vec<f64>,
    pub ak: Vec<f64>,
    pub ath: Vec<f64>,
    // water state
    pub h_new: Vec<f64>,
    pub h_old: Vec<f64>,
    pub h_temp: Vec<f64>,
    pub th_new: Vec<f64>,
    pub th_old: Vec<f64>,
    pub th_eq: Vec<f64>,
    pub con: Vec<f64>,
    pub cap: Vec<f64>,
    pub sink: Vec<f64>,
    pub con_o: Vec<f64>,
    pub kappa: Vec<i32>,
    pub kappa_o: Vec<i32>,
    pub ath_s: Vec<f64>,
    pub th_rr: Vec<f64>,
    pub con_r: Vec<f64>,
    pub ak_s: Vec<f64>,
    pub v_old: Vec<f64>,
    pub v_new: Vec<f64>,
	// Vapor transport state
	pub l_vapor: bool,
    pub con_lt: Vec<f64>,
    pub con_vt: Vec<f64>,
    pub con_vh: Vec<f64>,
    pub th_v_new: Vec<f64>,
    pub th_v_old: Vec<f64>,
    pub v_v_new: Vec<f64>,
    pub v_v_old: Vec<f64>,
    // boundary conditions
    pub kod_top: i32,
    pub kod_bot: i32,
    pub k_top_old: i32,
    pub k_bot_old: i32,
    pub r_top: f64,
    pub r_bot: f64,
    pub r_root: f64,
    pub h_top: f64,
    pub h_bot: f64,
    pub h_crit_a: f64,
    pub h_crit_s: f64,
    pub gwl0l: f64,
    pub aqh: f64,
    pub bqh: f64,
    pub prec: f64,
    pub r_soil: f64,
    pub cos_alf: f64,
    pub v_top: f64,
    pub v_bot: f64,
    pub h_seep: f64,
    // time control
    pub t: f64,
    pub t_old: f64,
    pub dt: f64,
    pub dt_old: f64,
    pub dt_opt: f64,
    pub dt_min: f64,
    pub dt_max: f64,
    pub d_mul: f64,
    pub d_mul2: f64,
    pub dt_init: f64,
    pub it_min: usize,
    pub it_max: usize,
    pub t_init: f64,
    pub t_max: f64,
    pub t_atm: f64,
    pub t_atm1: f64,
    pub t_atm2: f64,
    pub t_atm_old: f64,
    pub t_print: Vec<f64>,
    pub t_print1: f64,
    pub print_step: usize,
    pub print_daily: bool,
    pub print_int: f64,
    pub p_level: usize,
    pub t_level: usize,
    pub min_step: bool,
    pub iter_w: usize,
    pub iter_c: usize,
    pub it_cum: usize,
    pub convg: bool,
    pub n_noconv: usize,
    pub dt_max_c: f64,
    pub dt_max_t: f64,
    // atmosphere
    pub atm_idx: usize,
	pub meteo_idx: usize,
    pub r_root_d: f64,
    pub r_soil_d: f64,
    pub prec_d: f64,
	pub snow_layer: f64,
    pub interc_state: crate::meteo::InterceptionState,
    // cumulative water quantities
    pub cum_q: [f64; 12],
    pub w_cum_t: f64,
    pub w_cum_a: f64,
    pub w_vol_i: f64,
    pub wat_in: Vec<f64>,
    pub h_root: f64,
    pub v_root: f64,
    pub x_root: f64,
    pub rg: f64,
    pub omega_c: f64,
    pub l_end: bool,
    pub done: bool,
    pub failed: bool,
	// Dual porosity state
    pub i_dual_por: i32,
    pub th_new_im: Vec<f64>,
    pub th_old_im: Vec<f64>,
    pub sink_im: Vec<f64>,
    pub w_transf: f64,
	// Dual permeability state (Model 8)
    pub l_dual_perm: bool,
    pub w_fracture: f64,
    pub alpha_dw: f64,
    pub h_matrix_new: Vec<f64>,
    pub h_matrix_old: Vec<f64>,
    pub th_matrix_new: Vec<f64>,
    pub th_matrix_old: Vec<f64>,
    pub con_matrix: Vec<f64>,
    pub cap_matrix: Vec<f64>,
    // sub-models
    pub heat: Option<HeatState>,
    pub sol: Option<SoluteState>,
    // Bob Lenhard air-entrapment hysteresis (iHyst = 3)
    pub lenhard_i_hyst: i32,
    pub lenhard: Option<crate::lenhard::LenhardState>,
}

impl Simulation {
    pub fn new(prj: Project) -> Result<Simulation, HydrusError> {
        let problems = prj.validate();
        if !problems.is_empty() {
            return Err(HydrusError::Invalid(problems.join("\n")));
        }
        let n = prj.profile.nodes.len();
        let n_mat = prj.water.materials.len();
        let model = prj.water.model;
		let i_hyst = prj.water.hysteresis.code();
		let i_dual_por = match model {
			SoilModel::DualPorosityW => 1,
			SoilModel::DualPorosityH => 2,
			_ => 0,
		};

        let x_conv = prj.units.x_conv();
        let t_conv = prj.units.t_conv();
        let bc = prj.water.bc.clone();
        let l_wat = prj.processes.water_flow;
        let l_chem = prj.processes.solute;
        let l_temp = prj.processes.heat;
        let sink_f = prj.processes.root_water_uptake;
        let l_root = prj.processes.root_growth;
        let ns = prj.n_solutes();

        // ---- materials
        let mut par_d = vec![];
        let mut par_w = vec![[0.0; 11]; n_mat];
        let mut ah_w = vec![1.0; n_mat];
        let mut ath_w = vec![1.0; n_mat];
        let mut ak_w = vec![1.0; n_mat];
        for (m, mat) in prj.water.materials.iter().enumerate() {
            let mut p = par_of(mat, model, x_conv);
            if i_hyst > 0 {
                p[6] = mat.qm.max(mat.qs);
                let mut w = [0.0; 11];
                w[1] = mat.qs_w;
                w[2] = mat.alpha_w;
                w[4] = mat.ks_w;
                w[0] = p[0];
                w[3] = p[3];
                ah_w[m] = p[2] / w[2];
                ath_w[m] = (w[1] - w[0]) / (p[1] - p[0]);
                ak_w[m] = 1.0;
                if i_hyst == 2 {
                    ak_w[m] = w[4] / p[4];
                }
                w[6] = w[0] + ath_w[m] * (p[6] - p[0]);
                p[7] = p[0];
                p[8] = p[1];
                p[9] = p[4];
                w[7] = w[0];
                w[8] = w[1];
                w[9] = w[4];
                w[5] = p[5];
                par_w[m] = w;
            }
            par_d.push(p);
        }
        // tables
        let mut h1 = -prj.water.h_tab1.abs().min(prj.water.h_tab_n.abs());
        let mut hn = -prj.water.h_tab1.abs().max(prj.water.h_tab_n.abs());
        let mut l_table = true;
        if (h1 - hn).abs() < 1e-5 {
            l_table = false;
            h1 = -0.0001 * x_conv;
            hn = -100.0 * x_conv;
        }
        let tabs = if let Some(ref external_tabs) = prj.water.tabs {
            external_tabs.clone()
        } else {
            par_d.iter().map(|p| MatTable::build(model, p, h1, hn)).collect()
        };

        let h_sat: Vec<f64> = par_d.iter().map(|p| fh(model, 1.0, p)).collect();
        let con_sat: Vec<f64> = par_d.iter().map(|p| p[4]).collect();
        let thr: Vec<f64> = par_d.iter().map(|p| p[0]).collect();
        let ths: Vec<f64> = par_d.iter().map(|p| p[1]).collect();
        let con_s_max = con_sat.iter().cloned().fold(0.0, f64::max);

        // ---- profile (flip to bottom -> top)
        let nodes: Vec<&ProfileNode> = prj.profile.nodes.iter().rev().collect();
        let x: Vec<f64> = nodes.iter().map(|nd| -nd.depth).collect();
        let mat: Vec<usize> = nodes.iter().map(|nd| nd.mat - 1).collect();
        let lay: Vec<usize> = nodes.iter().map(|nd| nd.layer.max(1)).collect();
        let mut beta: Vec<f64> = nodes.iter().map(|nd| nd.beta).collect();
        let ah: Vec<f64> = nodes.iter().map(|nd| nd.ah).collect();
        let ak: Vec<f64> = nodes.iter().map(|nd| nd.ak).collect();
        let ath: Vec<f64> = nodes.iter().map(|nd| nd.ath).collect();
        // normalise root distribution (as in NodInf)
        {
            let mut sb = 0.0;
            if beta[n - 1] > 0.0 {
                sb = beta[n - 1] * (x[n - 1] - x[n - 2]) / 2.0;
            }
            for i in 1..n - 1 {
                if beta[i] > 0.0 {
                    sb += beta[i] * (x[i + 1] - x[i - 1]) / 2.0;
                }
            }
            for b in beta.iter_mut().skip(1) {
                *b = if sb > 0.0 { *b / sb } else { 0.0 };
            }
        }
		let mut h0: Vec<f64> = nodes.iter().map(|nd| nd.h).collect();
		let kappa0 = prj.water.hysteresis.initial_branch_code(prj.water.init_kappa);
		let kappa = vec![kappa0; n];
		
		// Lower bound on the pressure head (only meaningful when h0 really is a head)
		if !prj.water.init_in_water_content {
			for i in 0..n {
				let m = mat[i];
				if model.code() < 10 {
					h0[i] = h0[i].max(ah[i] * fh(model, 0.00000001, &par_d[m]));
				}
			}
		}
        // Immobile water content from the initial condition (Fortran InitW, iDualPor > 0)
		let mut th_im_water = vec![0.0; n];

		if prj.water.init_in_water_content {
			for i in 0..n {
				let m = mat[i];
				// Values greater than 1.0 are positive pressure heads (saturated zone/ponding)
				if h0[i] > 1.0 {
					continue;
				}
				if i_dual_por > 0 {
					let th_total = h0[i];
					h0[i] = th_total * par_d[m][1] / (par_d[m][1] + par_d[m][7]);
					th_im_water[i] = th_total - h0[i];
				}
				let (qr, qs) = if kappa[i] == -1 {
					(par_d[m][0], par_d[m][1])
				} else {
					(par_w[m][0], par_w[m][1])
				};
				let d_xz = if ath[i] > 0.0 { ath[i] } else { 1.0 };
				let qe = ((h0[i] - qr) / (qs - qr) / d_xz).min(1.0);
				if qe <= 0.0 {
					return Err(HydrusError::Invalid("Initial water content is lower than θr.".into()));
				}
				if qe >= 1.0 {
					h0[i] = 0.0;
				} else if kappa[i] == -1 {
					h0[i] = fh(model, qe, &par_d[m]) * ah[i];
				} else {
					h0[i] = fh(model, qe, &par_w[m]) * ah[i] * ah_w[m];
				}
			}
		}
        // ---- boundary condition codes (BasInf "input modifications")
        let mut kod_top = bc.kod_top;
        let mut kod_bot = bc.kod_bot;
        let l_var_bc = kod_top == 0;
        if bc.top_time_variable {
            kod_top = if kod_top >= 0 { 3 } else { -3 };
        }
        if bc.bot_time_variable {
            kod_bot = if kod_bot >= 0 { 3 } else { -3 };
        }
        let h_crit_s = prj.atmosphere.h_crit_s;
        if bc.atmospheric && kod_top < 0 {
            kod_top = -4;
        }
        if bc.surface_layer {
            kod_top = -kod_top.abs();
        }
        if bc.gwl_flux {
            kod_bot = -7;
        }
        if bc.free_drainage {
            kod_bot = -5;
        }
        if bc.seepage_face {
            kod_bot = -2;
        }
        let const_flux_read = (!bc.top_time_variable && bc.kod_top == -1)
            || (!bc.bot_time_variable && bc.kod_bot == -1 && !bc.gwl_flux && !bc.free_drainage && !bc.seepage_face && bc.drains.is_none());
        let (r_top, r_bot, r_root) = if const_flux_read { (bc.r_top, bc.r_bot, bc.r_root.abs()) } else { (0.0, 0.0, 0.0) };
		
		// Initialize th_old_im and th_new_im with ths_im (par_d[m][7]) or initial condition
		let mut th_init_im = vec![0.0; n];
		for i in 0..n {
			let m = mat[i];
			th_init_im[i] = if prj.water.init_in_water_content && i_dual_por > 0 {
				th_im_water[i]
			} else {
				par_d[m][7] // overwritten in initialise() when the IC is a pressure head
			};
		}
		
		let l_dual_perm = model == SoilModel::DualPermeability;
		let w_fracture = if l_dual_perm { 0.1 } else { 0.0 };
		let alpha_dw = if l_dual_perm { 1e-4 } else { 0.0 };

        let tm = &prj.time;
        let mut t_print = tm.print_times.clone();
        t_print.push(tm.t_max);
		let th_rr: Vec<f64> = (0..n).map(|i| thr[mat[i]]).collect();

        let mut sim = Simulation {
            n,
            n_mat,
            model,
            i_hyst,
            x_conv,
            t_conv,
            l_wat,
            l_chem,
            l_temp,
            sink_f,
            l_root,
            short_o: prj.processes.short_output,
            top_inf: bc.top_time_variable,
            bot_inf: bc.bot_time_variable,
            atm_bc: bc.atmospheric,
            w_layer: bc.surface_layer,
            free_d: bc.free_drainage,
            seep_f: bc.seepage_face,
            gwl_f: bc.gwl_flux,
            q_drain: bc.drains.is_some(),
            l_var_bc,
            max_it: prj.water.max_iter,
            tol_th: prj.water.tol_th,
            tol_h: prj.water.tol_h,
            par_d,
            par_w,
            ah_w,
            ath_w,
            ak_w,
            tabs,
            l_table,
            h_tab_first: h1,
            h_tab_last: hn,
            h_sat,
            con_sat,
            thr,
            ths,
            con_s_max,
            x,
            mat,
            lay,
            beta,
            ah,
            ak,
            ath,
            h_new: h0.clone(),
            h_old: h0.clone(),
            h_temp: h0.clone(),
            th_new: vec![0.0; n],
            th_old: vec![0.0; n],
            th_eq: vec![0.0; n],
            con: vec![0.0; n],
            cap: vec![0.0; n],
            sink: vec![0.0; n],
            con_o: vec![0.0; n],
            kappa: kappa.clone(),
            kappa_o: kappa.clone(),
            ath_s: vec![1.0; n],
            th_rr,
            con_r: vec![0.0; n],
            ak_s: vec![1.0; n],
            v_old: vec![0.0; n],
            v_new: vec![0.0; n],
            kod_top,
            kod_bot,
            k_top_old: kod_top,
            k_bot_old: kod_bot,
            r_top,
            r_bot,
            r_root,
            h_top: h0[n - 1],
            h_bot: h0[0],
            h_crit_a: -1e10,
            h_crit_s,
            gwl0l: bc.gwl0l,
            aqh: bc.aqh,
            bqh: bc.bqh,
            prec: 0.0,
            r_soil: 0.0,
            cos_alf: prj.cos_alpha,
            v_top: 0.0,
            v_bot: 0.0,
            h_seep: bc.h_seep,
            t: tm.t_init + tm.dt,
            t_old: tm.t_init,
            dt: tm.dt,
            dt_old: tm.dt,
            dt_opt: tm.dt,
            dt_min: tm.dt_min,
            dt_max: tm.dt_max,
            d_mul: tm.d_mul,
            d_mul2: tm.d_mul2,
            dt_init: tm.dt,
            it_min: tm.it_min,
            it_max: tm.it_max,
            t_init: tm.t_init,
            t_max: tm.t_max,
            t_atm: tm.t_max,
            t_atm1: tm.t_max,
            t_atm2: tm.t_max,
            t_atm_old: tm.t_init,
            t_print,
            t_print1: tm.t_max,
            print_step: tm.print_step.max(1),
            print_daily: tm.print_at_interval,
            print_int: tm.print_interval,
            p_level: 0,
            t_level: 1,
            min_step: true,
            iter_w: 0,
            iter_c: 0,
            it_cum: 0,
            convg: true,
            n_noconv: 0,
            dt_max_c: 1e30,
            dt_max_t: 1e30,
            atm_idx: 0,
            meteo_idx: 0,
            r_root_d: 0.0,
            r_soil_d: 0.0,
            prec_d: 0.0,
            snow_layer: 0.0,
            interc_state: crate::meteo::InterceptionState::default(),
            cum_q: [0.0; 12],
            w_cum_t: 0.0,
            w_cum_a: 0.0,
            w_vol_i: 0.0,
            wat_in: vec![0.0; n],
            h_root: 0.0,
            v_root: 0.0,
            x_root: 0.0,
            rg: 0.0,
            omega_c: prj.root.omega_c,
            l_end: false,
            done: false,
            failed: false,
            lenhard_i_hyst: 0,
            lenhard: None,
            l_vapor: prj.processes.vapor,
            con_lt: vec![0.0; n],
            con_vt: vec![0.0; n],
            con_vh: vec![0.0; n],
            th_v_new: vec![0.0; n],
            th_v_old: vec![0.0; n],
            v_v_new: vec![0.0; n],
            v_v_old: vec![0.0; n],
			i_dual_por,
            th_new_im: th_init_im.clone(),
            th_old_im: th_init_im,
            sink_im: vec![0.0; n],
            w_transf: 0.0,
            l_dual_perm,
            w_fracture,
            alpha_dw,
            h_matrix_new: h0.clone(),
            h_matrix_old: h0.clone(),
            th_matrix_new: vec![0.0; n],
            th_matrix_old: vec![0.0; n],
            con_matrix: vec![0.0; n],
            cap_matrix: vec![0.0; n],
            heat: None,
            sol: None,
            res: Results { n_solutes: ns, ..Default::default() },
            prj,
        };
        if sim.print_daily {
            sim.t_print1 = sim.t_init + sim.print_int;
        }
        sim.res.obs_nodes = sim.prj.profile.observation_nodes.clone();
        sim.initialise()?;
        Ok(sim)
    }

	fn initialise(&mut self) -> Result<(), HydrusError> {
		// hTop / hBot come from the initial condition at the boundary nodes
		self.h_top = self.h_new[self.n - 1];
		self.h_bot = self.h_new[0];

		// Temperature first: set_mat (vapor flow) reads it
		if self.l_temp {
			self.heat_init()?;
		}

		self.init_water_state();

		if self.l_chem {
			self.solute_init()?; // uses th_old
		}

		if self.top_inf || self.bot_inf || self.atm_bc {
			self.init_atmosphere()?;
		}

		if self.l_root {
			self.set_rg();
		}
		if self.sink_f {
			self.set_snk();
		}

		// Initial output
		self.profile_out(self.t_init);
		self.sub_reg(0);

		if self.l_chem || self.l_temp || self.l_vapor {
			let (ho, thn) = (self.h_old.clone(), self.th_old.clone());
			let temps: Vec<f64> = (0..self.n).map(|i| self.temp(i)).collect();
			let (v, vv) = self.veloc(&ho, &thn, &thn, &temps);
			self.v_old = v;
			self.v_v_old = vv;
		}
		self.v_new = self.v_old.clone();
		self.v_v_new = self.v_v_old.clone();
		Ok(())
	}

	/// Initial hydraulic properties, water contents and immobile water.
	fn init_water_state(&mut self) {
		if self.i_hyst == 3 {
			let ik = self.prj.water.init_kappa;
			self.lenhard_hyst(ik, 1, crate::lenhard::ThetaTarget::Old);
		} else {
			self.set_mat(0);
			if self.prj.water.init_in_water_content {
				for i in 0..self.n {
					let h_init = self.prj.profile.nodes[self.n - 1 - i].h;
					if h_init >= 0.0 && h_init <= 1.0 {
						self.th_eq[i] = h_init;
					}
				}
			}
			self.th_old.copy_from_slice(&self.th_eq);
		}
		self.th_new.copy_from_slice(&self.th_old);

		if self.l_vapor {
			self.th_v_old = self.th_v_new.clone();
			// th_v_new stays equal to th_v_old at the start (as in HYDRUS.FOR)
		}

		// InitDualPor: only if the IC was a pressure head (otherwise set in new())
		if self.i_dual_por > 0 && !self.prj.water.init_in_water_content {
			for i in 0..self.n {
				let p = self.par_d[self.mat[i]];
				self.th_new_im[i] = if self.i_dual_por == 1 {
					let se = (self.th_old[i] - p[0]) / (p[1] - p[0]);
					p[6] + se * (p[7] - p[6])
				} else {
					let par_im: Par = [p[6], p[7], p[8], p[9], p[10], p[5], 0.0, 0.0, 0.0, 0.0, 0.0];
					fq(SoilModel::VanGenuchten, self.h_new[i], &par_im)
				};
			}
			self.th_old_im = self.th_new_im.clone();
		}
	}

	/// First atmospheric record and derived surface fluxes.
	fn init_atmosphere(&mut self) -> Result<(), HydrusError> {
		self.atm_idx = 0;
		self.meteo_idx = 0;
		self.snow_layer = 0.0;
		self.interc_state = crate::meteo::InterceptionState::default();
		self.t_atm2 = self.t_max;
		self.set_bc()?; // also applies snow / interception (see #10)

		let mut next = self.t_atm1.min(self.t_atm2);
		if let Some(mp) = &self.prj.atmosphere.meteo {
			if let Some(rec) = mp.records.get(self.meteo_idx) {
				if rec.t > self.t_init {
					next = next.min(rec.t);
				}
			}
		}
		self.t_atm = next;

		if self.l_chem {
			self.set_chem_bc();
		}
		if self.prj.atmosphere.daily_variation {
			self.r_root_d = self.r_root;
			self.r_soil_d = self.r_soil;
			self.daily_var_root_soil();
		}
		if self.prj.atmosphere.sinusoidal_precip {
			self.prec_d = self.prec;
			self.t_atm_old = self.t_init;
			self.sin_prec();
		}
		if self.kod_top == -4 || !self.l_var_bc {
			self.r_top = self.r_soil.abs() - self.prec.abs();
		}
		Ok(())
	}
    /// Evaluates dynamic snow accumulation/melt and canopy interception/LAI partitioning.
    pub fn apply_snow_and_interception(&mut self) {
        if !self.prj.atmosphere.snow && !self.prj.atmosphere.interception && !self.prj.atmosphere.lai_partitioning {
            return;
        }

        let temp_air = if let Some(ref mp) = self.prj.atmosphere.meteo {
            if self.meteo_idx < mp.records.len() {
                let r = &mp.records[self.meteo_idx];
                (r.t_max + r.t_min) / 2.0
            } else {
                0.0
            }
        } else if self.atm_idx < self.prj.atmosphere.records.len() {
            self.prj.atmosphere.records[self.atm_idx].t_top
        } else {
            0.0
        };

        // 1. Snow accumulation and melt
        if self.prj.atmosphere.snow {
            let mut c_top_dummy = vec![0.0; self.n_species()];
            let c_t_dummy = vec![0.0; self.n_species()];

            let (p_eff, s_pack, evap_eff, min_step) = crate::meteo::calculate_snow(
                self.prec,
                temp_air,
                self.dt,
                self.prj.atmosphere.snow_mf,
                self.snow_layer,
                self.r_soil,
                self.x_conv,
                &mut c_top_dummy,
                &c_t_dummy,
            );

            self.prec = p_eff;
            self.snow_layer = s_pack;
            self.r_soil = evap_eff;
            if min_step {
                self.min_step = true;
            }
        }

        // 2. Canopy interception & LAI partitioning
        if self.prj.atmosphere.lai_partitioning || self.prj.atmosphere.interception {
            let lai = if let Some(ref mp) = self.prj.atmosphere.meteo {
                if self.meteo_idx < mp.records.len() {
                    mp.records[self.meteo_idx].lai.unwrap_or(mp.lai)
                } else {
                    mp.lai
                }
            } else {
                0.0
            };

            let scf = (1.0 - (-self.prj.atmosphere.extinction.max(0.1) * lai).exp()).max(0.0);

            if self.prj.atmosphere.lai_partitioning && lai > 0.0 {
                let r_pet = self.r_soil;
                self.r_root = r_pet * scf;
                self.r_soil = r_pet - self.r_root;
            }

            if self.prj.atmosphere.interception {
                let (net_p, net_tr, _act_int) = crate::meteo::calculate_interception(
                    self.prec,
                    self.r_root,
                    lai,
                    self.prj.atmosphere.interception_a,
                    scf,
                    &mut self.interc_state,
                );
                self.prec = net_p;
                self.r_root = net_tr;
            }
        }

        // Recompute net surface flux ONLY when an atmospheric BC is active.
        // For non-atmospheric simulations (e.g. constant boundary flux), r_top must NOT be touched.
        if self.atm_bc || self.top_inf {
            if self.kod_top.abs() == 4 || !self.l_var_bc {
                self.r_top = self.r_soil.abs() - self.prec.abs();
            }
        }
    }

    // ---------------------------------------------------------------- main loop
    /// Advance the simulation by one time step.
    pub fn step(&mut self) -> StepStatus {
        if self.done {
            return if self.failed { StepStatus::Failed } else { StepStatus::Finished };
        }
        // ---- water flow
        if self.l_wat {
            self.wat_flow();
            if !self.convg {
                self.n_noconv += 1;
            } else {
                self.n_noconv = 0;
            }
            if !self.convg && self.n_noconv >= 10 {
                self.res.messages.push(format!("The numerical solution has not converged (t = {:.6}).", self.t));
                self.done = true;
                self.failed = true;
                self.res.failed = true;
                return StepStatus::Failed;
            }
        } else {
            self.iter_w = 1;
            self.it_cum += 1;
        }

        if self.l_wat && (self.l_temp || self.l_chem || self.l_vapor) {
            let (hn, thn, tho) = (self.h_new.clone(), self.th_new.clone(), self.th_old.clone());
            let temps: Vec<f64> = (0..self.n).map(|i| self.temp(i)).collect();
            let (v, vv) = self.veloc(&hn, &thn, &tho, &temps);
            self.v_new = v;
            self.v_v_new = vv;
        }
        if self.l_root {
            self.set_rg();
        }
        if self.sink_f {
            self.set_snk();
        }
        if self.l_temp {
            self.temper();
        }
        if self.l_chem {
            if let Err(e) = self.solute() {
                self.res.messages.push(e.to_string());
                self.done = true;
                self.failed = true;
                self.res.failed = true;
                return StepStatus::Failed;
            }
            let _ = self.solute_cumulate(0.0);
        }

        // ---- output
        if (self.t - self.t_max).abs() <= 0.5 * self.dt_min || self.t > self.t_max {
            self.l_end = true;
        }
        let j_print = self.j_print();
        self.tl_inf(j_print);
        if j_print {
            self.obs_nod();
        }
        if self.print_daily && (self.t_print1 - self.t).abs() < 0.001 * self.dt {
            self.t_print1 += self.print_int;
        }
        let at_profile_print = (self.p_level < self.t_print.len() && (self.t_print[self.p_level] - self.t).abs() <= 0.01 * self.dt)
            || self.l_end;

        if at_profile_print {
            self.profile_out(self.t);
            self.sub_reg(self.p_level + 1);
            if self.p_level < self.t_print.len() && (self.t_print[self.p_level] - self.t).abs() <= 0.01 * self.dt {
                self.p_level += 1;
            }
        }

        // ---- A-level
        if (self.t - self.t_atm).abs() <= 0.001 * self.dt && (self.top_inf || self.bot_inf || self.atm_bc) {
            self.a_level();
            if (self.t - self.t_atm1).abs() <= 0.001 * self.dt {
                self.t_atm_old = self.t_atm1;
                if self.set_bc().is_err() {
                    self.res.messages.push("Error reading atmospheric data.".into());
                    self.done = true;
                    self.failed = true;
                    return StepStatus::Failed;
                }
            }
            let mut next_atm_time = self.t_atm1.min(self.t_atm2);
            if let Some(ref mp) = self.prj.atmosphere.meteo {
                if self.meteo_idx < mp.records.len() {
                    let t_meteo = mp.records[self.meteo_idx].t;
                    if t_meteo > self.t {
                        next_atm_time = next_atm_time.min(t_meteo);
                    }
                }
            }
            self.t_atm = next_atm_time;
            if self.l_chem {
                self.set_chem_bc();
            }
            if self.prj.atmosphere.daily_variation {
                self.r_root_d = self.r_root;
                self.r_soil_d = self.r_soil;
            }
            if self.prj.atmosphere.sinusoidal_precip {
                self.prec_d = self.prec;
            }
            self.apply_snow_and_interception();
            if !self.l_var_bc {
                self.r_top = self.r_soil.abs() - self.prec.abs();
            }
        }
        if self.l_chem {
            self.sol_wlayer_mix();
        }

        // ---- time governing
        if (self.t - self.t_max).abs() <= 0.5 * self.dt_min || self.t > self.t_max {
            self.done = true;
            self.res.finished = true;
            return StepStatus::Finished;
        }
        self.t_old = self.t;
        self.dt_old = self.dt;
        self.k_top_old = self.kod_top;
        self.k_bot_old = self.kod_bot;
        if !self.l_wat {
            self.iter_w = 1;
        }
        let iter = self.iter_w.max(self.iter_c);
        let mut dt_max_a = self.dt_max_c.min(self.dt_max_t);
        if self.prj.atmosphere.sinusoidal_precip && self.prec_d > 0.0 {
            dt_max_a = dt_max_a.min((self.t_atm1 - self.t_atm_old) / 20.0);
        }
        let tp = if self.p_level < self.t_print.len() { self.t_print[self.p_level] } else { self.t_max };
        self.tm_cont(iter, tp.min(self.t_print1), dt_max_a);
        self.t += self.dt;

        // Apply sub-daily variations (only modify r_top if atmospheric BC is enabled)
        if self.prj.atmosphere.daily_variation {
            self.daily_var_root_soil();
            if self.atm_bc || self.top_inf {
                self.r_top = self.r_soil.abs() - self.prec.abs();
            }
        }
        if self.prj.atmosphere.sinusoidal_precip {
            self.sin_prec();
            if self.atm_bc || self.top_inf {
                self.r_top = self.r_soil.abs() - self.prec.abs();
            }
        }

        self.t_level += 1;
        if self.t_level > 999_999 {
            self.t_level = 2;
        }
        self.update();
        StepStatus::Running
    }

    /// Run to completion, calling `progress(t_now, t_max)` after each step.
    /// The callback returns `false` to cancel.
    pub fn run<F: FnMut(&Simulation) -> bool>(&mut self, mut progress: F) -> StepStatus {
        loop {
            let s = self.step();
            if s != StepStatus::Running {
                return s;
            }
            if !progress(self) {
                return StepStatus::Running;
            }
        }
    }

    fn j_print(&self) -> bool {
        let at_print = self.p_level < self.t_print.len() && (self.t_print[self.p_level] - self.t).abs() < 0.001 * self.dt;
        let at_daily = self.print_daily && (self.t_print1 - self.t).abs() < 0.001 * self.dt;
        let regular = !self.short_o && (self.t_level - 1) % self.print_step == 0;
        at_print || at_daily || regular
    }

    /// Update variables at the end of a time step (Update in Fortran).
    fn update(&mut self) {
        let n = self.n;
        let i_bot = if self.kod_bot > 0 { 1 } else { 0 };
        let i_top = if self.kod_top > 0 { n - 2 } else { n - 1 };
        let mut l_sat = true;
		if self.i_dual_por > 0 {
			self.th_old_im.copy_from_slice(&self.th_new_im);
		}
		if self.l_dual_perm {
            self.h_matrix_old.copy_from_slice(&self.h_matrix_new);
            self.th_matrix_old.copy_from_slice(&self.th_matrix_new);
        }
        if self.l_wat {
            for i in i_bot..=i_top {
                if self.h_new[i] < 0.0 && self.h_old[i] < 0.0 {
                    self.h_temp[i] = self.h_new[i] + (self.h_new[i] - self.h_old[i]) * self.dt / self.dt_old;
                } else {
                    self.h_temp[i] = self.h_new[i];
                }
                self.h_old[i] = self.h_new[i];
                self.h_new[i] = self.h_temp[i];
            }
            if self.kod_top > 0 {
                self.h_old[n - 1] = self.h_new[n - 1];
            }
            if self.kod_bot > 0 {
                self.h_old[0] = self.h_new[0];
            }
        }
        for i in 0..n {
            if self.l_wat {
                self.th_old[i] = self.th_new[i];
                if self.l_temp || self.l_chem {
                    self.v_old[i] = self.v_new[i];
                }
                if self.h_new[i] < 0.0 {
                    l_sat = false;
                }
            }
            if self.l_temp {
                if let Some(h) = self.heat.as_mut() {
                    h.temp_o[i] = h.temp_n[i];
                }
            }
        }
        if self.l_wat && l_sat && (self.kod_top == -1 || self.kod_top == 4) && self.r_top >= -self.con_s_max {
            if self.kod_top == 4 {
                self.kod_top = -4;
            }
            let v = -0.005 * self.x_conv;
            self.h_new[n - 1] = v;
            self.h_temp[n - 1] = v;
            self.h_old[n - 1] = v;
        }
        if self.l_vapor {
            self.th_v_old.copy_from_slice(&self.th_v_new);
            self.v_v_old.copy_from_slice(&self.v_v_new);
        }
    }

    // ------------------------------------------------------------- utilities
    pub fn n_species(&self) -> usize {
        self.prj.n_solutes()
    }
    pub fn depth_of(&self, i: usize) -> f64 {
        self.x[self.n - 1] - self.x[i]
    }
	pub fn theta_of(&self, m: usize, h: f64) -> f64 {
        if m < self.par_d.len() {
            crate::material::fq(self.model, h, &self.par_d[m])
        } else {
            0.0
        }
    }
}