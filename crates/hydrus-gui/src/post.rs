//! Post-processing pages.

use crate::util::*;
use crate::App;
use egui::{ComboBox, Ui};
use egui_plot::{Legend, Line, Plot, PlotPoints};

fn pl(pts: Vec<[f64; 2]>, name: &str, c: egui::Color32) -> Line<'static> {
    Line::new(PlotPoints::from(pts)).name(name).color(c)
}

impl App {
    fn no_results(&self, ui: &mut Ui) -> bool {
        if self.results.is_none() {
            ui.label("No results yet – run the calculation first (F5).");
            true
        } else {
            false
        }
    }

    pub fn page_obs(&mut self, ui: &mut Ui) {
        ui.heading("Observation points");
        if self.no_results(ui) {
            return;
        }
        let r = self.results.as_ref().unwrap();
        if r.obs.is_empty() {
            ui.label("No observation nodes were defined (add them in the soil profile editor).");
            return;
        }
        let nobs = r.obs[0].points.len();
        if self.post.obs_sel.len() != nobs {
            self.post.obs_sel = vec![true; nobs];
        }
        let mut vars = vec!["Pressure head".to_string(), "Water content".into(), "Temperature".into(), "Flux".into()];
        for j in 0..r.n_solutes {
            vars.push(format!("Concentration {}", j + 1));
        }
        self.post.obs_var = self.post.obs_var.min(vars.len() - 1);
        ui.horizontal(|ui| {
            ui.label("Variable");
            ComboBox::from_id_salt("ov").selected_text(&vars[self.post.obs_var]).show_ui(ui, |ui| {
                for (i, v) in vars.iter().enumerate() {
                    ui.selectable_value(&mut self.post.obs_var, i, v);
                }
            });
        });
        ui.horizontal_wrapped(|ui| {
            for (k, p) in r.obs[0].points.iter().enumerate() {
                let d = self.prj.profile.nodes.get(p.node - 1).map(|n| n.depth).unwrap_or(0.0);
                ui.checkbox(&mut self.post.obs_sel[k], format!("node {} ({:.1} {})", p.node, d, self.prj.units.length_str()));
            }
        });
        let v = self.post.obs_var;
        let tu = self.prj.units.time_str();
        Plot::new("obs").height((ui.available_height() - 24.0).max(400.0)).legend(Legend::default()).x_axis_label(format!("time ({})", tu)).y_axis_label(vars[v].clone()).show(ui, |pu| {
            for k in 0..nobs {
                if !self.post.obs_sel[k] {
                    continue;
                }
                let node = r.obs[0].points[k].node;
                let pts: Vec<[f64; 2]> = r
                    .obs
                    .iter()
                    .map(|o| {
                        let p = &o.points[k];
                        let y = match v {
                            0 => p.h,
                            1 => p.theta,
                            2 => p.temp,
                            3 => p.flux,
                            j => p.conc.get(j - 4).copied().unwrap_or(0.0),
                        };
                        [o.t, y]
                    })
                    .collect();
                pu.line(pl(pts, &format!("node {}", node), PALETTE[k % PALETTE.len()]));
            }
        });
    }

    pub fn page_profile_info(&mut self, ui: &mut Ui) {
        ui.heading("Profile information");
        if self.no_results(ui) {
            return;
        }
        let r = self.results.as_ref().unwrap();
        let np = r.profiles.len();
        if self.post.prof_sel.len() != np {
            self.post.prof_sel = (0..np).map(|i| np <= 12 || i == 0 || i + 1 == np || i % (np / 8 + 1) == 0).collect();
        }
        let mut vars: Vec<String> = ["Pressure head", "Water content", "Hydraulic conductivity", "Water capacity", "Flux", "Root water uptake (sink)", "v / Ks(top)", "Temperature"].iter().map(|s| s.to_string()).collect();
        for j in 0..r.n_solutes {
            vars.push(format!("Concentration {}", j + 1));
        }
        for j in 0..r.n_solutes {
            vars.push(format!("Sorbed / immobile conc. {}", j + 1));
        }
        self.post.prof_var = self.post.prof_var.min(vars.len() - 1);
        ui.horizontal(|ui| {
            ui.label("Variable");
            ComboBox::from_id_salt("pvv").selected_text(&vars[self.post.prof_var]).show_ui(ui, |ui| {
                for (i, v) in vars.iter().enumerate() {
                    ui.selectable_value(&mut self.post.prof_var, i, v);
                }
            });
            if ui.button("All times").clicked() {
                self.post.prof_sel.iter_mut().for_each(|s| *s = true);
            }
            if ui.button("None").clicked() {
                self.post.prof_sel.iter_mut().for_each(|s| *s = false);
            }
        });
        ui.horizontal_wrapped(|ui| {
            for (k, p) in r.profiles.iter().enumerate() {
                ui.checkbox(&mut self.post.prof_sel[k], format!("t = {}", p.t));
            }
        });
        let v = self.post.prof_var;
        let ns = r.n_solutes;
        let lu = self.prj.units.length_str();
        Plot::new("profplot").height((ui.available_height() - 24.0).max(400.0)).legend(Legend::default()).x_axis_label(vars[v].clone()).y_axis_label(format!("depth ({})", lu)).show(ui, |pu| {
            let mut ci = 0;
            for (k, p) in r.profiles.iter().enumerate() {
                if !self.post.prof_sel[k] {
                    continue;
                }
                let pts: Vec<[f64; 2]> = p
                    .nodes
                    .iter()
                    .map(|n| {
                        let x = match v {
                            0 => n.h,
                            1 => n.theta,
                            2 => n.k,
                            3 => n.c,
                            4 => n.flux,
                            5 => n.sink,
                            6 => n.v_over_ks,
                            7 => n.temp,
                            j if j < 8 + ns => n.conc.get(j - 8).copied().unwrap_or(0.0),
                            j => n.sorb.get(j - 8 - ns).copied().unwrap_or(0.0),
                        };
                        [x, -n.depth]
                    })
                    .collect();
                pu.line(pl(pts, &format!("t = {}", p.t), PALETTE[ci % PALETTE.len()]));
                ci += 1;
            }
        });
    }

    pub fn page_fluxes(&mut self, ui: &mut Ui) {
        ui.heading("Boundary fluxes, heads and cumulative quantities");
        if self.no_results(ui) {
            return;
        }
        let r = self.results.as_ref().unwrap();
        let mut groups: Vec<String> = vec![
            "Water fluxes: rTop, vTop, vBot, vRoot".into(),
            "Cumulative water fluxes".into(),
            "Pressure heads: top, root zone, bottom".into(),
            "Surface: runoff, cumulative infiltration/evaporation".into(),
            "Soil water storage".into(),
            "Temperature at top and bottom".into(),
        ];
        for j in 0..r.n_solutes {
            groups.push(format!("Solute {}: boundary fluxes", j + 1));
            groups.push(format!("Solute {}: cumulative fluxes", j + 1));
            groups.push(format!("Solute {}: concentrations top / bottom", j + 1));
        }
        self.post.flux_var = self.post.flux_var.min(groups.len() - 1);
        ComboBox::from_id_salt("fg").width(420.0).selected_text(&groups[self.post.flux_var]).show_ui(ui, |ui| {
            for (i, g) in groups.iter().enumerate() {
                ui.selectable_value(&mut self.post.flux_var, i, g);
            }
        });
        let g = self.post.flux_var;
        let tu = self.prj.units.time_str();
        let series = |f: &dyn Fn(&hydrus_core::output::TLevel) -> f64| -> Vec<[f64; 2]> { r.tlevel.iter().map(|t| [t.t, f(t)]).collect() };
        Plot::new("fluxplot").height((ui.available_height() - 24.0).max(400.0)).legend(Legend::default()).x_axis_label(format!("time ({})", tu)).show(ui, |pu| {
            if g == 0 {
                pu.line(pl(series(&|t| t.r_top), "rTop (potential)", PALETTE[0]));
                pu.line(pl(series(&|t| t.v_top), "vTop (actual)", PALETTE[1]));
                pu.line(pl(series(&|t| t.v_bot), "vBot", PALETTE[2]));
                pu.line(pl(series(&|t| t.v_root), "vRoot (actual)", PALETTE[3]));
                pu.line(pl(series(&|t| t.r_root), "rRoot (potential)", PALETTE[4]));
            } else if g == 1 {
                pu.line(pl(series(&|t| t.cum_r_top), "Σ rTop", PALETTE[0]));
                pu.line(pl(series(&|t| t.cum_v_top), "Σ vTop", PALETTE[1]));
                pu.line(pl(series(&|t| t.cum_v_bot), "Σ vBot", PALETTE[2]));
                pu.line(pl(series(&|t| t.cum_v_root), "Σ vRoot", PALETTE[3]));
                pu.line(pl(series(&|t| t.cum_r_root), "Σ rRoot", PALETTE[4]));
            } else if g == 2 {
                pu.line(pl(series(&|t| t.h_top), "h top", PALETTE[0]));
                pu.line(pl(series(&|t| t.h_root), "h root zone", PALETTE[1]));
                pu.line(pl(series(&|t| t.h_bot), "h bottom", PALETTE[2]));
            } else if g == 3 {
                pu.line(pl(series(&|t| t.run_off), "runoff rate", PALETTE[0]));
                pu.line(pl(series(&|t| t.cum_run_off), "Σ runoff", PALETTE[1]));
                pu.line(pl(series(&|t| t.cum_infil), "Σ infiltration", PALETTE[2]));
                pu.line(pl(series(&|t| t.cum_evap), "Σ evaporation", PALETTE[3]));
            } else if g == 4 {
                pu.line(pl(series(&|t| t.volume), "water in profile", PALETTE[0]));
            } else if g == 5 {
                pu.line(pl(series(&|t| t.temp_top), "T top", PALETTE[0]));
                pu.line(pl(series(&|t| t.temp_bot), "T bottom", PALETTE[1]));
            } else {
                let j = (g - 6) / 3;
                let kind = (g - 6) % 3;
                let sr = |f: &dyn Fn(&hydrus_core::output::SoluteTLevel) -> f64| -> Vec<[f64; 2]> { r.tlevel.iter().filter_map(|t| t.solutes.get(j).map(|s| [t.t, f(s)])).collect() };
                match kind {
                    0 => {
                        pu.line(pl(sr(&|s| s.cv_top), "top flux", PALETTE[0]));
                        pu.line(pl(sr(&|s| s.cv_bot), "bottom flux", PALETTE[1]));
                        pu.line(pl(sr(&|s| s.cv_root), "root uptake", PALETTE[2]));
                    }
                    1 => {
                        pu.line(pl(sr(&|s| s.cum_top), "Σ top", PALETTE[0]));
                        pu.line(pl(sr(&|s| s.cum_bot), "Σ bottom", PALETTE[1]));
                        pu.line(pl(sr(&|s| s.cum_ch0), "Σ zero-order production", PALETTE[2]));
                        pu.line(pl(sr(&|s| s.cum_ch1), "Σ first-order decay", PALETTE[3]));
                        pu.line(pl(sr(&|s| s.cum_root), "Σ root uptake", PALETTE[4]));
                        pu.line(pl(sr(&|s| s.cum_neq), "Σ non-equilibrium exchange", PALETTE[5]));
                    }
                    _ => {
                        pu.line(pl(sr(&|s| s.c_top), "c top node", PALETTE[0]));
                        pu.line(pl(sr(&|s| s.c_bot), "c bottom node", PALETTE[1]));
                        pu.line(pl(sr(&|s| s.c_root), "c root zone", PALETTE[2]));
                    }
                }
            }
        });
    }

    pub fn page_soil_props(&mut self, ui: &mut Ui) {
        ui.heading("Soil hydraulic properties");
        let nm = self.prj.water.materials.len();
        if nm == 0 {
            ui.label("No materials defined.");
            return;
        }
        ui.horizontal(|ui| {
            for i in 0..nm {
                let name = &self.prj.water.materials[i].name;
                if ui.selectable_label(self.post.prop_mat == i, format!("{} {}", i + 1, name)).clicked() {
                    self.post.prop_mat = i;
                }
            }
        });
        let i = self.post.prop_mat.min(nm - 1);
        self.soil_plots_pub(ui, i);
        ui.label("Curves are plotted from the current material parameters and hydraulic model (no hysteresis).");
    }

    pub fn page_runtime(&mut self, ui: &mut Ui) {
        ui.heading("Run-time information");
        if self.no_results(ui) {
            return;
        }
        let r = self.results.as_ref().unwrap();
        ui.label(format!(
            "{} time levels stored; wall-clock time {:.2} s; cumulative iterations: {}.",
            r.tlevel.len(),
            self.run_seconds,
            r.tlevel.last().map(|t| t.it_cum).unwrap_or(0)
        ));
        for m in &r.messages {
            ui.colored_label(egui::Color32::LIGHT_RED, m);
        }
        let tu = self.prj.units.time_str();
        let s = |f: &dyn Fn(&hydrus_core::output::TLevel) -> f64| -> Vec<[f64; 2]> { r.tlevel.iter().map(|t| [t.t, f(t)]).collect() };
        ui.columns(2, |c| {
            Plot::new("dtplot").height(260.0).x_axis_label(format!("time ({})", tu)).y_axis_label("time step Δt").show(&mut c[0], |pu| pu.line(pl(s(&|t| t.dt), "dt", PALETTE[0])));
            Plot::new("itplot").height(260.0).x_axis_label(format!("time ({})", tu)).y_axis_label("iterations per step").show(&mut c[1], |pu| pu.line(pl(s(&|t| t.iter_w as f64), "iterations", PALETTE[1])));
        });
        if r.n_solutes > 0 {
            ui.columns(2, |c| {
                Plot::new("peplot").height(260.0).x_axis_label(format!("time ({})", tu)).y_axis_label("max. Peclet number").show(&mut c[0], |pu| pu.line(pl(s(&|t| t.peclet), "Pe", PALETTE[2])));
                Plot::new("coplot").height(260.0).x_axis_label(format!("time ({})", tu)).y_axis_label("max. Courant number").show(&mut c[1], |pu| pu.line(pl(s(&|t| t.courant), "Cr", PALETTE[3])));
            });
        }
    }

    pub fn page_balance(&mut self, ui: &mut Ui) {
        ui.heading("Mass balance");
        if self.no_results(ui) {
            return;
        }
        let r = self.results.as_ref().unwrap();
        egui::Grid::new("bal").striped(true).show(ui, |ui| {
            for h in ["time", "water volume", "WatBalT", "WatBalR [%]"] {
                ui.strong(h);
            }
            for j in 0..r.n_solutes {
                ui.strong(format!("mass {}", j + 1));
                ui.strong(format!("CncBalT {}", j + 1));
                ui.strong(format!("CncBalR {} [%]", j + 1));
            }
            ui.end_row();
            for b in &r.balance {
                ui.label(format!("{}", b.t));
                ui.label(format!("{:.5}", b.total.volume));
                ui.label(format!("{:.3e}", b.wat_bal_t));
                ui.label(format!("{:.4}", b.wat_bal_r));
                for j in 0..r.n_solutes {
                    ui.label(format!("{:.5}", b.total.c_vol.get(j).copied().unwrap_or(0.0)));
                    ui.label(format!("{:.3e}", b.sol_bal_t.get(j).copied().unwrap_or(0.0)));
                    ui.label(format!("{:.4}", b.sol_bal_r.get(j).copied().unwrap_or(0.0)));
                }
                ui.end_row();
            }
        });
        if let Some(b) = r.balance.last() {
            section(ui, &format!("Sub-regions at t = {}", b.t));
            egui::Grid::new("sub").striped(true).show(ui, |ui| {
                for h in ["region", "length", "water volume", "inflow", "mean h"] {
                    ui.strong(h);
                }
                for j in 0..r.n_solutes {
                    ui.strong(format!("solute {} mass", j + 1));
                    ui.strong(format!("solute {} mean c", j + 1));
                }
                ui.end_row();
                for (i, s) in b.sub.iter().enumerate() {
                    ui.label(format!("{}", i + 1));
                    ui.label(format!("{:.4}", s.area));
                    ui.label(format!("{:.5}", s.volume));
                    ui.label(format!("{:.4e}", s.change));
                    ui.label(format!("{:.3}", s.h_mean));
                    for j in 0..r.n_solutes {
                        ui.label(format!("{:.5}", s.c_vol.get(j).copied().unwrap_or(0.0)));
                        ui.label(format!("{:.5}", s.c_mean.get(j).copied().unwrap_or(0.0)));
                    }
                    ui.end_row();
                }
            });
        }
    }
}
