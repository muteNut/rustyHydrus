//! Pre-processing pages.

use crate::util::*;
use crate::{App, Page};
use egui::{ComboBox, RichText, Ui};
use egui_extras::{Column, TableBuilder};
use egui_plot::{Line, Plot, PlotPoints};
use hydrus_core::*;

#[derive(PartialEq, Clone, Copy)]
enum TopKind {
    ConstHead,
    ConstFlux,
    AtmLayer,
    AtmRunoff,
    VarHead,
    VarFlux,
    VarHeadFlux,
}

#[derive(PartialEq, Clone, Copy)]
enum BotKind {
    ConstHead,
    ConstFlux,
    VarHead,
    VarFlux,
    Free,
    Gwl,
    Seep,
    Drains,
}

fn top_kind(b: &WaterBc) -> TopKind {
    if b.atmospheric {
        if b.surface_layer {
            TopKind::AtmLayer
        } else {
            TopKind::AtmRunoff
        }
    } else if b.top_time_variable {
        match b.kod_top {
            1 => TopKind::VarHead,
            -1 => TopKind::VarFlux,
            _ => TopKind::VarHeadFlux,
        }
    } else if b.kod_top == 1 {
        TopKind::ConstHead
    } else {
        TopKind::ConstFlux
    }
}

fn set_top_kind(b: &mut WaterBc, k: TopKind) {
    let (atm, tv, sl, kod) = match k {
        TopKind::ConstHead => (false, false, false, 1),
        TopKind::ConstFlux => (false, false, false, -1),
        TopKind::AtmLayer => (true, true, true, -1),
        TopKind::AtmRunoff => (true, true, false, -1),
        TopKind::VarHead => (false, true, false, 1),
        TopKind::VarFlux => (false, true, false, -1),
        TopKind::VarHeadFlux => (false, true, false, 0),
    };
    b.atmospheric = atm;
    b.top_time_variable = tv;
    b.surface_layer = sl;
    b.kod_top = kod;
}

fn bot_kind(b: &WaterBc) -> BotKind {
    if b.drains.is_some() {
        BotKind::Drains
    } else if b.seepage_face {
        BotKind::Seep
    } else if b.free_drainage {
        BotKind::Free
    } else if b.gwl_flux {
        BotKind::Gwl
    } else if b.bot_time_variable {
        if b.kod_bot >= 0 {
            BotKind::VarHead
        } else {
            BotKind::VarFlux
        }
    } else if b.kod_bot == 1 {
        BotKind::ConstHead
    } else {
        BotKind::ConstFlux
    }
}

fn set_bot_kind(b: &mut WaterBc, k: BotKind) {
    b.free_drainage = false;
    b.gwl_flux = false;
    b.seepage_face = false;
    b.drains = None;
    b.bot_time_variable = false;
    match k {
        BotKind::ConstHead => b.kod_bot = 1,
        BotKind::ConstFlux => b.kod_bot = -1,
        BotKind::VarHead => {
            b.bot_time_variable = true;
            b.kod_bot = 1;
        }
        BotKind::VarFlux => {
            b.bot_time_variable = true;
            b.kod_bot = -1;
        }
        BotKind::Free => {
            b.free_drainage = true;
            b.kod_bot = -1;
        }
        BotKind::Gwl => {
            b.gwl_flux = true;
            b.kod_bot = -1;
        }
        BotKind::Seep => {
            b.seepage_face = true;
            b.kod_bot = -1;
        }
        BotKind::Drains => {
            b.drains = Some(DrainSettings::default());
            b.kod_bot = -1;
        }
    }
}

impl App {
    /// Keep all per-material / per-species vectors consistent with the project.
    pub fn sync(&mut self) {
        let nm = self.prj.water.materials.len();
        let ns = self.prj.solute.species.len();
        let p = &mut self.prj;
        while p.heat.materials.len() < nm {
            p.heat.materials.push(p.heat.materials.last().cloned().unwrap_or_default());
        }
        while p.solute.materials.len() < nm {
            p.solute.materials.push(p.solute.materials.last().cloned().unwrap_or_default());
        }
        for sp in p.solute.species.iter_mut() {
            while sp.per_material.len() < nm {
                sp.per_material.push(sp.per_material.last().cloned().unwrap_or_default());
            }
        }
        if let RootStress::Feddes { p_optm, .. } = &mut p.root.stress {
            while p_optm.len() < nm {
                p_optm.push(*p_optm.last().unwrap_or(&-25.0));
            }
        }
        while p.root.c_root_max.len() < ns {
            p.root.c_root_max.push(0.0);
        }
        if let Some(s) = p.root.solute_stress.as_mut() {
            while s.a_osm.len() < ns {
                s.a_osm.push(0.0);
            }
        }
        for a in p.atmosphere.records.iter_mut() {
            while a.c_top.len() < ns {
                a.c_top.push(0.0);
            }
            while a.c_bot.len() < ns {
                a.c_bot.push(0.0);
            }
        }
        for n in p.profile.nodes.iter_mut() {
            while n.conc.len() < ns {
                n.conc.push(0.0);
            }
            while n.sorb.len() < ns {
                n.sorb.push(0.0);
            }
            n.mat = n.mat.clamp(1, nm);
        }
        self.sel_mat = self.sel_mat.min(nm.saturating_sub(1));
    }

    pub fn page_main(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Main processes & units");
        egui::Grid::new("main").num_columns(2).show(ui, |ui| {
            ui.label("Project title");
            if ui.text_edit_singleline(&mut self.prj.title).changed() {
                *ch = true;
            }
            ui.end_row();
            ui.label("Length unit");
            ComboBox::from_id_salt("lu").selected_text(self.prj.units.length_str()).show_ui(ui, |ui| {
                for (u, n) in [(LengthUnit::Mm, "mm"), (LengthUnit::Cm, "cm"), (LengthUnit::M, "m")] {
                    if ui.selectable_value(&mut self.prj.units.length, u, n).changed() {
                        *ch = true;
                    }
                }
            });
            ui.end_row();
            ui.label("Time unit");
            ComboBox::from_id_salt("tu").selected_text(self.prj.units.time_str()).show_ui(ui, |ui| {
                for (u, n) in [(TimeUnit::Sec, "seconds"), (TimeUnit::Min, "minutes"), (TimeUnit::Hours, "hours"), (TimeUnit::Days, "days"), (TimeUnit::Years, "years")] {
                    if ui.selectable_value(&mut self.prj.units.time, u, n).changed() {
                        *ch = true;
                    }
                }
            });
            ui.end_row();
            ui.label("Mass unit");
            if ui.text_edit_singleline(&mut self.prj.units.mass).changed() {
                *ch = true;
            }
            ui.end_row();
        });
        ui.label(RichText::new("Changing units does not convert existing parameter values.").weak());
        section(ui, "Processes");
        check(ui, ch, &mut self.prj.processes.water_flow, "Water flow (uncheck for transport under steady flow from the initial condition)");
        check(ui, ch, &mut self.prj.processes.vapor, "Vapor flow (coupled non-isothermal liquid & water vapor flow)");
        check(ui, ch, &mut self.prj.processes.root_water_uptake, "Root water uptake");
        check(ui, ch, &mut self.prj.processes.root_growth, "Root growth");
        check(ui, ch, &mut self.prj.processes.heat, "Heat transport");
        check(ui, ch, &mut self.prj.processes.solute, "Solute transport");
        if self.prj.processes.solute {
            ui.horizontal(|ui| {
                ui.label("Number of solutes:");
                let mut n = self.prj.solute.species.len();
                if ui.add(egui::DragValue::new(&mut n).range(1..=11)).changed() {
                    while self.prj.solute.species.len() < n {
                        let k = self.prj.solute.species.len() + 1;
                        let mut sp = Species::default();
                        sp.name = format!("Solute {}", k);
                        self.prj.solute.species.push(sp);
                    }
                    self.prj.solute.species.truncate(n);
                    *ch = true;
                }
            });
        }
        check(ui, ch, &mut self.prj.processes.short_output, "Short output: write time-level information only at print times");
        let mut has_meteo = self.prj.atmosphere.meteo.is_some();
        if ui.checkbox(&mut has_meteo, "Meteorological ET (Penman-Monteith / Hargreaves)").changed() {
            self.prj.atmosphere.meteo = if has_meteo {
                Some(MeteoSettings::default())
            } else {
                None
            };
            *ch = true;
        }
        self.sync();
        section(ui, "Modules in Development");
        ui.label("inverse (parameter estimation) optimization, and major ion chemistry (UnsatChem).");
    }

    pub fn page_geometry(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Geometry & materials");
        egui::Grid::new("geo").num_columns(2).show(ui, |ui| {
            let mut depth = self.prj.profile.depth();
            ui.label("Depth of the soil profile");
            if ui.add(egui::DragValue::new(&mut depth).speed(1.0).suffix(format!(" {}", self.prj.units.length_str()))).changed() && depth > 0.0 {
                let n = self.prj.profile.nodes.len();
                self.prj.regrid_uniform(depth, n);
                *ch = true;
            }
            ui.end_row();
            ui.label("Number of nodes");
            let mut n = self.prj.profile.nodes.len();
            if ui.add(egui::DragValue::new(&mut n).range(3..=1001)).changed() {
                let d = self.prj.profile.depth();
                self.prj.regrid_uniform(d, n);
                *ch = true;
            }
            ui.end_row();
            field(ui, ch, "cos α (angle between flow direction and vertical)", &mut self.prj.cos_alpha);
        });
        section(ui, "Soil materials");
        ui.horizontal(|ui| {
            if ui.button("➕ Add material").clicked() {
                let mut m = self.prj.water.materials.last().cloned().unwrap_or_default();
                m.name = format!("Material {}", self.prj.water.materials.len() + 1);
                self.prj.water.materials.push(m);
                self.sync();
                *ch = true;
            }
            if ui.add_enabled(self.prj.water.materials.len() > 1, egui::Button::new("➖ Remove last")).clicked() {
                self.prj.water.materials.pop();
                let nm = self.prj.water.materials.len();
                self.prj.solute.materials.truncate(nm);
                self.prj.heat.materials.truncate(nm);
                for s in self.prj.solute.species.iter_mut() {
                    s.per_material.truncate(nm);
                }
                self.sync();
                *ch = true;
            }
        });
        for (i, m) in self.prj.water.materials.iter_mut().enumerate() {
            let col = soil_color_for_name(&m.name, i + 1);
            ui.horizontal(|ui| {
                ui.colored_label(col, "■");
                ui.label(format!("{}:", i + 1));
                if ui.text_edit_singleline(&mut m.name).changed() {
                    *ch = true;
                }
            });
        }
        ui.label("Assign materials to depths on the “Soil profile editor” page; edit parameters under “Water flow: soil parameters”.");
    }

    pub fn page_nodes(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Nodal table");
        ui.label("Nodes are listed from the top (surface) downward. Root uptake distribution b(x) is normalised automatically.");
        let ns = if self.prj.processes.solute { self.prj.solute.species.len() } else { 0 };
        let heat = self.prj.processes.heat;
        let nm = self.prj.water.materials.len();
        let mut cols = vec!["#", "depth", "h / θ", "mat", "layer", "b(x)", "Ah", "AK", "ATh"];
        if heat {
            cols.push("T");
        }
        let n = self.prj.profile.nodes.len();
        let mut tb = TableBuilder::new(ui).striped(true).max_scroll_height(600.0);
        for _ in 0..(cols.len() + 2 * ns) {
            tb = tb.column(Column::auto().at_least(58.0));
        }
        let nodes = &mut self.prj.profile.nodes;
        tb.header(22.0, |mut h| {
            for c in &cols {
                h.col(|ui| {
                    ui.strong(*c);
                });
            }
            for j in 0..ns {
                h.col(|ui| {
                    ui.strong(format!("c{}", j + 1));
                });
            }
            for j in 0..ns {
                h.col(|ui| {
                    ui.strong(format!("s{}", j + 1));
                });
            }
        })
        .body(|body| {
            body.rows(20.0, n, |mut row| {
                let i = row.index();
                let nd = &mut nodes[i];
                row.col(|ui| {
                    ui.label(format!("{}", i + 1));
                });
                row.col(|ui| dv(ui, ch, &mut nd.depth));
                row.col(|ui| dv(ui, ch, &mut nd.h));
                row.col(|ui| {
                    if ui.add(egui::DragValue::new(&mut nd.mat).range(1..=nm)).changed() {
                        *ch = true;
                    }
                });
                row.col(|ui| di(ui, ch, &mut nd.layer));
                row.col(|ui| dv(ui, ch, &mut nd.beta));
                row.col(|ui| dv(ui, ch, &mut nd.ah));
                row.col(|ui| dv(ui, ch, &mut nd.ak));
                row.col(|ui| dv(ui, ch, &mut nd.ath));
                if heat {
                    row.col(|ui| dv(ui, ch, &mut nd.temp));
                }
                for j in 0..ns {
                    if let Some(c) = nd.conc.get_mut(j) {
                        row.col(|ui| dv(ui, ch, c));
                    } else {
                        row.col(|ui| { ui.label("0.0"); });
                    }
                }
                for j in 0..ns {
                    if let Some(s) = nd.sorb.get_mut(j) {
                        row.col(|ui| dv(ui, ch, s));
                    } else {
                        row.col(|ui| { ui.label("0.0"); });
                    }
                }
            });
        });
    }

    pub fn page_profile(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Soil profile & materials");
        self.sync();
        let len = self.prj.units.length_str();

        ui.horizontal_top(|ui| {
            let total_w = ui.available_width();
            let left_w = (total_w * 0.48).max(440.0);

            // ---- Left Column: Geometry, Materials & Form Controls ----
            ui.allocate_ui_with_layout(
                egui::vec2(left_w, ui.available_height()),
                egui::Layout::top_down(egui::Align::Min),
                |ui| {
                    ui.set_max_width(left_w);

                    // 1. Column Geometry
                    section(ui, "Profile geometry");
                    egui::Grid::new("geo").num_columns(2).show(ui, |ui| {
                        let mut depth = self.prj.profile.depth();
                        ui.label("Profile depth");
                        if ui.add(egui::DragValue::new(&mut depth).speed(1.0).suffix(format!(" {}", len))).changed() && depth > 0.0 {
                            let n = self.prj.profile.nodes.len();
                            self.prj.regrid_uniform(depth, n);
                            self.edit.to = depth;
                            *ch = true;
                        }
                        ui.end_row();

                        ui.label("Number of nodes");
                        let mut n = self.prj.profile.nodes.len();
                        if ui.add(egui::DragValue::new(&mut n).range(3..=1001)).changed() {
                            let d = self.prj.profile.depth();
                            self.prj.regrid_uniform(d, n);
                            self.edit.n_nodes = n;
                            *ch = true;
                        }
                        ui.end_row();

                        field(ui, ch, "Flow angle cos α (1 = vertical)", &mut self.prj.cos_alpha);
                    });

                    // 2. Soil Materials Manager
                    section(ui, "Soil materials");
                    ui.horizontal(|ui| {
                        if ui.button("➕ Add material").clicked() {
                            let mut m = self.prj.water.materials.last().cloned().unwrap_or_default();
                            m.name = format!("Material {}", self.prj.water.materials.len() + 1);
                            self.prj.water.materials.push(m);
                            self.sync();
                            *ch = true;
                        }
                        if ui.add_enabled(self.prj.water.materials.len() > 1, egui::Button::new("➖ Remove last")).clicked() {
                            self.prj.water.materials.pop();
                            let nm = self.prj.water.materials.len();
                            self.prj.solute.materials.truncate(nm);
                            self.prj.heat.materials.truncate(nm);
                            for s in self.prj.solute.species.iter_mut() {
                                s.per_material.truncate(nm);
                            }
                            self.sync();
                            *ch = true;
                        }
                    });

                    for (i, m) in self.prj.water.materials.iter_mut().enumerate() {
                        let col = soil_color_for_name(&m.name, i + 1);
                        ui.horizontal(|ui| {
                            ui.colored_label(col, "■");
                            ui.label(format!("{}:", i + 1));
                            if ui.text_edit_singleline(&mut m.name).changed() {
                                *ch = true;
                            }
                        });
                    }

                    // 3. Edit Depth Range & Assign Materials
                    section(ui, "Assign to depth range");
                    let depth = self.prj.profile.depth();
                    if self.edit.to <= self.edit.from {
                        self.edit.to = depth;
                    }
                    egui::Grid::new("edit").num_columns(2).show(ui, |ui| {
                        ui.label("From depth");
                        ui.add(egui::DragValue::new(&mut self.edit.from).speed(0.5).suffix(format!(" {}", len)));
                        ui.end_row();

                        ui.label("To depth");
                        ui.add(egui::DragValue::new(&mut self.edit.to).speed(0.5).suffix(format!(" {}", len)));
                        ui.end_row();

                        ui.label("Material");
                        let cur_m = self.edit.mat.clamp(1, self.prj.water.materials.len());
                        let cur_name = &self.prj.water.materials[cur_m - 1].name;
                        ComboBox::from_id_salt("prof_mat_dd")
                            .selected_text(format!("{}: {}", cur_m, cur_name))
                            .show_ui(ui, |ui| {
                                for (m_idx, mat) in self.prj.water.materials.iter().enumerate() {
                                    let label = format!("{}: {}", m_idx + 1, mat.name);
                                    let col = soil_mat_color(&self.prj, m_idx + 1);
                                    ui.horizontal(|ui| {
                                        ui.colored_label(col, "■");
                                        if ui.selectable_value(&mut self.edit.mat, m_idx + 1, label).changed() {
                                            *ch = true;
                                        }
                                    });
                                }
                            });
                        ui.end_row();
                    });

                    let (from, to) = (self.edit.from, self.edit.to);
                    let sel = |d: f64| d >= from - 1e-9 && d <= to + 1e-9;
                    ui.horizontal(|ui| {
                        if ui.button("Assign material").clicked() {
                            for n in self.prj.profile.nodes.iter_mut().filter(|n| sel(n.depth)) {
                                n.mat = self.edit.mat;
                                n.layer = self.edit.mat;
                            }
                            *ch = true;
                        }
                    });

                    ui.separator();
                    ui.horizontal(|ui| {
                        ui.label(if self.prj.water.init_in_water_content { "θ top" } else { "h top" });
                        ui.add(egui::DragValue::new(&mut self.edit.h_top).speed(1.0));
                        ui.label(if self.prj.water.init_in_water_content { "θ bottom" } else { "h bottom" });
                        ui.add(egui::DragValue::new(&mut self.edit.h_bot).speed(1.0));
                        if ui.button("Set initial condition (linear)").clicked() {
                            let (a, b) = (from, to);
                            for n in self.prj.profile.nodes.iter_mut().filter(|n| n.depth >= a - 1e-9 && n.depth <= b + 1e-9) {
                                let f = if b > a { (n.depth - a) / (b - a) } else { 0.0 };
                                n.h = self.edit.h_top + (self.edit.h_bot - self.edit.h_top) * f;
                            }
                            *ch = true;
                        }
                    });

                    ui.horizontal(|ui| {
                        ui.label("Root uptake weight b");
                        ui.add(egui::DragValue::new(&mut self.edit.beta).speed(0.01));
                        if ui.button("Set").clicked() {
                            for n in self.prj.profile.nodes.iter_mut().filter(|n| sel(n.depth)) {
                                n.beta = self.edit.beta;
                            }
                            *ch = true;
                        }
                    });

                    if self.prj.processes.heat {
                        ui.horizontal(|ui| {
                            ui.label("Temperature");
                            ui.add(egui::DragValue::new(&mut self.edit.temp).speed(0.5));
                            if ui.button("Set").clicked() {
                                for n in self.prj.profile.nodes.iter_mut().filter(|n| sel(n.depth)) {
                                    n.temp = self.edit.temp;
                                }
                                *ch = true;
                            }
                        });
                    }

                    if self.prj.processes.solute {
                        ui.horizontal(|ui| {
                            ui.label("Concentration (solute");
                            ui.add(egui::DragValue::new(&mut self.edit.view).range(1..=self.prj.solute.species.len().max(1)));
                            ui.label(")");
                            ui.add(egui::DragValue::new(&mut self.edit.conc).speed(0.1));
                            if ui.button("Set").clicked() {
                                let j = self.edit.view.max(1) - 1;
                                for n in self.prj.profile.nodes.iter_mut().filter(|n| sel(n.depth)) {
                                    if j < n.conc.len() {
                                        n.conc[j] = self.edit.conc;
                                    }
                                }
                                *ch = true;
                            }
                        });
                    }

                    check(ui, ch, &mut self.prj.water.init_in_water_content, "Initial condition given as water content θ");

                    // 4. Observation Nodes
                    section(ui, "Observation nodes");
                    ui.horizontal(|ui| {
                        ui.label("Depth");
                        ui.add(egui::DragValue::new(&mut self.edit.obs_depth).speed(0.5).suffix(format!(" {}", len)));
                        if ui.button("Add nearest node").clicked() {
                            let mut best = 1;
                            let mut bd = f64::MAX;
                            for (i, n) in self.prj.profile.nodes.iter().enumerate() {
                                let d = (n.depth - self.edit.obs_depth).abs();
                                if d < bd {
                                    bd = d;
                                    best = i + 1;
                                }
                            }
                            if !self.prj.profile.observation_nodes.contains(&best) && self.prj.profile.observation_nodes.len() < 100 {
                                self.prj.profile.observation_nodes.push(best);
                                self.prj.profile.observation_nodes.sort();
                                *ch = true;
                            }
                        }
                        if ui.button("Clear").clicked() {
                            self.prj.profile.observation_nodes.clear();
                            *ch = true;
                        }
                    });

                    let obs: Vec<String> = self.prj.profile.observation_nodes
                        .iter()
                        .filter_map(|&k| {
                            if k >= 1 && k <= self.prj.profile.nodes.len() {
                                Some(format!("{} ({:.1} {})", k, self.prj.profile.nodes[k - 1].depth, len))
                            } else {
                                None
                            }
                        })
                        .collect();
                    ui.label(format!("Observation nodes: {}", if obs.is_empty() { "none".into() } else { obs.join(", ") }));

                    // 5. Nodal Density & Regridding
                    section(ui, "Nodal density");
                    ui.label("Element length is proportional to density. Fixed points: depth, top density, bot density.");
                    let mut remove = None;
                    for (i, f) in self.edit.fixed.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.add(egui::DragValue::new(&mut f.0).speed(0.5).prefix("z "));
                            ui.add(egui::DragValue::new(&mut f.1).speed(0.05).prefix("top "));
                            ui.add(egui::DragValue::new(&mut f.2).speed(0.05).prefix("bot "));
                            if ui.small_button("✖").clicked() {
                                remove = Some(i);
                            }
                        });
                    }
                    if let Some(i) = remove {
                        self.edit.fixed.remove(i);
                    }
                    ui.horizontal(|ui| {
                        if ui.button("Add fixed point").clicked() {
                            self.edit.fixed.push((depth, 1.0, 1.0));
                        }
                        ui.add(egui::DragValue::new(&mut self.edit.n_nodes).range(3..=1001).prefix("nodes "));
                        if ui.button("Generate mesh").clicked() {
                            let z = generate_nodes(depth, self.edit.n_nodes, &self.edit.fixed);
                            let old = self.prj.profile.nodes.clone();
                            let mut nodes = vec![];
                            for &zz in &z {
                                let mut best = old[0].clone();
                                let mut bd = f64::MAX;
                                for o in &old {
                                    let d = (o.depth - zz).abs();
                                    if d < bd {
                                        bd = d;
                                        best = o.clone();
                                    }
                                }
                                best.depth = zz;
                                nodes.push(best);
                            }
                            self.prj.profile.nodes = nodes;
                            self.prj.profile.observation_nodes.clear();
                            *ch = true;
                        }
                    });
                },
            );

            ui.add_space(16.0);
            ui.separator();
            ui.add_space(16.0);

            // ---- Right Column: Plots & Strata ----
            ui.vertical(|ui| {
                let plot_w = ui.available_width().clamp(340.0, 720.0);
                let nodes = &self.prj.profile.nodes;

                // 1. Initial State Plot
                let series: Vec<(&str, Vec<[f64; 2]>, egui::Color32)> = vec![
                    ("initial h / θ", nodes.iter().map(|n| [n.h, -n.depth]).collect(), PALETTE[0]),
                ];
                Plot::new("prof_init")
                    .width(plot_w)
                    .height(210.0)
                    .x_axis_label("initial pressure head / water content")
                    .y_axis_label(format!("depth ({})", len))
                    .legend(egui_plot::Legend::default())
                    .show(ui, |pu| {
                        for (name, pts, c) in series {
                            pu.line(Line::new(PlotPoints::from(pts)).name(name).color(c));
                        }
                    });

                ui.add_space(8.0);
                ui.label(egui::RichText::new("Soil Strata Horizon").strong());

                // 2. Earth Stratum Horizon Plot
                struct Stratum {
                    mat_idx: usize,
                    z_top: f64, // positive depth downward
                    z_bot: f64,
                }
                let mut strata: Vec<Stratum> = vec![];
                if !nodes.is_empty() {
                    let mut cur_mat = nodes[0].mat;
                    let mut z_top = nodes[0].depth;
                    for i in 1..nodes.len() {
                        let d = nodes[i].depth;
                        if nodes[i].mat != cur_mat {
                            strata.push(Stratum { mat_idx: cur_mat.saturating_sub(1), z_top, z_bot: d });
                            cur_mat = nodes[i].mat;
                            z_top = d;
                        }
                    }
                    if let Some(last) = nodes.last() {
                        strata.push(Stratum { mat_idx: cur_mat.saturating_sub(1), z_top, z_bot: last.depth });
                    }
                }

                let mut hovered_coord: Option<egui_plot::PlotPoint> = None;

                let plot_resp = Plot::new("prof_strata")
                    .width(plot_w)
                    .height(180.0)
                    .show_x(false)
                    .show_axes([false, true])
                    .y_axis_label(format!("depth ({})", len))
                    .include_x(0.0)
                    .include_x(1.0)
                    .allow_zoom(false)
                    .allow_drag(false)
                    .legend(egui_plot::Legend::default())
                    .show(ui, |pu| {
                        // Capture pointer coordinates in plot-space
                        hovered_coord = pu.pointer_coordinate();

                        for s in &strata {
                            let base_col = soil_mat_color(&self.prj, s.mat_idx + 1);
							let col = egui::Color32::from_rgba_unmultiplied(base_col.r(), base_col.g(), base_col.b(), 140);
                            let name = self.prj.water.materials.get(s.mat_idx)
                                .map(|m| m.name.clone())
                                .unwrap_or_else(|| format!("Material {}", s.mat_idx + 1));

                            let is_selected = (self.edit.from - s.z_top).abs() < 1e-4 
                                && (self.edit.to - s.z_bot).abs() < 1e-4;
                            let stroke_width = if is_selected { 2.5_f32 } else { 1.0_f32 };
                            let stroke_col = if is_selected {
                                egui::Color32::WHITE
                            } else {
                                col.linear_multiply(0.7_f32)
                            };

                            let poly = egui_plot::Polygon::new(PlotPoints::new(vec![
                                [0.0, -s.z_top],
                                [1.0, -s.z_top],
                                [1.0, -s.z_bot],
                                [0.0, -s.z_bot],
                            ]))
                            .fill_color(col)
                            .stroke(egui::Stroke::new(stroke_width, stroke_col))
                            .name(name);

                            pu.polygon(poly);
                        }
                    });

                // Detect click from the plot response and update edit fields
                if plot_resp.response.clicked() {
                    if let Some(pos) = hovered_coord {
                        let depth_clicked = -pos.y;
                        for s in &strata {
                            let min_z = s.z_top.min(s.z_bot);
                            let max_z = s.z_top.max(s.z_bot);
                            if depth_clicked >= min_z && depth_clicked <= max_z {
                                self.edit.from = min_z;
                                self.edit.to = max_z;
                                self.edit.mat = s.mat_idx + 1;
                                break;
                            }
                        }
                    }
                }

                ui.add_space(8.0);

                // 3. Root Uptake Distribution Plot
                Plot::new("prof_beta")
                    .width(plot_w)
                    .height(170.0)
                    .x_axis_label("root distribution b(x)")
                    .y_axis_label(format!("depth ({})", len))
                    .show(ui, |pu| {
                        let pts: Vec<[f64; 2]> = nodes.iter().map(|n| [n.beta, -n.depth]).collect();
                        pu.line(Line::new(PlotPoints::from(pts)).color(PALETTE[2]));
                    });
            });
        });
    }

    pub fn page_time(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Time & printing");
        let tu = self.prj.units.time_str().to_string();
        egui::Grid::new("time").num_columns(2).show(ui, |ui| {
            let t = &mut self.prj.time;
            field_u(ui, ch, "Initial time step Δt", &mut t.dt, &tu);
            field_u(ui, ch, "Minimum time step", &mut t.dt_min, &tu);
            field_u(ui, ch, "Maximum time step", &mut t.dt_max, &tu);
            field(ui, ch, "Multiplication factor (few iterations)", &mut t.d_mul);
            field(ui, ch, "Multiplication factor (many iterations)", &mut t.d_mul2);
            ui.label("Lower optimal iteration range");
            di(ui, ch, &mut t.it_min);
            ui.end_row();
            ui.label("Upper optimal iteration range");
            di(ui, ch, &mut t.it_max);
            ui.end_row();
            field_u(ui, ch, "Initial time", &mut t.t_init, &tu);
            field_u(ui, ch, "Final time", &mut t.t_max, &tu);
        });
        section(ui, "Printing");
        if self.print_text_for != self.prj.time.print_times {
            self.print_text = self.prj.time.print_times.iter().map(|v| format!("{}", v)).collect::<Vec<_>>().join(", ");
            self.print_text_for = self.prj.time.print_times.clone();
        }
        ui.label("Print times (profiles and mass balance are stored at these times), comma-separated:");
        if ui.add(egui::TextEdit::multiline(&mut self.print_text).desired_width(500.0).desired_rows(2)).changed() {
            let v: Vec<f64> = self.print_text.split(|c: char| c == ',' || c.is_whitespace() || c == ';').filter_map(|s| s.trim().parse().ok()).collect();
            self.prj.time.print_times = v.clone();
            self.print_text_for = v;
            *ch = true;
        }
        ui.horizontal(|ui| {
            if ui.button("Generate equidistant print times").clicked() {
                let n = 10;
                let (a, b) = (self.prj.time.t_init, self.prj.time.t_max);
                self.prj.time.print_times = (1..=n).map(|i| a + (b - a) * i as f64 / n as f64).collect();
                *ch = true;
            }
        });
        egui::Grid::new("time2").num_columns(2).show(ui, |ui| {
            ui.label("Store time-level info every n-th step");
            di(ui, ch, &mut self.prj.time.print_step);
            ui.end_row();
        });
        check(ui, ch, &mut self.prj.time.print_at_interval, "Additionally print at regular intervals");
        if self.prj.time.print_at_interval {
            egui::Grid::new("time3").show(ui, |ui| field_u(ui, ch, "Interval", &mut self.prj.time.print_interval, &tu));
        }
    }

    pub fn page_iteration(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Water flow: iteration criteria & lookup table");
        let lu = self.prj.units.length_str().to_string();
        egui::Grid::new("it").num_columns(2).show(ui, |ui| {
            let w = &mut self.prj.water;
            ui.label("Maximum number of iterations");
            di(ui, ch, &mut w.max_iter);
            ui.end_row();
            field(ui, ch, "Water content tolerance", &mut w.tol_th);
            field_u(ui, ch, "Pressure head tolerance", &mut w.tol_h, &lu);
            field_u(ui, ch, "Lower limit of the tension interval (table)", &mut w.h_tab1, &lu);
            field_u(ui, ch, "Upper limit of the tension interval (table)", &mut w.h_tab_n, &lu);
        });
        ui.label(RichText::new("Hydraulic properties are pre-tabulated (100 log-spaced points) and interpolated within this interval, as in HYDRUS-1D.").weak());
    }

    pub fn page_hydmodel(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Water flow: hydraulic model");
        ui.horizontal(|ui| {
            ui.label("Soil hydraulic model");
            ComboBox::from_id_salt("model").selected_text(self.prj.water.model.name()).show_ui(ui, |ui| {
                for m in [
                    SoilModel::VanGenuchten,
                    SoilModel::ModifiedVG,
                    SoilModel::BrooksCorey,
                    SoilModel::VGAirEntry,
                    SoilModel::Kosugi,
                    SoilModel::Durner,
                    SoilModel::DualPorosityW,
                    SoilModel::DualPorosityH,
                    SoilModel::DualPermeability,
                    SoilModel::Tabular,
                ] {
                    if ui.selectable_value(&mut self.prj.water.model, m, m.name()).changed() {
                        *ch = true;
                    }
                }
            });
        });
        if self.prj.water.model == SoilModel::Tabular {
            ui.horizontal(|ui| {
                let status = if self.prj.water.tabs.is_some() {
                    RichText::new("✔ External tables (Mater.in) loaded").color(egui::Color32::LIGHT_GREEN)
                } else {
                    RichText::new("⚠ No external tables loaded (Mater.in required)").color(egui::Color32::LIGHT_RED)
                };
                ui.label(status);
            });
        }
        ui.horizontal(|ui| {
            ui.label("Hysteresis");
            let cur = match self.prj.water.hysteresis {
				Hysteresis::None => "No hysteresis",
				Hysteresis::Retention => "Hysteresis in the retention curve",
				Hysteresis::RetentionAndConductivity => "Hysteresis in retention and conductivity",
				Hysteresis::Lenhard => "Lenhard et al. air-entrapment hysteresis",
			};
			ComboBox::from_id_salt("hyst").selected_text(cur).show_ui(ui, |ui| {
				for (h, n) in [
					(Hysteresis::None, "No hysteresis"),
					(Hysteresis::Retention, "Hysteresis in the retention curve"),
					(Hysteresis::RetentionAndConductivity, "Hysteresis in retention and conductivity"),
					(Hysteresis::Lenhard, "Lenhard et al. air-entrapment hysteresis"),
				] {
					if ui.selectable_value(&mut self.prj.water.hysteresis, h, n).changed() {
						*ch = true;
					}
				}
			});
        });
        if self.prj.water.hysteresis != Hysteresis::None {
            ui.horizontal(|ui| {
                ui.label("Initial condition lies on the");
                if ui.selectable_value(&mut self.prj.water.init_kappa, -1, "main drying curve").changed() {
                    *ch = true;
                }
                if ui.selectable_value(&mut self.prj.water.init_kappa, 1, "main wetting curve").changed() {
                    *ch = true;
                }
            });
        }
    }

    fn soil_plots(&self, ui: &mut Ui, mi: usize, id: &str) {
        let p = &self.prj;
        let xc = p.units.x_conv();
        let m = &p.water.materials[mi];
        let par = material::par_of(m, p.water.model, xc);
        let model = p.water.model;
        let lu = p.units.length_str();

        if model == SoilModel::Tabular {
            if let Some(ref tabs) = p.water.tabs {
                if let Some(tb) = tabs.get(mi) {
                    let mut th = vec![];
                    let mut k = vec![];
                    let mut c = vec![];
                    for idx in 0..tb.h.len() {
                        let x = (-tb.h[idx].min(-1e-6)).log10();
                        th.push([x, tb.the[idx]]);
                        k.push([x, tb.con[idx].max(1e-30).log10()]);
                        c.push([x, tb.cap[idx].max(1e-30).log10()]);
                    }
                    ui.columns(3, |cols| {
                        Plot::new(format!("{}th", id)).height(240.0).x_axis_label(format!("log10 |h| ({})", lu)).y_axis_label("θ").show(&mut cols[0], |pu| pu.line(Line::new(PlotPoints::from(th)).color(PALETTE[0])));
                        Plot::new(format!("{}k", id)).height(240.0).x_axis_label(format!("log10 |h| ({})", lu)).y_axis_label("log10 K").show(&mut cols[1], |pu| pu.line(Line::new(PlotPoints::from(k)).color(PALETTE[1])));
                        Plot::new(format!("{}c", id)).height(240.0).x_axis_label(format!("log10 |h| ({})", lu)).y_axis_label("log10 C").show(&mut cols[2], |pu| pu.line(Line::new(PlotPoints::from(c)).color(PALETTE[2])));
                    });
                    return;
                }
            }
        }

        let hs = log_space(1e-2 * xc, 1e5 * xc, 200);
        let (mut th, mut k, mut c) = (vec![], vec![], vec![]);
        for h in hs {
            let x = h.log10();
            th.push([x, material::fq(model, -h, &par)]);
            k.push([x, material::fk(model, -h, &par).max(1e-30).log10()]);
            c.push([x, material::fc(model, -h, &par).max(1e-30).log10()]);
        }
        ui.columns(3, |cols| {
            Plot::new(format!("{}th", id)).height(240.0).x_axis_label(format!("log10 |h| ({})", lu)).y_axis_label("θ").show(&mut cols[0], |pu| pu.line(Line::new(PlotPoints::from(th)).color(PALETTE[0])));
            Plot::new(format!("{}k", id)).height(240.0).x_axis_label(format!("log10 |h| ({})", lu)).y_axis_label("log10 K").show(&mut cols[1], |pu| pu.line(Line::new(PlotPoints::from(k)).color(PALETTE[1])));
            Plot::new(format!("{}c", id)).height(240.0).x_axis_label(format!("log10 |h| ({})", lu)).y_axis_label("log10 C").show(&mut cols[2], |pu| pu.line(Line::new(PlotPoints::from(c)).color(PALETTE[2])));
        });
    }

    pub fn soil_plots_pub(&self, ui: &mut Ui, mi: usize) {
        self.soil_plots(ui, mi, "sp");
    }

    pub fn page_soil(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Water flow: soil hydraulic parameters");
        self.sync();
        let model = self.prj.water.model;
        let hyst = self.prj.water.hysteresis != Hysteresis::None;
        let (lf, tf) = (cm_to_unit(&self.prj), day_to_unit(&self.prj));
        let lu = self.prj.units.length_str();
        ui.label(format!("Units: α in 1/{}, Ks in {}/{}", lu, lu, self.prj.units.time_str()));
        let nm = self.prj.water.materials.len();
        ui.horizontal_wrapped(|ui| {
            for i in 0..nm {
                let sel = self.sel_mat == i;
                let col = soil_mat_color(&self.prj, i + 1);
                if ui.selectable_label(sel, RichText::new(format!("{} {}", i + 1, self.prj.water.materials[i].name)).color(col)).clicked() {
                    self.sel_mat = i;
                }
            }
        });
        let i = self.sel_mat;
        ui.separator();
        ui.horizontal(|ui| {
            ui.label("Load from catalog (Carsel & Parrish 1988):");
            ComboBox::from_id_salt("cat").selected_text("choose…").show_ui(ui, |ui| {
                for (k, c) in SOIL_CATALOG.iter().enumerate() {
                    if ui.selectable_label(false, c.0).clicked() {
                        apply_catalog(&mut self.prj.water.materials[i], k, lf, tf);
                        *ch = true;
                    }
                }
            });
        });
        egui::Grid::new("soilp").num_columns(2).show(ui, |ui| {
            let m = &mut self.prj.water.materials[i];

            field(ui, ch, "θr – residual water content", &mut m.qr);
            field(ui, ch, "θs – saturated water content", &mut m.qs);
            match model {
                SoilModel::BrooksCorey => {
                    field(ui, ch, "α – inverse of air-entry value", &mut m.alpha);
                    field(ui, ch, "n – pore-size distribution index (λ)", &mut m.n);
                }
                SoilModel::Kosugi => {
                    field(ui, ch, "hm – median pressure head (as α)", &mut m.alpha);
                    field(ui, ch, "σ – log-standard deviation (as n)", &mut m.n);
                }
                _ => {
                    field(ui, ch, "α", &mut m.alpha);
                    field(ui, ch, "n", &mut m.n);
                }
            }
            field(ui, ch, "Ks – saturated hydraulic conductivity", &mut m.ks);
            field(ui, ch, "l – pore connectivity", &mut m.l);
            for (k, lab) in model.extra_labels().iter().enumerate() {
                if k < m.extra.len() {
                    field(ui, ch, lab, &mut m.extra[k]);
                }
            }
            if hyst {
                field(ui, ch, "θm (drying)", &mut m.qm);
                field(ui, ch, "θs – wetting branch", &mut m.qs_w);
                field(ui, ch, "α – wetting branch", &mut m.alpha_w);
                field(ui, ch, "Ks – wetting branch", &mut m.ks_w);
            }
        });
        ui.separator();
        self.soil_plots(ui, i, "soilp");
    }

    pub fn page_waterbc(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Water flow: boundary conditions");
        let lu = self.prj.units.length_str().to_string();
        let fu = format!("{}/{}", lu, self.prj.units.time_str());
        let bc = &mut self.prj.water.bc;
        section(ui, "Upper boundary");
        let mut tk = top_kind(bc);
        let names = [
            (TopKind::ConstHead, "Constant pressure head (from the initial condition at the top node)"),
            (TopKind::ConstFlux, "Constant flux"),
            (TopKind::AtmLayer, "Atmospheric BC with surface layer"),
            (TopKind::AtmRunoff, "Atmospheric BC with surface runoff"),
            (TopKind::VarHead, "Variable pressure head"),
            (TopKind::VarFlux, "Variable flux (Evap − Prec)"),
            (TopKind::VarHeadFlux, "Variable pressure head / flux (code in the evaporation column)"),
        ];
        let cur = names.iter().find(|(k, _)| *k == tk).map(|x| x.1).unwrap_or("");
        ComboBox::from_id_salt("topbc").width(480.0).selected_text(cur).show_ui(ui, |ui| {
            for (k, n) in names {
                if ui.selectable_value(&mut tk, k, n).changed() {
                    set_top_kind(bc, k);
                    *ch = true;
                }
            }
        });
        if tk == TopKind::ConstFlux {
            egui::Grid::new("q").show(ui, |ui| field_u(ui, ch, "Top flux (negative = infiltration)", &mut bc.r_top, &fu));
        }
        if matches!(tk, TopKind::AtmLayer | TopKind::AtmRunoff) {
            ui.label("Precipitation, potential evaporation/transpiration and critical head are read from the atmospheric data table.");
        }
        section(ui, "Lower boundary");
        let mut bk = bot_kind(bc);
        let names = [
            (BotKind::ConstHead, "Constant pressure head (initial condition at the bottom node)"),
            (BotKind::ConstFlux, "Constant flux"),
            (BotKind::VarHead, "Variable pressure head"),
            (BotKind::VarFlux, "Variable flux"),
            (BotKind::Free, "Free drainage"),
            (BotKind::Gwl, "Deep drainage: q = A·exp(B·|h − GWL0|)"),
            (BotKind::Seep, "Seepage face"),
            (BotKind::Drains, "Horizontal drains"),
        ];
        let cur = names.iter().find(|(k, _)| *k == bk).map(|x| x.1).unwrap_or("");
        ComboBox::from_id_salt("botbc").width(480.0).selected_text(cur).show_ui(ui, |ui| {
            for (k, n) in names {
                if ui.selectable_value(&mut bk, k, n).changed() {
                    set_bot_kind(bc, k);
                    *ch = true;
                }
            }
        });
        match bk {
            BotKind::ConstFlux => {
                egui::Grid::new("qb").show(ui, |ui| field_u(ui, ch, "Bottom flux (positive = upward)", &mut bc.r_bot, &fu));
            }
            BotKind::Gwl => {
                egui::Grid::new("gwl").show(ui, |ui| {
                    field_u(ui, ch, "Reference groundwater level GWL0", &mut bc.gwl0l, &lu);
                    field_u(ui, ch, "Aqh", &mut bc.aqh, &fu);
                    field(ui, ch, "Bqh", &mut bc.bqh);
                });
            }
            BotKind::Seep => {
                egui::Grid::new("seep").show(ui, |ui| {
                    field_u(ui, ch, "Pressure head at the seepage face", &mut bc.h_seep, &lu);
                });
                ui.label(RichText::new("Dynamic switching between zero flux (unsaturated) and Dirichlet head (saturated).").weak());
            }
            BotKind::Drains => {
                if let Some(d) = bc.drains.as_mut() {
                    ui.horizontal(|ui| {
                        ui.label("Drain geometry / configuration");
                        ComboBox::from_id_salt("dpos").selected_text(["", "1: on impervious layer", "2: above impervious layer", "3: at layer interface", "4: in bottom layer", "5: in top layer"][d.position.clamp(1, 5) as usize]).show_ui(ui, |ui| {
                            for k in 1..=5 {
                                if ui.selectable_value(&mut d.position, k, ["", "1: homogeneous, on impervious layer", "2: homogeneous, above impervious layer", "3: heterogeneous, at layer interface", "4: heterogeneous, in bottom layer", "5: heterogeneous, in top layer"][k as usize]).changed() {
                                    *ch = true;
                                }
                            }
                        });
                    });
                    egui::Grid::new("drain").show(ui, |ui| {
                        field_u(ui, ch, "Elevation of the drain bottom (negative)", &mut d.z_bot, &lu);
                        field_u(ui, ch, "Drain spacing", &mut d.spacing, &lu);
                        field(ui, ch, "Entrance resistance", &mut d.entrance_resistance);
                        if d.position >= 2 {
                            field_u(ui, ch, "Elevation of the impervious layer (negative)", &mut d.base_gw, &lu);
                            field_u(ui, ch, "Wet perimeter", &mut d.wet_perimeter, &lu);
                        }
                        field_u(ui, ch, "Kh above drain", &mut d.kh_top, &fu);
                        if matches!(d.position, 3 | 4 | 5) {
                            field_u(ui, ch, "Kh below drain", &mut d.kh_bot, &fu);
                        }
                        if matches!(d.position, 4 | 5) {
                            field_u(ui, ch, "Kv top", &mut d.kv_top, &fu);
                            field_u(ui, ch, "Depth of layer interface", &mut d.z_interface, &lu);
                        }
                        if d.position == 4 {
                            field_u(ui, ch, "Kv bottom", &mut d.kv_bot, &fu);
                        }
                        if d.position == 5 {
                            field(ui, ch, "Geometry factor", &mut d.geo_factor);
                        }
                    });
                }
            }
            _ => {}
        }
        if bc.top_time_variable || bc.bot_time_variable || bc.atmospheric {
            egui::Grid::new("hc").show(ui, |ui| field_u(ui, ch, "hCritS – max. surface layer height before runoff", &mut self.prj.atmosphere.h_crit_s, &lu));
        }
        if self.prj.processes.root_water_uptake && !(bc.top_time_variable || bc.atmospheric) {
            egui::Grid::new("rr").show(ui, |ui| field_u(ui, ch, "Constant potential transpiration", &mut bc.r_root, &fu));
        }
    }

    pub fn import_atmosphere_text(&mut self, path: &std::path::Path) -> Result<usize, String> {
        let s = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
        let ns = if self.prj.processes.solute { self.prj.solute.species.len() } else { 0 };
        let heat = self.prj.processes.heat;
        let root = self.prj.atmosphere.has_root_depth;
        let mut recs = vec![];
        for l in s.lines() {
            let v: Vec<f64> = l.split(|c: char| c.is_whitespace() || c == ',' || c == ';').filter(|t| !t.is_empty()).filter_map(|t| t.replace(['d', 'D'], "e").parse().ok()).collect();
            if v.len() < 8 {
                continue;
            }
            let mut r = AtmRecord { t: v[0], prec: v[1], evap: v[2], transp: v[3], h_crit_a: v[4], r_bot: v[5], h_bot: v[6], h_top: v[7], ..Default::default() };
            let mut k = 8;
            if heat && v.len() >= k + 3 {
                r.t_top = v[k];
                r.t_bot = v[k + 1];
                r.ampl = v[k + 2];
                k += 3;
            }
            for _ in 0..ns {
                r.c_top.push(v.get(k).copied().unwrap_or(0.0));
                r.c_bot.push(v.get(k + 1).copied().unwrap_or(0.0));
                k += 2;
            }
            if root {
                r.x_root = v.get(k).copied().unwrap_or(0.0);
            }
            recs.push(r);
        }
        if recs.is_empty() {
            return Err("No numeric records with at least 8 columns found.".into());
        }
        let n = recs.len();
        self.prj.atmosphere.records = recs;
        self.dirty = true;
        Ok(n)
    }

    pub fn page_atmosphere(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Atmospheric data");
        self.sync();
        let tu = self.prj.units.time_str().to_string();
        ui.label(format!("Rates are per {}. Columns: time, precipitation, potential evaporation, potential transpiration, hCritA (minimum surface pressure head), bottom flux, bottom head, top head{}.", tu, if self.prj.processes.heat { ", T top, T bottom, amplitude" } else { "" }));
        ui.horizontal(|ui| {
            check(ui, ch, &mut self.prj.atmosphere.daily_variation, "Daily sinusoidal variation of evaporation/transpiration");
            check(ui, ch, &mut self.prj.atmosphere.sinusoidal_precip, "Sinusoidal precipitation within each interval");
        });
        ui.horizontal(|ui| {
            check(ui, ch, &mut self.prj.atmosphere.snow, "Snow accumulation and melt");
            if self.prj.atmosphere.snow {
                ui.label("Snow melt factor:");
                dv(ui, ch, &mut self.prj.atmosphere.snow_mf);
            }
            check(ui, ch, &mut self.prj.atmosphere.interception, "Canopy interception");
            if self.prj.atmosphere.interception {
                ui.label("Interception a:");
                dv(ui, ch, &mut self.prj.atmosphere.interception_a);
            }
        });
        ui.horizontal(|ui| {
            if ui.button("➕ Add record").clicked() {
                let mut r = self.prj.atmosphere.records.last().cloned().unwrap_or(AtmRecord { h_crit_a: 10000.0, ..Default::default() });
                r.t += 1.0;
                self.prj.atmosphere.records.push(r);
                self.sync();
                *ch = true;
            }
            if ui.button("➖ Remove last").clicked() {
                self.prj.atmosphere.records.pop();
                *ch = true;
            }
            if ui.button("Import from text file…").clicked() {
                self.dialog = Some(crate::files::FileDialog::new(crate::files::Purpose::ImportAtmosphere, crate::files::Mode::OpenFile, "Import atmospheric records (whitespace/CSV columns)", self.last_dir.clone(), "", &[]));
            }
            ui.label(format!("{} records", self.prj.atmosphere.records.len()));
        });
        let ns = if self.prj.processes.solute { self.prj.solute.species.len() } else { 0 };
        let heat = self.prj.processes.heat;
        let root = self.prj.atmosphere.has_root_depth;
        let mut heads = vec!["t".to_string(), "Prec".into(), "Evap".into(), "Transp".into(), "hCritA".into(), "rBot".into(), "hBot".into(), "hTop".into()];
        if heat {
            heads.extend(["Ttop".to_string(), "Tbot".into(), "Ampl".into()]);
        }
        for j in 0..ns {
            heads.push(format!("cTop{}", j + 1));
            heads.push(format!("cBot{}", j + 1));
        }
        if root {
            heads.push("xRoot".into());
        }
        let n = self.prj.atmosphere.records.len();
        let mut tb = TableBuilder::new(ui).striped(true).max_scroll_height(320.0);
        for _ in 0..heads.len() {
            tb = tb.column(Column::exact(82.0));
        }
        let recs = &mut self.prj.atmosphere.records;
        tb.header(22.0, |mut h| {
            for t in &heads {
                h.col(|ui| {
                    ui.strong(t);
                });
            }
        })
        .body(|body| {
            body.rows(20.0, n, |mut row| {
                let r = &mut recs[row.index()];
                row.col(|ui| dv(ui, ch, &mut r.t));
                row.col(|ui| dv(ui, ch, &mut r.prec));
                row.col(|ui| dv(ui, ch, &mut r.evap));
                row.col(|ui| dv(ui, ch, &mut r.transp));
                row.col(|ui| dv(ui, ch, &mut r.h_crit_a));
                row.col(|ui| dv(ui, ch, &mut r.r_bot));
                row.col(|ui| dv(ui, ch, &mut r.h_bot));
                row.col(|ui| dv(ui, ch, &mut r.h_top));
                if heat {
                    row.col(|ui| dv(ui, ch, &mut r.t_top));
                    row.col(|ui| dv(ui, ch, &mut r.t_bot));
                    row.col(|ui| dv(ui, ch, &mut r.ampl));
                }
                for j in 0..ns {
                    row.col(|ui| dv(ui, ch, &mut r.c_top[j]));
                    row.col(|ui| dv(ui, ch, &mut r.c_bot[j]));
                }
                if root {
                    row.col(|ui| dv(ui, ch, &mut r.x_root));
                }
            });
        });
        ui.add_space(8.0);
        let recs = &self.prj.atmosphere.records;
        if !recs.is_empty() {
            Plot::new("atmplot").height(240.0).legend(egui_plot::Legend::default()).x_axis_label(format!("time ({})", tu)).show(ui, |pu| {
                let step = |f: &dyn Fn(&AtmRecord) -> f64| -> Vec<[f64; 2]> {
                    let mut pts = vec![];
                    let mut t0 = self.prj.time.t_init;
                    for r in recs {
                        pts.push([t0, f(r)]);
                        pts.push([r.t, f(r)]);
                        t0 = r.t;
                    }
                    pts
                };
                pu.line(Line::new(PlotPoints::from(step(&|r| r.prec))).name("precipitation").color(PALETTE[0]));
                pu.line(Line::new(PlotPoints::from(step(&|r| r.evap))).name("pot. evaporation").color(PALETTE[1]));
                pu.line(Line::new(PlotPoints::from(step(&|r| r.transp))).name("pot. transpiration").color(PALETTE[2]));
            });
        }
    }

    pub fn page_meteo(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Meteorological parameters");

        let Some(mp) = self.prj.atmosphere.meteo.as_mut() else {
            ui.label("Meteorological ET is not enabled. Turn it on in 'Main processes & units'.");
            return;
        };

        let tu = self.prj.units.time_str().to_string();

        section(ui, "Radiation & Cloudiness Input Options");
        ui.columns(2, |cols| {
            cols[0].group(|ui| {
                ui.strong("Radiation Input");
                ui.radio_value(&mut mp.i_radiation, 0, "Potential Radiation (Ra computed)");
                ui.radio_value(&mut mp.i_radiation, 1, "Solar Radiation (measured Rs)");
                ui.radio_value(&mut mp.i_radiation, 2, "Net Radiation (measured Rn)");
            });

            cols[1].group(|ui| {
                ui.strong("Cloudiness / Sunshine Input");
                ui.radio_value(&mut mp.i_sun_sh, 0, "Sunshine Hours");
                ui.radio_value(&mut mp.i_sun_sh, 1, "Cloudiness Index");
                ui.radio_value(&mut mp.i_sun_sh, 2, "Transmission Coeff. / Cloud Cover");
                ui.radio_value(&mut mp.i_sun_sh, 3, "Derived from Solar Radiation");
            });
        });

        ui.horizontal(|ui| {
            ui.label("Humidity specification:");
            ui.radio_value(&mut mp.i_rel_hum, 0, "Relative Humidity (%)");
            ui.radio_value(&mut mp.i_rel_hum, 1, "Vapor Pressure (kPa)");
            check(ui, ch, &mut mp.hargreaves, "Use Hargreaves equation instead of Penman-Monteith");
        });

        section(ui, "Geographical & Radiation Parameters");
        egui::Grid::new("meteo_geo").num_columns(4).striped(true).show(ui, |ui| {
            ui.label("Latitude (° N+, S-)");
            dv(ui, ch, &mut mp.latitude);
            ui.label("Altitude (m)");
            dv(ui, ch, &mut mp.altitude);
            ui.end_row();

            ui.label("Ångström a (shortwave)");
            dv(ui, ch, &mut mp.short_wave_a);
            ui.label("Ångström b (shortwave)");
            dv(ui, ch, &mut mp.short_wave_b);
            ui.end_row();

            ui.label("Longwave cloud factor a1");
            dv(ui, ch, &mut mp.long_wave_a);
            ui.label("Longwave cloud factor b1");
            dv(ui, ch, &mut mp.long_wave_b);
            ui.end_row();

            ui.label("Longwave emissivity al");
            dv(ui, ch, &mut mp.long_wave_a1);
            ui.label("Longwave emissivity bl");
            dv(ui, ch, &mut mp.long_wave_b1);
            ui.end_row();

            if mp.i_sun_sh == 3 {
                ui.label("Cloudiness from Solar ac");
                dv(ui, ch, &mut mp.cloud_fact_ac);
                ui.label("Cloudiness from Solar bc");
                dv(ui, ch, &mut mp.cloud_fact_bc);
                ui.end_row();
            }

            ui.label("Wind measurement height (cm)");
            dv(ui, ch, &mut mp.wind_height);
            ui.label("Temp. measurement height (cm)");
            dv(ui, ch, &mut mp.temp_height);
            ui.end_row();
        });

        section(ui, "Crop & Surface Parameters");
        ui.horizontal(|ui| {
            ui.label("Crop mode:");
            ui.radio_value(&mut mp.i_crop, 0, "Bare Soil");
            ui.radio_value(&mut mp.i_crop, 1, "Active Crop");
        });

        if mp.i_crop != 0 {
            ui.horizontal(|ui| {
                ui.label("Leaf Area Index (LAI):");
                ui.radio_value(&mut mp.i_lai, 1, "From Height (Clipped Grass)");
                ui.radio_value(&mut mp.i_lai, 2, "From Height (Alfalfa)");
                ui.radio_value(&mut mp.i_lai, 3, "Direct Input / Surface Fraction");
            });

            egui::Grid::new("meteo_crop").num_columns(4).striped(true).show(ui, |ui| {
                ui.label("Crop height (cm)");
                dv(ui, ch, &mut mp.crop_height);
                ui.label("Albedo");
                dv(ui, ch, &mut mp.albedo);
                ui.end_row();

                if mp.i_lai == 3 {
                    ui.label("Input LAI / Surface Fraction");
                    dv(ui, ch, &mut mp.lai);
                }
                ui.label("Radiation Extinction Coeff.");
                dv(ui, ch, &mut mp.r_extinct);
                ui.end_row();
            });
        }

        section(ui, "Meteorological Records");
        ui.horizontal(|ui| {
            if ui.button("➕ Add record").clicked() {
                let mut r = mp.records.last().cloned().unwrap_or(MeteoRecord {
                    t: 0.0,
                    rad: 15.0,
                    t_max: 20.0,
                    t_min: 10.0,
                    rh_mean: 50.0,
                    wind_kmd: 150.0,
                    sun_hours: 8.0,
                    ..Default::default()
                });
                r.t += 1.0;
                mp.records.push(r);
                *ch = true;
            }
            if ui.button("➖ Remove last").clicked() {
                mp.records.pop();
                *ch = true;
            }
            ui.label(format!("{} record(s)", mp.records.len()));
        });

        let heads = [
            format!("t ({})", tu),
            "Rad (MJ/m²)".into(),
            "T Max (°C)".into(),
            "T Min (°C)".into(),
            if mp.i_rel_hum == 0 { "RH Mean (%)".into() } else { "Vapor (kPa)".into() },
            "Wind (km/d)".into(),
            match mp.i_sun_sh {
                0 => "Sun Hours (h)",
                1 => "Cloudiness",
                2 => "Transm. Coeff",
                _ => "Solar Factor",
            }.into(),
        ];

        let num_records = mp.records.len();
        let mut tb = TableBuilder::new(ui).striped(true).max_scroll_height(260.0);
        for _ in 0..heads.len() {
            tb = tb.column(Column::exact(110.0));
        }

        let recs = &mut mp.records;
        tb.header(22.0, |mut h| {
            for title in &heads {
                h.col(|ui| {
                    ui.strong(title);
                });
            }
        })
        .body(|body| {
            body.rows(20.0, num_records, |mut row| {
                let r = &mut recs[row.index()];
                row.col(|ui| dv(ui, ch, &mut r.t));
                row.col(|ui| dv(ui, ch, &mut r.rad));
                row.col(|ui| dv(ui, ch, &mut r.t_max));
                row.col(|ui| dv(ui, ch, &mut r.t_min));
                row.col(|ui| dv(ui, ch, &mut r.rh_mean));
                row.col(|ui| dv(ui, ch, &mut r.wind_kmd));
                row.col(|ui| dv(ui, ch, &mut r.sun_hours));
            });
        });
    }

    pub fn page_root(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Root water uptake");
        self.sync();
        let lu = self.prj.units.length_str().to_string();
        let is_feddes = matches!(self.prj.root.stress, RootStress::Feddes { .. });
        ui.horizontal(|ui| {
            ui.label("Water stress response function");
            if ui.selectable_label(is_feddes, "Feddes et al. (1978)").clicked() && !is_feddes {
                self.prj.root.stress = RootUptake::default().stress;
                self.sync();
                *ch = true;
            }
            if ui.selectable_label(!is_feddes, "van Genuchten (1987) S-shaped").clicked() && is_feddes {
                self.prj.root.stress = RootStress::SShaped { p50: -800.0, exponent: 3.0 };
                *ch = true;
            }
        });
        egui::Grid::new("root").num_columns(2).show(ui, |ui| match &mut self.prj.root.stress {
            RootStress::Feddes { p0, p2h, p2l, p3, r2h, r2l, p_optm } => {
                field_u(ui, ch, "P0 (anaerobiosis)", p0, &lu);
                for (i, po) in p_optm.iter_mut().enumerate() {
                    field_u(ui, ch, &format!("POptm – material {}", i + 1), po, &lu);
                }
                field_u(ui, ch, "P2H (high demand)", p2h, &lu);
                field_u(ui, ch, "P2L (low demand)", p2l, &lu);
                field_u(ui, ch, "P3 (wilting)", p3, &lu);
                field(ui, ch, "r2H (potential transpiration, high)", r2h);
                field(ui, ch, "r2L (potential transpiration, low)", r2l);
            }
            RootStress::SShaped { p50, exponent } => {
                field_u(ui, ch, "P50", p50, &lu);
                field(ui, ch, "Exponent p", exponent);
            }
        });
        egui::Grid::new("root2").num_columns(2).show(ui, |ui| {
            field(ui, ch, "ωc – critical stress index for compensated uptake (1 = none)", &mut self.prj.root.omega_c);
        });
        section(ui, "Root distribution / growth");
        if !self.prj.processes.root_growth {
            ui.label("Root growth is off: the nodal distribution b(x) from the profile is used.");
        } else {
            let mode = match &self.prj.root.growth {
                RootGrowthMode::FromAtmosphere => 0,
                RootGrowthMode::Table(_) => 1,
                RootGrowthMode::Logistic { .. } => 2,
            };
            let mut m2 = mode;
            ui.horizontal(|ui| {
                ui.radio_value(&mut m2, 0, "Rooting depth from atmospheric data");
                ui.radio_value(&mut m2, 1, "Table (time, depth)");
                ui.radio_value(&mut m2, 2, "Logistic (Verhulst–Pearl)");
            });
            if m2 != mode {
                self.prj.root.growth = match m2 {
                    0 => RootGrowthMode::FromAtmosphere,
                    1 => RootGrowthMode::Table(vec![(self.prj.time.t_init, 1.0), (self.prj.time.t_max, 30.0)]),
                    _ => RootGrowthMode::Logistic { t_min: 0.0, t_med: 30.0, t_harv: self.prj.time.t_max, x_min: 1.0, x_med: 15.0, x_max: 30.0, period: 1e30, fit_from_half: false },
                };
                self.prj.atmosphere.has_root_depth = m2 == 0;
                *ch = true;
            }
            match &mut self.prj.root.growth {
                RootGrowthMode::Table(t) => {
                    for r in t.iter_mut() {
                        ui.horizontal(|ui| {
                            dv_unit(ui, ch, &mut r.0, "time");
                            dv_unit(ui, ch, &mut r.1, &lu);
                        });
                    }
                    ui.horizontal(|ui| {
                        if ui.button("Add row").clicked() {
                            let l = t.last().cloned().unwrap_or((0.0, 1.0));
                            t.push((l.0 + 1.0, l.1));
                            *ch = true;
                        }
                        if ui.button("Remove last").clicked() {
                            t.pop();
                            *ch = true;
                        }
                    });
                }
                RootGrowthMode::Logistic { t_min, t_med, t_harv, x_min, x_med, x_max, period, .. } => {
                    egui::Grid::new("log").show(ui, |ui| {
                        field(ui, ch, "Sowing time", t_min);
                        field(ui, ch, "Time of known root depth", t_med);
                        field(ui, ch, "Harvest time", t_harv);
                        field(ui, ch, "Initial root depth", x_min);
                        field(ui, ch, "Root depth at the known time", x_med);
                        field(ui, ch, "Maximum root depth", x_max);
                        field(ui, ch, "Period (1e30 = no cycling)", period);
                    });
                }
                _ => {}
            }
        }
        if self.prj.processes.solute {
            section(ui, "Solute stress");
            let mut on = self.prj.root.solute_stress.is_some();
            if ui.checkbox(&mut on, "Reduce uptake by osmotic / solute stress").changed() {
                self.prj.root.solute_stress = if on { Some(SoluteStress { additive: true, a_osm: vec![0.0; self.prj.solute.species.len()], c50: 3000.0, p3c: 3.0, s_shaped: true }) } else { None };
                *ch = true;
            }
            if let Some(s) = self.prj.root.solute_stress.as_mut() {
                check(ui, ch, &mut s.additive, "Additive (osmotic head added to pressure head)");
                if !s.additive {
                    check(ui, ch, &mut s.s_shaped, "S-shaped (else threshold-slope)");
                    egui::Grid::new("ss").show(ui, |ui| {
                        field(ui, ch, "c50 / threshold", &mut s.c50);
                        field(ui, ch, "exponent / slope", &mut s.p3c);
                    });
                }
                for (j, a) in s.a_osm.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        ui.label(format!("Osmotic coefficient, solute {}", j + 1));
                        dv(ui, ch, a);
                    });
                }
            }
            for (j, c) in self.prj.root.c_root_max.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(format!("Maximum concentration for passive uptake, solute {}", j + 1));
                    dv(ui, ch, c);
                });
            }
        }
    }

    pub fn page_sol_general(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Solute transport: general");
        let s = &mut self.prj.solute;
        egui::Grid::new("sg").num_columns(2).show(ui, |ui| {
            field(ui, ch, "Time weighting ε (0 explicit, 0.5 Crank–Nicolson, 1 implicit)", &mut s.epsi);
            ui.label("Max. iterations (non-linear)");
            di(ui, ch, &mut s.max_iter);
            ui.end_row();
            field(ui, ch, "Absolute concentration tolerance", &mut s.tol_abs);
            field(ui, ch, "Relative concentration tolerance", &mut s.tol_rel);
            field(ui, ch, "Stability criterion Pe_crit", &mut s.pe_cr);
        });
        check(ui, ch, &mut s.upstream_weighting, "Upstream weighting finite elements (else Galerkin)");
        if !s.upstream_weighting {
            check(ui, ch, &mut s.artificial_dispersion, "Add artificial dispersion to stabilise (Pe > Pe_crit)");
        }
        check(ui, ch, &mut s.tortuosity, "Tortuosity factor");
        if s.tortuosity {
            ui.horizontal(|ui| {
                ui.radio_value(&mut s.tort_model, 0, "Millington & Quirk");
                ui.radio_value(&mut s.tort_model, 1, "Moldrup & Olesen");
            });
        }
        check(ui, ch, &mut s.equil_init, "Non-equilibrium phase initially in equilibrium with the liquid phase");
        check(ui, ch, &mut s.mass_init, "Initial conditions given as total solute mass (per volume of soil)");
        check(ui, ch, &mut s.l_moist, "Water content dependent reaction rates (Walker / exponential)");
        check(ui, ch, &mut self.prj.root.l_act_rsu, "Active root solute uptake (Michaelis-Menten kinetics)");

        if self.prj.root.l_act_rsu {
            ui.indent("act_rsu", |ui| {
                ui.label("Active solute uptake parameters (OmegaAct, Km, cMin):");
                let ns = self.prj.solute.species.len();
                while self.prj.root.omega_act.len() < ns { self.prj.root.omega_act.push(1.0); }
                while self.prj.root.r_km.len() < ns { self.prj.root.r_km.push(0.5); }
                while self.prj.root.c_min.len() < ns { self.prj.root.c_min.push(0.0); }
                for j in 0..ns {
                    ui.horizontal(|ui| {
                        ui.label(format!("Solute {}:", j + 1));
                        ui.add(egui::DragValue::new(&mut self.prj.root.omega_act[j]).speed(0.05).prefix("ω_act: "));
                        ui.add(egui::DragValue::new(&mut self.prj.root.r_km[j]).speed(0.05).prefix("Km: "));
                        ui.add(egui::DragValue::new(&mut self.prj.root.c_min[j]).speed(0.01).prefix("cMin: "));
                    });
                }
            });
        }
        
        let l = self.prj.solute.species.len();
        section(ui, "Solutes");
        for i in 0..l {
            ui.horizontal(|ui| {
                ui.label(format!("{}:", i + 1));
                if ui.text_edit_singleline(&mut self.prj.solute.species[i].name).changed() {
                    *ch = true;
                }
            });
        }
        ui.label("Decay chain: solute k is produced from solute k−1 by the γ coefficients on the reaction page.");
    }

    pub fn page_sol_materials(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Solute transport: soil parameters");
        self.sync();
        let lu = self.prj.units.length_str().to_string();
        let nm = self.prj.water.materials.len();
        egui::Grid::new("sm").striped(true).show(ui, |ui| {
            for h in ["Material", "Bulk density ρ", &format!("Dispersivity Dl ({})", lu), "Fraction f (equilibrium sites)", "Immobile water θim"] {
                ui.strong(h);
            }
            ui.end_row();
            for i in 0..nm {
                let m = &mut self.prj.solute.materials[i];
                ui.colored_label(mat_color(i + 1), format!("{} {}", i + 1, self.prj.water.materials[i].name));
                dv(ui, ch, &mut m.bulk_density);
                dv(ui, ch, &mut m.disp_l);
                dv(ui, ch, &mut m.frac);
                dv(ui, ch, &mut m.th_immobile);
                ui.end_row();
            }
        });
        ui.label("Non-equilibrium transport is switched on automatically when f < 1 (two-site sorption) or θim > 0 (mobile–immobile water).");
        let ns = self.prj.solute.species.len();
        section(ui, "Molecular diffusion");
        egui::Grid::new("diff").striped(true).show(ui, |ui| {
            ui.strong("Solute");
            ui.strong("Dw (liquid)");
            ui.strong("Dg (gas)");
            ui.end_row();
            for j in 0..ns {
                let s = &mut self.prj.solute.species[j];
                ui.label(&s.name);
                dv(ui, ch, &mut s.diff_w);
                dv(ui, ch, &mut s.diff_g);
                ui.end_row();
            }
        });
    }

    pub fn page_sol_reactions(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Solute transport: reaction parameters");
        self.sync();
        let nm = self.prj.water.materials.len();
        let ns = self.prj.solute.species.len();
        ui.label("Adsorption s = Ks·c^β/(1+η·c^β); μ: first-order decay in liquid (w), solid (s), gas (g); γ: production from the parent solute; μ0: zero-order production; α: first-order mass-transfer rate for non-equilibrium.");
        for j in 0..ns {
            ui.add_space(6.0);
            ui.label(RichText::new(self.prj.solute.species[j].name.clone()).strong());
            egui::Grid::new(format!("re{}", j)).striped(true).show(ui, |ui| {
                for h in ["Mat", "Ks", "η (Nu)", "β", "Henry kH", "μw", "μs", "μg", "γw", "γs", "γg", "μ0w", "μ0s", "μ0g", "α"] {
                    ui.strong(h);
                }
                ui.end_row();
                for m in 0..nm {
                    let p = &mut self.prj.solute.species[j].per_material[m];
                    ui.colored_label(mat_color(m + 1), format!("{}", m + 1));
                    for v in [&mut p.ks, &mut p.nu, &mut p.beta, &mut p.henry, &mut p.mu_w, &mut p.mu_s, &mut p.mu_g, &mut p.gam_w, &mut p.gam_s, &mut p.gam_g, &mut p.mu0_w, &mut p.mu0_s, &mut p.mu0_g, &mut p.omega] {
                        dv(ui, ch, v);
                    }
                    ui.end_row();
                }
            });
        }
    }

    pub fn page_sol_bc(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Solute transport: boundary conditions");
        let tu = self.prj.units.time_str().to_string();
        let s = &mut self.prj.solute;
        ui.horizontal(|ui| {
            ui.label("Upper boundary");
            ui.radio_value(&mut s.k_top, 1, "Concentration BC");
            ui.radio_value(&mut s.k_top, -1, "Concentration flux BC");
        });
        ui.horizontal(|ui| {
            ui.label("Lower boundary");
            ui.radio_value(&mut s.k_bot, 1, "Concentration BC");
            ui.radio_value(&mut s.k_bot, -1, "Concentration flux BC");
            ui.radio_value(&mut s.k_bot, 0, "Zero gradient");
        });
        egui::Grid::new("sbc").show(ui, |ui| field_u(ui, ch, "Pulse duration (BC concentrations set to 0 afterwards, unless atmospheric)", &mut s.t_pulse, &tu));
        egui::Grid::new("sbc2").striped(true).show(ui, |ui| {
            ui.strong("Solute");
            ui.strong("c top");
            ui.strong("c bottom");
            ui.end_row();
            for sp in s.species.iter_mut() {
                ui.label(&sp.name);
                dv(ui, ch, &mut sp.c_top);
                dv(ui, ch, &mut sp.c_bot);
                ui.end_row();
            }
        });
        ui.label("With time-variable / atmospheric boundary conditions the top and bottom concentrations are taken from the atmospheric data table.");
    }

    pub fn page_heat(&mut self, ui: &mut Ui, ch: &mut bool) {
        ui.heading("Heat transport");
        self.sync();
        let h = &mut self.prj.heat;
        egui::Grid::new("ht").striped(true).show(ui, |ui| {
            for t in ["Material", "θn", "θo", "λL", "b1", "b2", "b3", "Cn", "Co", "Cw"] {
                ui.strong(t);
            }
            ui.end_row();
            for (i, m) in h.materials.iter_mut().enumerate() {
                ui.colored_label(mat_color(i + 1), format!("{}", i + 1));
                for v in [&mut m.qn, &mut m.qo, &mut m.disp, &mut m.b1, &mut m.b2, &mut m.b3, &mut m.cn, &mut m.co, &mut m.cw] {
                    dv(ui, ch, v);
                }
                ui.end_row();
            }
        });
        check(ui, ch, &mut h.campbell, "Chung & Horton (Campbell) thermal conductivity: θn = solid fraction, b1/b2/b3 = quartz, other minerals, clay fractions");
        ui.label("Thermal parameters must be given in the project's units (energy per length³ etc., as in HYDRUS-1D).");
        section(ui, "Boundary conditions");
        ui.horizontal(|ui| {
            ui.label("Upper");
            ui.radio_value(&mut h.k_top, 1, "Temperature");
            ui.radio_value(&mut h.k_top, -1, "Heat flux");
            ui.label("Lower");
            ui.radio_value(&mut h.k_bot, 1, "Temperature");
            ui.radio_value(&mut h.k_bot, -1, "Heat flux");
            ui.radio_value(&mut h.k_bot, 0, "Zero gradient");
        });
        egui::Grid::new("ht2").show(ui, |ui| {
            field(ui, ch, "Upper temperature", &mut h.t_top);
            field(ui, ch, "Lower temperature", &mut h.t_bot);
            field(ui, ch, "Amplitude of the sinusoidal top temperature", &mut h.amplitude);
            field(ui, ch, "Period", &mut h.period);
        });
        ui.label("If an atmospheric data table is used, its temperature columns override these values.");
    }

    pub fn page_summary(&mut self, ui: &mut Ui) {
        ui.heading("Check & summary");
        self.problems = self.prj.validate();
        if self.problems.is_empty() {
            ui.colored_label(egui::Color32::LIGHT_GREEN, "✔ No problems found – the project can be run.");
        } else {
            for p in &self.problems {
                ui.colored_label(egui::Color32::LIGHT_RED, format!("✖ {}", p));
            }
        }
        section(ui, "Overview");
        let p = &self.prj;
        ui.label(format!("Title: {}", p.title));
        ui.label(format!("Units: {}, {}, {}", p.units.length_str(), p.units.time_str(), p.units.mass));
        ui.label(format!("Profile: {} nodes over {:.3} {}; {} material(s); {} observation node(s)", p.profile.nodes.len(), p.profile.depth(), p.units.length_str(), p.water.materials.len(), p.profile.observation_nodes.len()));
        ui.label(format!("Hydraulic model: {}", p.water.model.name()));
        ui.label(format!("Simulation time: {} – {} {}", p.time.t_init, p.time.t_max, p.units.time_str()));
        let mut pr = vec![];
        if p.processes.water_flow {
            pr.push("water flow");
        }
        if p.processes.vapor {
            pr.push("water vapor flow");
        }
        if p.processes.root_water_uptake {
            pr.push("root water uptake");
        }
        if p.processes.solute {
            pr.push("solute transport");
        }
        if p.processes.heat {
            pr.push("heat transport");
        }
        if p.atmosphere.meteo.is_some() {
            pr.push("meteorological ET");
        }
        ui.label(format!("Processes: {}", pr.join(", ")));
        if ui.button("Go to hydraulic model page").clicked() {
            self.page = Page::HydModel;
        }
    }
}