#![cfg_attr(all(not(debug_assertions), windows), windows_subsystem = "windows")]

mod files;
mod post;
mod pre;
mod util;

use eframe::egui;
use files::{DialogResult, FileDialog, Mode, Purpose};
use hydrus_core::{Project, Results, Simulation};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::mpsc::{channel, Receiver};
use std::sync::Arc;
use std::time::Instant;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Page {
    Main,
    Geometry,
    Profile,
    NodeTable,
    Time,
    Iteration,
    HydModel,
    SoilParams,
    WaterBc,
    Atmosphere,
	Meteo,
    RootUptake,
    SoluteGeneral,
    SoluteMaterials,
    SoluteReactions,
    SoluteBc,
    Heat,
    Summary,
    Obs,
    ProfileInfo,
    Fluxes,
    SoilProps,
    RunTime,
    MassBalance,
}

struct RunHandle {
    rx: Receiver<Result<Results, String>>,
    progress: Arc<AtomicU64>,
    cancel: Arc<AtomicBool>,
    t_init: f64,
    t_max: f64,
    started: Instant,
}

#[derive(Default)]
pub struct EditState {
    pub from: f64,
    pub to: f64,
    pub mat: usize,
    pub h_top: f64,
    pub h_bot: f64,
    pub beta: f64,
    pub temp: f64,
    pub conc: f64,
    pub obs_depth: f64,
    pub n_nodes: usize,
    pub fixed: Vec<(f64, f64, f64)>,
    pub view: usize,
}

pub struct PostState {
    pub obs_var: usize,
    pub obs_sel: Vec<bool>,
    pub prof_var: usize,
    pub prof_sel: Vec<bool>,
    pub flux_var: usize,
    pub solute: usize,
    pub run_var: usize,
    pub prop_mat: usize,
    pub prop_kind: usize,
}
impl Default for PostState {
    fn default() -> Self {
        PostState { obs_var: 0, obs_sel: vec![], prof_var: 0, prof_sel: vec![], flux_var: 0, solute: 0, run_var: 0, prop_mat: 0, prop_kind: 0 }
    }
}

pub struct App {
    pub prj: Project,
    pub path: Option<PathBuf>,
    pub dirty: bool,
    pub page: Page,
    pub results: Option<Results>,
    run: Option<RunHandle>,
    pub status: String,
    pub dialog: Option<FileDialog>,
    pub problems: Vec<String>,
    pub sel_mat: usize,
    pub edit: EditState,
    pub post: PostState,
    pub print_text: String,
    pub print_text_for: Vec<f64>,
    pub show_about: bool,
    pub last_dir: Option<PathBuf>,
    pub run_seconds: f64,
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        cc.egui_ctx.set_visuals(egui::Visuals::dark());
        let mut style = (*cc.egui_ctx.style()).clone();
        style.spacing.item_spacing = egui::vec2(8.0, 6.0);
        style.spacing.button_padding = egui::vec2(8.0, 4.0);
        cc.egui_ctx.set_style(style);
        let prj = hydrus_core::examples::infiltration_evaporation();
        let mut app = App {
            prj,
            path: None,
            dirty: false,
            page: Page::Main,
            results: None,
            run: None,
            status: "Ready. Load an example, import a HYDRUS-1D project folder, or edit the current project.".into(),
            dialog: None,
            problems: vec![],
            sel_mat: 0,
            edit: EditState { to: 100.0, mat: 1, h_top: -100.0, h_bot: -100.0, beta: 1.0, temp: 20.0, n_nodes: 101, fixed: vec![(0.0, 1.0, 1.0), (100.0, 4.0, 4.0)], ..Default::default() },
            post: PostState::default(),
            print_text: String::new(),
            print_text_for: vec![],
            show_about: false,
            last_dir: None,
            run_seconds: 0.0,
        };
        app.edit.n_nodes = app.prj.profile.nodes.len();
        // optional command-line argument: project file (.h1dr) or legacy HYDRUS-1D folder
        if let Some(arg) = std::env::args().nth(1) {
            let p = PathBuf::from(&arg);
            let r = if p.is_dir() { hydrus_io::read_legacy_project(&p) } else { hydrus_io::load_project(&p) };
            match r {
                Ok(prj) => {
                    let keep = if p.is_dir() { None } else { Some(p.clone()) };
                    app.set_project(prj, keep);
                    app.status = format!("Opened {}", p.display());
                }
                Err(e) => app.status = format!("Could not open {}: {}", arg, e),
            }
        }
        app
    }

    pub fn set_project(&mut self, p: Project, path: Option<PathBuf>) {
        self.prj = p;
        self.path = path;
        self.dirty = false;
        self.results = None;
        self.sel_mat = 0;
        self.problems.clear();
        self.print_text_for.clear();
        self.edit.n_nodes = self.prj.profile.nodes.len();
        self.edit.to = self.prj.profile.depth();
        self.post = PostState::default();
        self.page = Page::Main;
		self.sync();
    }

    fn start_run(&mut self) {
        self.problems = self.prj.validate();
        if !self.problems.is_empty() {
            self.status = "The project has problems – see the Summary page.".into();
            self.page = Page::Summary;
            return;
        }
        let prj = self.prj.clone();
        let (tx, rx) = channel();
        let progress = Arc::new(AtomicU64::new(prj.time.t_init.to_bits()));
        let cancel = Arc::new(AtomicBool::new(false));
        let (p2, c2) = (progress.clone(), cancel.clone());
        let (t_init, t_max) = (prj.time.t_init, prj.time.t_max);
        std::thread::spawn(move || {
            let res = match Simulation::new(prj) {
                Ok(mut sim) => {
                    sim.run(|s| {
                        p2.store(s.t.to_bits(), Ordering::Relaxed);
                        !c2.load(Ordering::Relaxed)
                    });
                    Ok(std::mem::take(&mut sim.res))
                }
                Err(e) => Err(e.to_string()),
            };
            let _ = tx.send(res);
        });
        self.run = Some(RunHandle { rx, progress, cancel, t_init, t_max, started: Instant::now() });
        self.status = "Running…".into();
    }

    fn poll_run(&mut self, ctx: &egui::Context) {
        let mut finished = None;
        if let Some(r) = &self.run {
            ctx.request_repaint_after(std::time::Duration::from_millis(80));
            if let Ok(res) = r.rx.try_recv() {
                finished = Some((res, r.started.elapsed().as_secs_f64(), r.cancel.load(Ordering::Relaxed)));
            }
        }
        if let Some((res, secs, cancelled)) = finished {
            self.run = None;
            self.run_seconds = secs;
            match res {
                Ok(r) => {
                    let msg = if r.failed {
                        format!("Calculation stopped: {}", r.messages.join(" "))
                    } else if cancelled {
                        "Calculation cancelled (partial results).".to_string()
                    } else {
                        format!("Calculation finished in {:.2} s – {} time levels.", secs, r.tlevel.len())
                    };
                    self.status = msg;
                    self.post.obs_sel.clear();
                    self.post.prof_sel.clear();
                    self.results = Some(r);
                    if self.page == Page::Main {
                        self.page = Page::ProfileInfo;
                    }
                }
                Err(e) => self.status = format!("Error: {}", e),
            }
        }
    }

    fn menu(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("New project").clicked() {
                        let mut p = Project::default();
                        p.profile = hydrus_core::examples::infiltration_evaporation().profile;
                        for n in p.profile.nodes.iter_mut() {
                            n.beta = 0.0;
                            n.mat = 1;
                            n.layer = 1;
                            n.h = -100.0;
                        }
                        self.set_project(p, None);
                        ui.close_menu();
                    }
                    ui.menu_button("Examples", |ui| {
                        for (name, f) in hydrus_core::examples::all() {
                            if ui.button(name).clicked() {
                                self.set_project(f(), None);
                                ui.close_menu();
                            }
                        }
                    });
                    ui.separator();
                    if ui.button("Open project (.h1dr)…").clicked() {
                        self.dialog = Some(FileDialog::new(Purpose::OpenProject, Mode::OpenFile, "Open project", self.last_dir.clone(), "", &["h1dr", "json"]));
                        ui.close_menu();
                    }
                    if ui.button("Import HYDRUS-1D project folder…").clicked() {
                        self.dialog = Some(FileDialog::new(Purpose::ImportLegacy, Mode::PickDir, "Select a HYDRUS-1D project folder (contains Selector.in)", self.last_dir.clone(), "", &[]));
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.add_enabled(self.path.is_some(), egui::Button::new("Save")).clicked() {
                        self.save();
                        ui.close_menu();
                    }
                    if ui.button("Save as…").clicked() {
                        let name = self.path.as_ref().and_then(|p| p.file_name()).map(|s| s.to_string_lossy().to_string()).unwrap_or_else(|| "project.h1dr".into());
                        self.dialog = Some(FileDialog::new(Purpose::SaveProject, Mode::SaveFile, "Save project", self.last_dir.clone(), &name, &["h1dr"]));
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.add_enabled(self.results.is_some(), egui::Button::new("Export results (CSV)…")).clicked() {
                        self.dialog = Some(FileDialog::new(Purpose::ExportCsv, Mode::PickDir, "Choose a folder for the CSV files", self.last_dir.clone(), "", &[]));
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
                ui.menu_button("Calculation", |ui| {
                    if ui.add_enabled(self.run.is_none(), egui::Button::new("▶ Run   (F5)")).clicked() {
                        self.start_run();
                        ui.close_menu();
                    }
                    if ui.add_enabled(self.run.is_some(), egui::Button::new("■ Stop")).clicked() {
                        if let Some(r) = &self.run {
                            r.cancel.store(true, Ordering::Relaxed);
                        }
                        ui.close_menu();
                    }
                });
                ui.menu_button("View", |ui| {
                    if ui.button("Dark theme").clicked() {
                        ctx.set_visuals(egui::Visuals::dark());
                        ui.close_menu();
                    }
                    if ui.button("Light theme").clicked() {
                        ctx.set_visuals(egui::Visuals::light());
                        ui.close_menu();
                    }
                });
                ui.menu_button("Help", |ui| {
                    if ui.button("About").clicked() {
                        self.show_about = true;
                        ui.close_menu();
                    }
                });
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let running = self.run.is_some();
                    if running {
                        if ui.button("■ Stop").clicked() {
                            if let Some(r) = &self.run {
                                r.cancel.store(true, Ordering::Relaxed);
                            }
                        }
                    } else if ui.button("▶ Run").clicked() {
                        self.start_run();
                    }
                    let title = format!("{}{}", self.prj.title, if self.dirty { " *" } else { "" });
                    ui.label(egui::RichText::new(title).weak());
                });
            });
        });
        if ctx.input(|i| i.key_pressed(egui::Key::F5)) && self.run.is_none() {
            self.start_run();
        }
    }

    fn save(&mut self) {
        if let Some(p) = self.path.clone() {
            match hydrus_io::save_project(&self.prj, &p) {
                Ok(_) => {
                    self.dirty = false;
                    self.status = format!("Saved {}", p.display());
                }
                Err(e) => self.status = e.to_string(),
            }
        }
    }

    fn dialogs(&mut self, ctx: &egui::Context) {
        if let Some(mut d) = self.dialog.take() {
            match d.show(ctx) {
                DialogResult::Open => self.dialog = Some(d),
                DialogResult::Cancel => {}
                DialogResult::Chosen(path) => {
                    self.last_dir = if path.is_dir() { Some(path.clone()) } else { path.parent().map(|p| p.to_path_buf()) };
                    match d.purpose {
                        Purpose::OpenProject => match hydrus_io::load_project(&path) {
                            Ok(p) => {
                                self.set_project(p, Some(path.clone()));
                                self.status = format!("Opened {}", path.display());
                            }
                            Err(e) => self.status = e.to_string(),
                        },
                        Purpose::ImportLegacy => {
                            match hydrus_io::read_legacy_project(&path) {
                                Ok(mut p) => {
                                    if p.title.trim().is_empty() {
                                        p.title = path
                                            .file_name()
                                            .map(|s| s.to_string_lossy().to_string())
                                            .unwrap_or_default();
                                    }
                                    self.set_project(p, None);
                                    self.dirty = true;
                                    self.status = format!("Imported HYDRUS-1D project from {}", path.display());
                                }
                                Err(e) => {
                                    self.status = format!("Import failed: {}", e);
                                    self.problems = vec![format!("Failed to import from {}: {}", path.display(), e)];
                                }
                            }
                        },
                        Purpose::SaveProject => {
                            let mut p = path.clone();
                            if p.extension().is_none() {
                                p.set_extension("h1dr");
                            }
                            self.path = Some(p);
                            self.save();
                        }
                        Purpose::ExportCsv => {
                            if let Some(r) = &self.results {
                                self.status = match hydrus_io::writer::write_all(r, &path) {
                                    Ok(_) => format!("Results written to {}", path.display()),
                                    Err(e) => e.to_string(),
                                };
                            }
                        }
                        Purpose::ImportAtmosphere => {
                            self.status = match self.import_atmosphere_text(&path) {
                                Ok(n) => format!("Imported {} atmospheric records", n),
                                Err(e) => e,
                            };
                        }
                    }
                }
            }
        }
        if self.show_about {
            egui::Window::new("About").open(&mut self.show_about).collapsible(false).show(ctx, |ui| {
                ui.heading("HYDRUS-1D (Rust port)");
                ui.label("Numerical engine ported from the HYDRUS-1D Fortran source:\nvariably saturated water flow (Richards equation), root water uptake,\nsolute transport with sorption and decay chains, and heat transport.");
                ui.label("Validated against the original code (see README).");
                ui.hyperlink_to("Original HYDRUS-1D: Šimůnek, van Genuchten & Šejna (PC-Progress)", "https://www.pc-progress.com/en/Default.aspx?hydrus-1d");
            });
        }
    }

    fn nav(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("nav").resizable(true).default_width(250.0).show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.add_space(4.0);
                ui.label(egui::RichText::new("Pre-processing").strong().size(16.0));
                ui.separator();
                let pr = self.prj.processes.clone();
                let mut items: Vec<(Page, &str)> = vec![
                    (Page::Main, "Main processes & units"),
                    (Page::Geometry, "Geometry & materials"),
                    (Page::Profile, "Soil profile editor"),
                    (Page::NodeTable, "Nodal table"),
                    (Page::Time, "Time & printing"),
                    (Page::Iteration, "Water flow: iteration criteria"),
                    (Page::HydModel, "Water flow: hydraulic model"),
                    (Page::SoilParams, "Water flow: soil parameters"),
                    (Page::WaterBc, "Water flow: boundary conditions"),
                    (Page::Atmosphere, "Atmospheric data"),
                ];
				if self.prj.atmosphere.meteo.is_some() {
					items.push((Page::Meteo, "Meteorological parameters"));
				}
                if pr.root_water_uptake {
                    items.push((Page::RootUptake, "Root water uptake"));
                }
                if pr.solute {
                    items.push((Page::SoluteGeneral, "Solute: general"));
                    items.push((Page::SoluteMaterials, "Solute: transport parameters"));
                    items.push((Page::SoluteReactions, "Solute: reaction parameters"));
                    items.push((Page::SoluteBc, "Solute: boundary conditions"));
                }
                if pr.heat {
                    items.push((Page::Heat, "Heat transport"));
                }
                items.push((Page::Summary, "Check & summary"));
                for (p, label) in items {
                    if ui.selectable_label(self.page == p, label).clicked() {
                        self.page = p;
                    }
                }
                ui.add_space(12.0);
                ui.label(egui::RichText::new("Post-processing").strong().size(16.0));
                ui.separator();
                let has = self.results.is_some();
                for (p, label) in [
                    (Page::Obs, "Observation points"),
                    (Page::ProfileInfo, "Profile information"),
                    (Page::Fluxes, "Boundary fluxes & heads"),
                    (Page::SoilProps, "Soil hydraulic properties"),
                    (Page::RunTime, "Run-time information"),
                    (Page::MassBalance, "Mass balance"),
                ] {
                    let enabled = has || p == Page::SoilProps;
                    if ui.add_enabled(enabled, egui::SelectableLabel::new(self.page == p, label)).clicked() {
                        self.page = p;
                    }
                }
            });
        });
    }

    fn status_bar(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if let Some(r) = &self.run {
                    let t = f64::from_bits(r.progress.load(Ordering::Relaxed));
                    let f = ((t - r.t_init) / (r.t_max - r.t_init)).clamp(0.0, 1.0) as f32;
                    ui.add(egui::ProgressBar::new(f).desired_width(220.0).show_percentage());
                    ui.label(format!("t = {:.4}", t));
                }
                ui.label(&self.status);
            });
        });
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_run(ctx);
        self.menu(ctx);
        self.status_bar(ctx);
        self.nav(ctx);
        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
                self.page_ui(ui);
            });
        });
        self.dialogs(ctx);
    }
}

impl App {
    fn page_ui(&mut self, ui: &mut egui::Ui) {
        let before = self.dirty;
        let mut ch = false;
        match self.page {
            Page::Main => self.page_main(ui, &mut ch),
            Page::Geometry => self.page_geometry(ui, &mut ch),
            Page::Profile => self.page_profile(ui, &mut ch),
            Page::NodeTable => self.page_nodes(ui, &mut ch),
            Page::Time => self.page_time(ui, &mut ch),
            Page::Iteration => self.page_iteration(ui, &mut ch),
            Page::HydModel => self.page_hydmodel(ui, &mut ch),
            Page::SoilParams => self.page_soil(ui, &mut ch),
            Page::WaterBc => self.page_waterbc(ui, &mut ch),
            Page::Atmosphere => self.page_atmosphere(ui, &mut ch),
			Page::Meteo => self.page_meteo(ui, &mut ch),
            Page::RootUptake => self.page_root(ui, &mut ch),
            Page::SoluteGeneral => self.page_sol_general(ui, &mut ch),
            Page::SoluteMaterials => self.page_sol_materials(ui, &mut ch),
            Page::SoluteReactions => self.page_sol_reactions(ui, &mut ch),
            Page::SoluteBc => self.page_sol_bc(ui, &mut ch),
            Page::Heat => self.page_heat(ui, &mut ch),
            Page::Summary => self.page_summary(ui),
            Page::Obs => self.page_obs(ui),
            Page::ProfileInfo => self.page_profile_info(ui),
            Page::Fluxes => self.page_fluxes(ui),
            Page::SoilProps => self.page_soil_props(ui),
            Page::RunTime => self.page_runtime(ui),
            Page::MassBalance => self.page_balance(ui),
        }
        if ch {
            self.dirty = true;
        }
        let _ = before;
    }
}

fn main() -> eframe::Result<()> {
    let opts = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([1440.0, 900.0]).with_min_inner_size([900.0, 600.0]).with_title("HYDRUS-1D (Rust)"),
        ..Default::default()
    };
    eframe::run_native("HYDRUS-1D", opts, Box::new(|cc| Ok(Box::new(App::new(cc)))))
}
