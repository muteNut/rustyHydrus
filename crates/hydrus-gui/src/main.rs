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
    Profile,
    NodeTable,
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
        PostState {
            obs_var: 0,
            obs_sel: vec![],
            prof_var: 0,
            prof_sel: vec![],
            flux_var: 0,
            solute: 0,
            run_var: 0,
            prop_mat: 0,
            prop_kind: 0,
        }
    }
}

fn load_system_font(ctx: &egui::Context) {
    let mut fonts = egui::FontDefinitions::default();

    let candidate_paths: &[&str] = &[
        // Linux / Ubuntu / GNOME
        "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf",
        "/usr/share/fonts/opentype/cantarell/Cantarell-Regular.otf",
        "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/TTF/DejaVuSans.ttf",
        // Windows
        "C:\\Windows\\Fonts\\segoeui.ttf",
        "C:\\Windows\\Fonts\\arial.ttf",
        // macOS
        "/System/Library/Fonts/SFProText-Regular.otf",
        "/System/Library/Fonts/HelveticaNeue.ttc",
        "/Library/Fonts/Arial.ttf",
    ];

    for path in candidate_paths {
        if let Ok(data) = std::fs::read(path) {
            fonts.font_data.insert(
                "os_default".to_owned(),
                std::sync::Arc::new(egui::FontData::from_owned(data)),
            );

            if let Some(family) = fonts.families.get_mut(&egui::FontFamily::Proportional) {
                family.insert(0, "os_default".to_owned());
            }
            break;
        }
    }

    ctx.set_fonts(fonts);
}

pub fn apply_custom_dark_theme(ctx: &egui::Context) {
    let mut visuals = egui::Visuals::dark();
    visuals.panel_fill = egui::Color32::from_rgb(26, 28, 32);
    visuals.window_fill = egui::Color32::from_rgb(32, 34, 39);
    visuals.faint_bg_color = egui::Color32::from_rgb(38, 41, 48);
    visuals.extreme_bg_color = egui::Color32::from_rgb(18, 19, 22);

    visuals.window_corner_radius = egui::CornerRadius::same(8);
    visuals.widgets.noninteractive.corner_radius = egui::CornerRadius::same(4);
    visuals.widgets.inactive.corner_radius = egui::CornerRadius::same(5);
    visuals.widgets.hovered.corner_radius = egui::CornerRadius::same(5);
    visuals.widgets.active.corner_radius = egui::CornerRadius::same(5);

    visuals.widgets.noninteractive.bg_stroke =
        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(45, 48, 56));
    visuals.widgets.inactive.bg_stroke =
        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(52, 56, 65));
    visuals.widgets.hovered.bg_stroke =
        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(80, 140, 220));

    ctx.set_visuals(visuals);
}

pub fn apply_custom_light_theme(ctx: &egui::Context) {
    let mut visuals = egui::Visuals::light();
    visuals.panel_fill = egui::Color32::from_rgb(245, 246, 248);
    visuals.window_fill = egui::Color32::from_rgb(255, 255, 255);
    visuals.faint_bg_color = egui::Color32::from_rgb(238, 240, 244);
    visuals.extreme_bg_color = egui::Color32::from_rgb(255, 255, 255);

    visuals.window_corner_radius = egui::CornerRadius::same(8);
    visuals.widgets.noninteractive.corner_radius = egui::CornerRadius::same(4);
    visuals.widgets.inactive.corner_radius = egui::CornerRadius::same(5);
    visuals.widgets.hovered.corner_radius = egui::CornerRadius::same(5);
    visuals.widgets.active.corner_radius = egui::CornerRadius::same(5);

    visuals.widgets.noninteractive.bg_stroke =
        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(215, 220, 228));
    visuals.widgets.inactive.bg_stroke =
        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(200, 205, 215));
    visuals.widgets.hovered.bg_stroke =
        egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(60, 120, 210));

    ctx.set_visuals(visuals);
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
    pub first_frame: bool,
    pub show_close_confirm: bool,
    pub pending_action: Option<PendingAction>,
}

#[derive(Clone, Copy, Debug)]
pub enum PendingAction {
    Quit,
    NewProject,
    OpenProject,
}

impl App {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        load_system_font(&cc.egui_ctx);

        #[cfg(target_os = "linux")]
        {
            if cc.egui_ctx.pixels_per_point() < 1.1 {
                cc.egui_ctx.set_pixels_per_point(1.15);
            }
        }

        apply_custom_dark_theme(&cc.egui_ctx);

        let mut style = (*cc.egui_ctx.style()).clone();
        style.spacing.item_spacing = egui::vec2(10.0, 8.0);
        style.spacing.button_padding = egui::vec2(10.0, 6.0);
        style.spacing.window_margin = egui::Margin::same(12);
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
            edit: EditState {
                to: 100.0,
                mat: 1,
                h_top: -100.0,
                h_bot: -100.0,
                beta: 1.0,
                temp: 20.0,
                n_nodes: 101,
                fixed: vec![(0.0, 1.0, 1.0), (100.0, 4.0, 4.0)],
                ..Default::default()
            },
            post: PostState::default(),
            print_text: String::new(),
            print_text_for: vec![],
            show_about: false,
            last_dir: None,
            run_seconds: 0.0,
            first_frame: true,
            show_close_confirm: false,
            pending_action: None,
        };
        app.edit.n_nodes = app.prj.profile.nodes.len();
        if let Some(arg) = std::env::args().nth(1) {
            let p = PathBuf::from(&arg);
            let r = if p.is_dir() {
                hydrus_io::read_legacy_project(&p)
            } else {
                hydrus_io::load_project(&p)
            };
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
        self.run = Some(RunHandle {
            rx,
            progress,
            cancel,
            t_init,
            t_max,
            started: Instant::now(),
        });
        self.status = "Running…".into();
    }

    fn poll_run(&mut self, ctx: &egui::Context) {
        let mut finished = None;
        if let Some(r) = &self.run {
            ctx.request_repaint_after(std::time::Duration::from_millis(80));
            if let Ok(res) = r.rx.try_recv() {
                finished = Some((
                    res,
                    r.started.elapsed().as_secs_f64(),
                    r.cancel.load(Ordering::Relaxed),
                ));
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
                        format!(
                            "Calculation finished in {:.2} s – {} time levels.",
                            secs,
                            r.tlevel.len()
                        )
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

    fn request_new_project(&mut self) {
        if self.dirty {
            self.pending_action = Some(PendingAction::NewProject);
            self.show_close_confirm = true;
        } else {
            self.do_new_project();
        }
    }

    fn do_new_project(&mut self) {
        let mut p = Project::default();
        p.profile = hydrus_core::examples::infiltration_evaporation().profile;
        for n in p.profile.nodes.iter_mut() {
            n.beta = 0.0;
            n.mat = 1;
            n.layer = 1;
            n.h = -100.0;
        }
        self.set_project(p, None);
    }

    fn request_open_project(&mut self) {
        if self.dirty {
            self.pending_action = Some(PendingAction::OpenProject);
            self.show_close_confirm = true;
        } else {
            self.dialog = Some(FileDialog::new(
                Purpose::OpenProject,
                Mode::OpenFile,
                "Open project",
                self.last_dir.clone(),
                "",
                &["h1dr", "json"],
            ));
        }
    }

    fn request_quit(&mut self, ctx: &egui::Context) {
        if self.dirty {
            self.pending_action = Some(PendingAction::Quit);
            self.show_close_confirm = true;
        } else {
            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
        }
    }

    fn handle_shortcuts(&mut self, ctx: &egui::Context) {
        ctx.input(|i| {
            let ctrl = i.modifiers.command || i.modifiers.ctrl;
            if ctrl && !i.modifiers.shift && i.key_pressed(egui::Key::N) {
                self.request_new_project();
            }
            if ctrl && !i.modifiers.shift && i.key_pressed(egui::Key::O) {
                self.request_open_project();
            }
            if ctrl && i.modifiers.shift && i.key_pressed(egui::Key::S) {
                let name = self
                    .path
                    .as_ref()
                    .and_then(|p| p.file_name())
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_else(|| "project.h1dr".into());
                self.dialog = Some(FileDialog::new(
                    Purpose::SaveProject,
                    Mode::SaveFile,
                    "Save project",
                    self.last_dir.clone(),
                    &name,
                    &["h1dr"],
                ));
            } else if ctrl && !i.modifiers.shift && i.key_pressed(egui::Key::S) {
                if self.path.is_some() {
                    self.save();
                } else {
                    self.dialog = Some(FileDialog::new(
                        Purpose::SaveProject,
                        Mode::SaveFile,
                        "Save project",
                        self.last_dir.clone(),
                        "project.h1dr",
                        &["h1dr"],
                    ));
                }
            }
        });
    }

    fn menu(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("New project   (Ctrl+N)").clicked() {
                        ui.close_menu();
                        self.request_new_project();
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
                    if ui.button("Open project (.h1dr)…   (Ctrl+O)").clicked() {
                        ui.close_menu();
                        self.request_open_project();
                    }
                    if ui.button("Import HYDRUS-1D project folder…").clicked() {
                        self.dialog = Some(FileDialog::new(
                            Purpose::ImportLegacy,
                            Mode::PickDir,
                            "Select a HYDRUS-1D project folder (contains Selector.in)",
                            self.last_dir.clone(),
                            "",
                            &[],
                        ));
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui
                        .add_enabled(self.path.is_some(), egui::Button::new("Save   (Ctrl+S)"))
                        .clicked()
                    {
                        self.save();
                        ui.close_menu();
                    }
                    if ui.button("Save as…   (Ctrl+Shift+S)").clicked() {
                        let name = self
                            .path
                            .as_ref()
                            .and_then(|p| p.file_name())
                            .map(|s| s.to_string_lossy().to_string())
                            .unwrap_or_else(|| "project.h1dr".into());
                        self.dialog = Some(FileDialog::new(
                            Purpose::SaveProject,
                            Mode::SaveFile,
                            "Save project",
                            self.last_dir.clone(),
                            &name,
                            &["h1dr"],
                        ));
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui
                        .add_enabled(
                            self.results.is_some(),
                            egui::Button::new("Export results (CSV)…"),
                        )
                        .clicked()
                    {
                        self.dialog = Some(FileDialog::new(
                            Purpose::ExportCsv,
                            Mode::PickDir,
                            "Choose a folder for the CSV files",
                            self.last_dir.clone(),
                            "",
                            &[],
                        ));
                        ui.close_menu();
                    }
                    ui.separator();
                    if ui.button("Quit").clicked() {
                        ui.close_menu();
                        self.request_quit(ctx);
                    }
                });
                ui.menu_button("Calculation", |ui| {
                    if ui
                        .add_enabled(self.run.is_none(), egui::Button::new("▶ Run   (F5)"))
                        .clicked()
                    {
                        self.start_run();
                        ui.close_menu();
                    }
                    if ui
                        .add_enabled(self.run.is_some(), egui::Button::new("■ Stop"))
                        .clicked()
                    {
                        if let Some(r) = &self.run {
                            r.cancel.store(true, Ordering::Relaxed);
                        }
                        ui.close_menu();
                    }
                });
                ui.menu_button("View", |ui| {
                    if ui.button("Dark theme").clicked() {
                        apply_custom_dark_theme(ctx);
                        ui.close_menu();
                    }
                    if ui.button("Light theme").clicked() {
                        apply_custom_light_theme(ctx);
                        ui.close_menu();
                    }
                });
                ui.menu_button("Help", |ui| {
                    if ui.button("About Rusteau-1D").clicked() {
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
        if self.show_close_confirm {
            egui::Window::new("Unsaved Changes")
                .collapsible(false)
                .resizable(false)
                .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
                .show(ctx, |ui| {
                    ui.label("The current project has unsaved changes.\nDo you want to discard them and continue?");
                    ui.add_space(8.0);
                    ui.horizontal(|ui| {
                        if ui.button("Discard & Continue").clicked() {
                            self.show_close_confirm = false;
                            self.dirty = false;
                            match self.pending_action.take() {
                                Some(PendingAction::Quit) => {
                                    ctx.send_viewport_cmd(egui::ViewportCommand::Close)
                                }
                                Some(PendingAction::NewProject) => self.do_new_project(),
                                Some(PendingAction::OpenProject) => {
                                    self.dialog = Some(FileDialog::new(
                                        Purpose::OpenProject,
                                        Mode::OpenFile,
                                        "Open project",
                                        self.last_dir.clone(),
                                        "",
                                        &["h1dr", "json"],
                                    ));
                                }
                                None => {}
                            }
                        }
                        if ui.button("Cancel").clicked() {
                            self.show_close_confirm = false;
                            self.pending_action = None;
                        }
                    });
                });
        }

        if let Some(mut d) = self.dialog.take() {
            match d.show(ctx) {
                DialogResult::Open => self.dialog = Some(d),
                DialogResult::Cancel => {}
                DialogResult::Chosen(path) => {
                    self.last_dir = if path.is_dir() {
                        Some(path.clone())
                    } else {
                        path.parent().map(|p| p.to_path_buf())
                    };
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
                                    self.status = format!(
                                        "Imported HYDRUS-1D project from {}",
                                        path.display()
                                    );
                                }
                                Err(e) => {
                                    self.status = format!("Import failed: {}", e);
                                    self.problems = vec![format!(
                                        "Failed to import from {}: {}",
                                        path.display(),
                                        e
                                    )];
                                }
                            }
                        }
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
            egui::Window::new("About Rusteau-1D")
                .open(&mut self.show_about)
                .collapsible(false)
                .show(ctx, |ui| {
                    ui.heading("Rusteau-1D");
                    ui.label("Modern, performant porous media & vadose zone simulator.\nNumerical engine ported and verified against HYDRUS-1D Fortran routines:\nvariably saturated flow (Richards), solute transport, and heat flow.");
                    ui.add_space(4.0);
                    ui.hyperlink_to(
                        "Original HYDRUS-1D: Šimůnek, van Genuchten & Šejna",
                        "https://www.pc-progress.com/en/Default.aspx?hydrus-1d",
                    );
                });
        }
    }

    fn nav(&mut self, ctx: &egui::Context) {
        egui::SidePanel::left("nav")
            .resizable(true)
            .default_width(260.0)
            .min_width(220.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Pre-processing").strong().size(16.0));
                    ui.separator();
                    let pr = self.prj.processes.clone();
                    let mut items: Vec<(Page, &str)> = vec![
                        (Page::Main, "Main processes & units"),
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
                    items.push((Page::Profile, "Soil profile editor"));
                    items.push((Page::NodeTable, "Nodal table"));
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
                        if ui
                            .add_enabled(enabled, egui::SelectableLabel::new(self.page == p, label))
                            .clicked()
                        {
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
        if self.first_frame {
            self.first_frame = false;
            let screen_rect = ctx.screen_rect();
            let screen_w = screen_rect.width();
            let screen_h = screen_rect.height();
            let target_w = (screen_w * 0.80).clamp(960.0, 1920.0);
            let target_h = (screen_h * 0.82).clamp(620.0, 1200.0);
            ctx.send_viewport_cmd(egui::ViewportCommand::InnerSize(egui::vec2(
                target_w, target_h,
            )));
        }

        if ctx.input(|i| i.viewport().close_requested()) && self.dirty {
            ctx.send_viewport_cmd(egui::ViewportCommand::CancelClose);
            self.pending_action = Some(PendingAction::Quit);
            self.show_close_confirm = true;
        }

        self.handle_shortcuts(ctx);
        self.poll_run(ctx);
        self.menu(ctx);
        self.status_bar(ctx);
        self.nav(ctx);
        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::both()
                .auto_shrink([false, false])
                .show(ui, |ui| {
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
            Page::Profile => self.page_profile(ui, &mut ch),
            Page::NodeTable => self.page_nodes(ui, &mut ch),
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
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1280.0, 800.0])
            .with_min_inner_size([920.0, 600.0])
            .with_resizable(true)
            .with_title("Rusteau-1D"),
        ..Default::default()
    };
    eframe::run_native("Rusteau-1D", opts, Box::new(|cc| Ok(Box::new(App::new(cc)))))
}