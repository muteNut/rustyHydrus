//! A small self-contained file/folder browser (works identically on Linux and Windows
//! without native dialog dependencies).

use egui::{Context, ScrollArea, Window};
use std::path::PathBuf;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mode {
    OpenFile,
    SaveFile,
    PickDir,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Purpose {
    OpenProject,
    SaveProject,
    ImportLegacy,
    ExportCsv,
    ImportAtmosphere,
}

pub enum DialogResult {
    Open,
    Cancel,
    Chosen(PathBuf),
}

pub struct FileDialog {
    pub purpose: Purpose,
    pub mode: Mode,
    pub title: String,
    pub dir: PathBuf,
    pub name: String,
    pub ext: Vec<String>,
    dir_text: String,
    entries: Vec<(String, bool)>,
    err: String,
}

pub fn home_dir() -> PathBuf {
    std::env::var_os("HOME")
        .or_else(|| std::env::var_os("USERPROFILE"))
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."))
}

impl FileDialog {
    pub fn new(purpose: Purpose, mode: Mode, title: &str, start: Option<PathBuf>, name: &str, ext: &[&str]) -> Self {
        let dir = start.filter(|p| p.is_dir()).unwrap_or_else(|| std::env::current_dir().unwrap_or_else(|_| home_dir()));
        let mut d = FileDialog {
            purpose,
            mode,
            title: title.to_string(),
            dir_text: dir.display().to_string(),
            dir,
            name: name.to_string(),
            ext: ext.iter().map(|s| s.to_string()).collect(),
            entries: vec![],
            err: String::new(),
        };
        d.refresh();
        d
    }

    fn refresh(&mut self) {
        self.entries.clear();
        self.err.clear();
        match std::fs::read_dir(&self.dir) {
            Ok(rd) => {
                for e in rd.flatten() {
                    let name = e.file_name().to_string_lossy().to_string();
                    if name.starts_with('.') {
                        continue;
                    }
                    let is_dir = e.path().is_dir();
                    if !is_dir && self.mode == Mode::PickDir {
                        continue;
                    }
                    if !is_dir && !self.ext.is_empty() {
                        let lower = name.to_ascii_lowercase();
                        if !self.ext.iter().any(|x| lower.ends_with(&format!(".{}", x.to_ascii_lowercase()))) {
                            continue;
                        }
                    }
                    self.entries.push((name, is_dir));
                }
                self.entries.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.to_lowercase().cmp(&b.0.to_lowercase())));
            }
            Err(e) => self.err = e.to_string(),
        }
        self.dir_text = self.dir.display().to_string();
    }

    fn go(&mut self, p: PathBuf) {
        if p.is_dir() {
            self.dir = p;
            self.refresh();
        }
    }

    pub fn show(&mut self, ctx: &Context) -> DialogResult {
        let mut result = DialogResult::Open;
        let mut open = true;
        let mut nav: Option<PathBuf> = None;
        Window::new(self.title.clone()).open(&mut open).collapsible(false).default_size([620.0, 460.0]).show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui.button("⬆ Up").clicked() {
                    if let Some(p) = self.dir.parent() {
                        nav = Some(p.to_path_buf());
                    }
                }
                if ui.button("🏠 Home").clicked() {
                    nav = Some(home_dir());
                }
                #[cfg(windows)]
                {
                    for c in b'C'..=b'Z' {
                        let d = format!("{}:\\", c as char);
                        if std::path::Path::new(&d).exists() && ui.button(format!("{}:", c as char)).clicked() {
                            nav = Some(PathBuf::from(d));
                        }
                    }
                }
                let r = ui.add(egui::TextEdit::singleline(&mut self.dir_text).desired_width(f32::INFINITY));
                if r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                    nav = Some(PathBuf::from(self.dir_text.clone()));
                }
            });
            ui.separator();
            ScrollArea::vertical().max_height(300.0).auto_shrink([false, false]).show(ui, |ui| {
                for (name, is_dir) in self.entries.clone() {
                    let label = if is_dir { format!("📁 {}", name) } else { format!("📄 {}", name) };
                    let r = ui.selectable_label(self.name == name && !is_dir, label);
                    if r.clicked() {
                        if is_dir {
                            if self.mode == Mode::PickDir {
                                self.name = name.clone();
                            }
                        } else {
                            self.name = name.clone();
                        }
                    }
                    if r.double_clicked() && is_dir {
                        nav = Some(self.dir.join(&name));
                    }
                }
                if !self.err.is_empty() {
                    ui.colored_label(egui::Color32::LIGHT_RED, &self.err);
                }
            });
            ui.separator();
            ui.horizontal(|ui| {
                match self.mode {
                    Mode::PickDir => {
                        ui.label("Selected folder:");
                        let sel = if self.name.is_empty() { self.dir.clone() } else { self.dir.join(&self.name) };
                        ui.monospace(sel.display().to_string());
                    }
                    _ => {
                        ui.label("File name:");
                        ui.add(egui::TextEdit::singleline(&mut self.name).desired_width(280.0));
                    }
                }
            });
            ui.horizontal(|ui| {
                let label = match self.mode {
                    Mode::OpenFile => "Open",
                    Mode::SaveFile => "Save",
                    Mode::PickDir => "Select this folder",
                };
                if ui.button(label).clicked() {
                    let p = match self.mode {
                        Mode::PickDir => {
                            if self.name.is_empty() {
                                self.dir.clone()
                            } else {
                                self.dir.join(&self.name)
                            }
                        }
                        _ => self.dir.join(&self.name),
                    };
                    if self.mode == Mode::OpenFile && !p.is_file() {
                        self.err = "Please select an existing file.".into();
                    } else {
                        result = DialogResult::Chosen(p);
                    }
                }
                if ui.button("Cancel").clicked() {
                    result = DialogResult::Cancel;
                }
            });
        });
        if let Some(p) = nav {
            self.name.clear();
            self.go(p);
        }
        if !open {
            return DialogResult::Cancel;
        }
        result
    }
}
