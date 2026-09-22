use egui::{Color32, DragValue, Ui};
use hydrus_core::{LengthUnit, Project, SoilMaterial, TimeUnit};

pub const PALETTE: [Color32; 10] = [
    Color32::from_rgb(31, 119, 180),
    Color32::from_rgb(255, 127, 14),
    Color32::from_rgb(44, 160, 44),
    Color32::from_rgb(214, 39, 40),
    Color32::from_rgb(148, 103, 189),
    Color32::from_rgb(140, 86, 75),
    Color32::from_rgb(227, 119, 194),
    Color32::from_rgb(127, 127, 127),
    Color32::from_rgb(188, 189, 34),
    Color32::from_rgb(23, 190, 207),
];

pub fn mat_color(m: usize) -> Color32 {
    PALETTE[(m.max(1) - 1) % PALETTE.len()]
}

/// USDA soil textural classes, Carsel & Parrish (1988) as in the HYDRUS soil catalog (cm, day).
pub const SOIL_CATALOG: [(&str, f64, f64, f64, f64, f64); 12] = [
    ("Sand", 0.045, 0.43, 0.145, 2.68, 712.8),
    ("Loamy Sand", 0.057, 0.41, 0.124, 2.28, 350.2),
    ("Sandy Loam", 0.065, 0.41, 0.075, 1.89, 106.1),
    ("Loam", 0.078, 0.43, 0.036, 1.56, 24.96),
    ("Silt", 0.034, 0.46, 0.016, 1.37, 6.0),
    ("Silt Loam", 0.067, 0.45, 0.020, 1.41, 10.8),
    ("Sandy Clay Loam", 0.100, 0.39, 0.059, 1.48, 31.44),
    ("Clay Loam", 0.095, 0.41, 0.019, 1.31, 6.24),
    ("Silty Clay Loam", 0.089, 0.43, 0.010, 1.23, 1.68),
    ("Sandy Clay", 0.100, 0.38, 0.027, 1.23, 2.88),
    ("Silty Clay", 0.070, 0.36, 0.005, 1.09, 0.48),
    ("Clay", 0.068, 0.38, 0.008, 1.09, 4.8),
];

/// Factor converting a length given in cm to the project's length unit.
pub fn cm_to_unit(p: &Project) -> f64 {
    match p.units.length {
        LengthUnit::Mm => 10.0,
        LengthUnit::Cm => 1.0,
        LengthUnit::M => 0.01,
    }
}
/// Factor converting a time given in days to the project's time unit.
pub fn day_to_unit(p: &Project) -> f64 {
    match p.units.time {
        TimeUnit::Sec => 86400.0,
        TimeUnit::Min => 1440.0,
        TimeUnit::Hours => 24.0,
        TimeUnit::Days => 1.0,
        TimeUnit::Years => 1.0 / 365.0,
    }
}

pub fn apply_catalog(m: &mut SoilMaterial, idx: usize, prj_l: f64, prj_t: f64) {
    let c = SOIL_CATALOG[idx];
    m.name = c.0.to_string();
    m.qr = c.1;
    m.qs = c.2;
    m.alpha = c.3 / prj_l;
    m.n = c.4;
    m.ks = c.5 * prj_l * prj_t;
    m.l = 0.5;
    m.extra = [c.2, c.1, c.2, m.ks];
    m.qm = c.2;
    m.qs_w = c.2;
    m.alpha_w = 2.0 * m.alpha;
    m.ks_w = m.ks;
}

pub fn log_space(a: f64, b: f64, n: usize) -> Vec<f64> {
    let (la, lb) = (a.log10(), b.log10());
    (0..n).map(|i| 10f64.powf(la + (lb - la) * i as f64 / (n as f64 - 1.0))).collect()
}

/// Fixed points (depth, density_top, density_bottom) -> node depths (top to bottom).
/// Element length is proportional to the (linearly interpolated) density.
pub fn generate_nodes(depth: f64, n: usize, fixed: &[(f64, f64, f64)]) -> Vec<f64> {
    // density profile as piecewise-linear function of depth
    let mut pts: Vec<(f64, f64)> = vec![];
    let mut fp = fixed.to_vec();
    fp.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for (z, dt, db) in fp {
        pts.push((z.clamp(0.0, depth), dt));
        if (dt - db).abs() > 1e-12 {
            pts.push((z.clamp(0.0, depth) + 1e-9, db));
        }
    }
    if pts.is_empty() {
        pts.push((0.0, 1.0));
    }
    let dens = |z: f64| -> f64 {
        if z <= pts[0].0 {
            return pts[0].1;
        }
        for w in pts.windows(2) {
            if z <= w[1].0 {
                let f = (z - w[0].0) / (w[1].0 - w[0].0).max(1e-12);
                return w[0].1 + (w[1].1 - w[0].1) * f;
            }
        }
        pts[pts.len() - 1].1
    };
    let build = |c: f64| -> Vec<f64> {
        let mut z = vec![0.0];
        while *z.last().unwrap() < depth && z.len() < n + 2 {
            let zc = *z.last().unwrap();
            let h1 = c * dens(zc);
            let h2 = c * dens(zc + h1);
            z.push(zc + 0.5 * (h1 + h2));
        }
        z
    };
    let (mut lo, mut hi) = (1e-6 * depth, depth);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let z = build(mid);
        if z.len() - 1 > n - 1 || (z.len() - 1 == n - 1 && *z.last().unwrap() < depth) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let z = build(hi);
    let mut z: Vec<f64> = z.into_iter().take(n).collect();
    while z.len() < n {
        z.push(depth);
    }
    // rescale so that the last node is exactly at depth
    let last = *z.last().unwrap();
    if last > 0.0 {
        for v in z.iter_mut() {
            *v *= depth / last;
        }
    }
    z
}

fn fmt_num(v: f64, _r: std::ops::RangeInclusive<usize>) -> String {
    let a = v.abs();
    if a != 0.0 && (a >= 1e6 || a < 1e-4) {
        format!("{:e}", v)
    } else {
        let s = format!("{:.6}", v);
        let s = s.trim_end_matches('0').trim_end_matches('.').to_string();
        if s.is_empty() || s == "-" { "0".into() } else { s }
    }
}

fn parse_num(s: &str) -> Option<f64> {
    s.trim().replace(',', ".").replace(['d', 'D'], "e").parse::<f64>().ok()
}

pub fn dv(ui: &mut Ui, ch: &mut bool, v: &mut f64) {
    let speed = (v.abs() * 0.01).max(1e-6);
    if ui.add(DragValue::new(v).speed(speed).custom_formatter(fmt_num).custom_parser(parse_num)).changed() {
        *ch = true;
    }
}

pub fn dv_unit(ui: &mut Ui, ch: &mut bool, v: &mut f64, unit: &str) {
    let speed = (v.abs() * 0.01).max(1e-6);
    if ui.add(DragValue::new(v).speed(speed).custom_formatter(fmt_num).custom_parser(parse_num).suffix(format!(" {}", unit))).changed() {
        *ch = true;
    }
}

pub fn di(ui: &mut Ui, ch: &mut bool, v: &mut usize) {
    if ui.add(DragValue::new(v).speed(0.2)).changed() {
        *ch = true;
    }
}

pub fn field(ui: &mut Ui, ch: &mut bool, label: &str, v: &mut f64) {
    ui.label(label);
    dv(ui, ch, v);
    ui.end_row();
}

pub fn field_u(ui: &mut Ui, ch: &mut bool, label: &str, v: &mut f64, unit: &str) {
    ui.label(label);
    dv_unit(ui, ch, v, unit);
    ui.end_row();
}

pub fn check(ui: &mut Ui, ch: &mut bool, v: &mut bool, label: &str) {
    if ui.checkbox(v, label).changed() {
        *ch = true;
    }
}

pub fn section(ui: &mut Ui, title: &str) {
    ui.add_space(6.0);
    ui.label(egui::RichText::new(title).strong().size(15.0));
    ui.separator();
}
