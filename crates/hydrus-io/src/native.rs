use hydrus_core::{HydrusError, Project, Results};
use std::path::Path;

pub fn save_project(p: &Project, path: &Path) -> Result<(), HydrusError> {
    let s = serde_json::to_string_pretty(p).map_err(|e| HydrusError::Io(e.to_string()))?;
    std::fs::write(path, s).map_err(|e| HydrusError::Io(format!("{}: {}", path.display(), e)))
}

pub fn load_project(path: &Path) -> Result<Project, HydrusError> {
    let s = std::fs::read_to_string(path).map_err(|e| HydrusError::Io(format!("{}: {}", path.display(), e)))?;
    serde_json::from_str(&s).map_err(|e| HydrusError::Parse(e.to_string()))
}

pub fn save_results(r: &Results, path: &Path) -> Result<(), HydrusError> {
    let s = serde_json::to_string(r).map_err(|e| HydrusError::Io(e.to_string()))?;
    std::fs::write(path, s).map_err(|e| HydrusError::Io(format!("{}: {}", path.display(), e)))
}

pub fn load_results(path: &Path) -> Result<Results, HydrusError> {
    let s = std::fs::read_to_string(path).map_err(|e| HydrusError::Io(format!("{}: {}", path.display(), e)))?;
    serde_json::from_str(&s).map_err(|e| HydrusError::Parse(e.to_string()))
}
