//! File input/output: HYDRUS-1D legacy project import (Selector.in, Profile.dat,
//! Atmosph.in), native JSON projects and result export.

pub mod legacy;
pub mod native;
pub mod writer;

pub use legacy::read_legacy_project;
pub use native::{load_project, save_project};
