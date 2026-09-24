//! File input/output: HYDRUS-1D legacy project import (Selector.in, Profile.dat,
//! Atmosph.in), native JSON projects and result export.

pub mod legacy;
pub mod moist_dep;
pub mod native;
pub mod writer;

pub use legacy::{
    read_legacy_project, read_nod_inf, sanitize_fortran_line, LegacyProfileNode,
    LegacyProfileTimeBlock,
};
pub use moist_dep::{read_moist_dep, write_moist_dep};
pub use native::{load_project, save_project};
