//! # hydrus-core
//! A Rust port of the HYDRUS-1D numerical engine for one-dimensional variably
//! saturated water flow, heat and solute transport.

pub mod atm;
pub mod error;
pub mod examples;
pub mod heat;
pub mod lenhard;
pub mod material;
pub mod meteo;
pub mod model;
pub mod output;
pub mod report;
pub mod root;
pub mod sim;
pub mod solute;
pub mod water;
pub mod vapor;

pub use error::HydrusError;
pub use model::*;
pub use output::Results;
pub use sim::{Simulation, StepStatus};
