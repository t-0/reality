//! Model amendments: declare new field content on top of the Standard Model,
//! check it for consistency, and derive predictions without fitting.

pub mod field;
pub mod gates;
pub mod predict;
pub mod proton;
pub mod running;
pub mod sm;
pub mod spec;
pub mod two_loop;

pub use predict::{evaluate, Prediction, Report};
pub use spec::Amendment;
