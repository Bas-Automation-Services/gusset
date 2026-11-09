//! Sheet metal part model, flat layout and fold map.

pub mod check;
pub mod clash;
mod cone;
mod cut;
pub mod flat;
pub mod fold;
mod foldline;
mod hem;
mod profile;
pub mod section;
pub mod spec;

/// Why a part can't be built.
#[derive(Clone, Debug, PartialEq)]
pub enum SheetError {
    /// The spec contradicts itself or describes something unbuildable.
    Invalid(String),
    /// Valid sheet metal that Gusset doesn't build yet.
    Unsupported(String),
    /// The kernel refused the model (a bug in Gusset or the kernel).
    Kernel(String),
}

impl std::fmt::Display for SheetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(m) => write!(f, "invalid: {m}"),
            Self::Unsupported(m) => write!(f, "unsupported: {m}"),
            Self::Kernel(m) => write!(f, "kernel: {m}"),
        }
    }
}

impl std::error::Error for SheetError {}
