//! LP status error
use std::fmt;


/// Defines the Errors occur in linear program.
#[derive(Debug, PartialEq, Eq)]
pub enum LpError {
    /// The feasible region is unbounded
    Unbounded,
    /// There is no solution that satisfies the constraints
    Infeasible,
}


impl fmt::Display for LpError {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            LpError::Unbounded
                => write!(f, "Unbounded"),
            LpError::Infeasible
                => write!(f, "Infeasible"),
        }
    }
}
