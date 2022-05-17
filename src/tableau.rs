//! Defines a model of a linear programming.
//! 
use crate::rational::Rational;
use crate::error::LpError;


use std::fmt;
use std::collections::{HashMap, HashSet};
use std::ops::RangeBounds;


/// Simplex tableau.
/// The first lines of `matrix` and `rhs` represents
/// the relative cost coefficients.
#[derive(Debug)]
pub(crate) struct Tableau {
    matrix: Vec<Vec<Rational>>,
    rhs:    Vec<Rational>,
    basis:  HashMap<usize, usize>
}


impl Tableau {
    /// Construct an empty tableau.
    /// The first row of the tableau is 
    ///  reserved for an objective function.
    #[inline]
    pub(super) fn new() -> Self {
        let obj = Vec::new();
        Self {
            matrix: vec![obj],
            rhs:    vec![Rational::from(0.0)],
            basis:  HashMap::new(),
        }
    }


    #[inline]
    pub(super) fn relative_costs(&self) -> &Vec<Rational> {
        &self.matrix[0]
    }


    #[inline]
    pub(super) fn extend<T>(&mut self, coef: &[T], rhs: T)
        where Rational: From<T>,
              T: Clone
    {
        let coef = coef.into_iter()
            .map(|c| Rational::from(c.clone()))
            .collect::<Vec<_>>();
        self.matrix.push(coef);
        self.rhs.push(Rational::from(rhs));
    }


    #[inline]
    pub(super) fn set_objective<T>(&mut self, coef: &[T])
        where Rational: From<T>,
              T: Clone
    {
        let coef = coef.into_iter()
            .map(|c| Rational::from(c.clone()))
            .collect::<Vec<_>>();
        self.matrix[0] = coef;
    }


    /// Returns true if current solution is optimal.
    #[inline]
    pub(crate) fn is_optimal(&self) -> bool {
        self.matrix[0]
            .iter()
            .all(|c| !c.is_negative())
    }


    /// Solve the linear program.
    #[inline]
    pub(crate) fn solve(&mut self) -> Result<(), LpError> {
        // Initialize the tableau in the canonical form.
        self.to_canonical()?;


        Ok(())
    }


    /// Convert the tableau into the canonical form.
    #[inline]
    fn to_canonical(&mut self) -> Result<(), LpError> {
        let m = self.matrix.len();
        let n = self.matrix.iter()
            .map(|row| row.len())
            .max()
            .unwrap();


        // DEBUG
        assert!(self.matrix.iter().all(|row| row.len() == n));
        assert_eq!(self.matrix.len(), self.rhs.len());


        // let dummy_vars = (0..m).into_iter()
        //     .map(|i| {
        //         let mut v = vec![Rational::from(0.0); m];
        //         v[i] = Rational::from(1.0);
        //         v
        //     })
        //     .collect::<Vec<_>>();

        let mut rcc = vec![Rational::from(0.0); m+n];

        for row in self.matrix.iter().skip(1) {
            for (r, &v) in rcc.iter_mut().zip(row.iter()) {
                *r -= v;
            }
        }


        // TODO
        // The following code assumes `self.matrix` is full rank,
        // which implies the non-degeneracy of the solution.
        // I need to review the code to adapt the degenerate case.
        // loop {
        for _ in 0..m {
            let basis = rcc.iter()
                .enumerate()
                .filter_map(|(i, r)| {
                    if r.is_negative() {
                        Some((i, r))
                    } else {
                        None
                    }
                })
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap());


            let basis = match basis {
                None => {
                    break;
                },
                Some((i, _)) => i
            };


            // Find povit
            let pivot = self.matrix.iter()
                .zip(&self.rhs)
                .enumerate()
                .skip(1)
                .filter_map(|(i, (m, &r))| {
                    if m[basis] > Rational::zero() {
                        Some((i, r / m[basis]))
                    } else {
                        None
                    }
                })
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap());


            let pivot = match pivot {
                None => {
                    // TODO
                    return Err(LpError::Unbounded);
                },
                Some((i, _)) => i
            };


            let divisor = self.matrix[pivot][basis];
            assert!(!divisor.is_zero());
            for m in self.matrix[pivot].iter_mut() {
                *m /= divisor;
            }
            self.rhs[pivot] /= divisor;


            let pivot_row = std::mem::replace(
                &mut self.matrix[pivot],
                Vec::with_capacity(0)
            );
            let pivot_rhs = self.rhs[pivot];


            // Update the tableau
            let iter = self.matrix.iter_mut()
                .zip(self.rhs.iter_mut())
                .enumerate()
                .skip(1);
            for (i, (row, rhs)) in iter {
                if i == pivot { continue; }

                let multiplier = row[basis];

                for (val, p) in row.iter_mut().zip(pivot_row.iter()) {
                    *val -= multiplier * *p;
                }

                *rhs -= multiplier * pivot_rhs;
            }


            self.matrix[pivot] = pivot_row;


            // Update relative cost coefficients
            let multiplier = rcc[basis];
            for (r, p) in rcc.iter_mut().zip(self.matrix[pivot].iter()) {
                *r -= multiplier * *p;
            }

            let column = self.basis.entry(pivot).or_insert(0);
            *column = basis;
        }


        let mut obj_coef = std::mem::replace(
            &mut self.matrix[0],
            Vec::with_capacity(0)
        );


        let mut obj_val = Rational::zero();
        for (&row, &column) in self.basis.iter() {
            let multiplier = obj_coef[column];

            let iter = obj_coef.iter_mut()
                .zip(self.matrix[row].iter());
            for (c, r) in iter {
                *c -= multiplier * *r;
            }


            obj_val -= multiplier * self.rhs[row];
        }

        self.matrix[0] = obj_coef;
        self.rhs[0] = obj_val;


        Ok(())
    }
}


impl fmt::Display for Tableau {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let formatted = self.matrix[1..]
            .iter()
            .zip(&self.rhs[1..])
            .map(|(coef, rhs)| {
                let lhs = coef.iter()
                    .map(|c| format!("{c}"))
                    .collect::<Vec<_>>()
                    .join("\t");
                format!("{lhs}\t{rhs}")
            })
            .collect::<Vec<_>>()
            .join("\n");

        write!(f, "{formatted}\n")?;


        write!(f, "{lhs}\t{rhs}",
            lhs = self.matrix[0].iter()
                .map(|coef| format!("{coef}"))
                .collect::<Vec<_>>()
                .join("\t"),
            rhs = self.rhs[0]
        )
    }
}
