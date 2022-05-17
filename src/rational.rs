//! 
//! Defines the rational number class.
//! 
//! The struct `Rational` is used for avoiding the numerical errors.
//! 
//! The conversion from `T` to `Rational` is based on
//! the expansion of `T` by continued fractions.
//! 
//! This conversion assures the `1e-9` accuracy of the
//! rational representaion of the original real number.
//! 
use std::fmt;
use std::convert::From;
use std::cmp::Ordering;
use std::ops::{
    Add, Sub, Mul, Div, Neg,
    AddAssign, SubAssign, MulAssign, DivAssign
};

/// The maximum numerical error allowed in conversion
/// from a real number to a rational number.
const MAXIMUM_DENOMINATOR: i64 = 1_000_000_000_000;

/// Defines the rational numbers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Rational {
    numer: i64,
    denom: i64,
}


impl Rational {
    /// Maximum value that the `Rational` can deal with.
    #[inline(always)]
    pub const fn infinity() -> Rational {
        Rational { numer: i64::MAX, denom: 1 }
    }


    /// Minimum value that the `Rational` can deal with.
    #[inline(always)]
    pub const fn neg_infinity() -> Rational {
        Rational { numer: i64::MIN, denom: 1 }
    }


    /// Zero value for the `Rational`.
    #[inline(always)]
    pub const fn zero() -> Rational {
        Rational { numer: 0, denom: 1 }
    }


    /// Construct the `Rational` from raw representation
    #[inline(always)]
    pub fn raw(numer: i64, denom: i64) -> Self {

        if denom == 0 {
            panic!("Zero division occurred.");
        }

        let rational = Rational { numer, denom };
        rational.reduction()
    }

    /// Convert the rational number into an `f64`.
    #[inline(always)]
    pub fn as_f64(&self) -> f64 {
        self.numer as f64 / self.denom as f64
    }


    /// Reduce the fraction of `self`
    #[inline(always)]
    fn reduction(mut self) -> Self {
        if self.numer == 0 {
            Rational { numer: 0, denom: 1 }
        } else {
            let divisor = gcd(self.denom.abs(), self.numer.abs());
            self.numer /= divisor;
            self.denom /= divisor;

            self
        }
    }


    #[inline(always)]
    pub(crate) fn inv(mut self) -> Self {
        if self.numer == 0 {
            panic!("Zero division occurred.");
        }
        std::mem::swap(&mut self.numer, &mut self.denom);
        if self.denom < 0 {
            self.denom *= -1;
            self.numer *= -1;
        }
        self
    }


    #[inline(always)]
    pub(crate) fn is_zero(&self) -> bool {
        self.numer == 0
    }


    #[inline(always)]
    pub(crate) fn is_negative(&self) -> bool {
        self.numer < 0
    }


    #[inline(always)]
    pub(crate) fn is_positive(&self) -> bool {
        self.numer > 0
    }
}


/// Finds the greatest commondivisor by the Euclidean algorithm
fn gcd(m: i64, n: i64) -> i64 {
    if n == 0 {
        return m;
    }

    gcd(n, m % n)
}


impl<T> From<T> for Rational
    where T: Into<f64>
{
    /// Convert the real number to `Rational`.
    fn from(num: T) -> Self {
        let mut num = num.into();
        let mut sign: i64 = 1;

        // For simplicity, we set the `num` to be non-negative.
        if num < 0.0 {
            num  = num.abs();
            sign = -1;
        }

        // Reserve a 2x2 matrix of the form
        // 
        //  +-                       -+
        //  | p := p_n   q := p_{n-1} |
        //  | r := q_n   s := q_{n-1} |
        //  +-                       -+
        // 
        // where
        //  * p_n = a_n * p_{n-1} + p_{n-2}    if n >= 1
        //        = a_0                        if n == 0
        //        = 1                          otherwise
        //  * q_n = a_n * q_{n-1} + q_{n-2}    if n >= 1
        //        = 1                          if n == 0
        //        = 0                          otherwise
        let mut p: i64 = 1;
        let mut q: i64 = 0;
        let mut r: i64 = 0;
        let mut s: i64 = 1;

        let mut ai = num as i64;

        // LOOP until the denominator reaches some limit.
        // Note that this stopping criterion guarantees that
        // the resulting rational number is
        // `1.0 / MAXIMUM_DENOMINATOR.pow(2)`-close to the original num.
        // Check the number theory.
        while MAXIMUM_DENOMINATOR > r * ai + s {
            ai = num as i64;

            // Update the matrix
            let temp = p;
            p = p * ai + q;
            q = temp;

            let temp = r;
            r = r * ai + s;
            s = temp;

            if (num - ai as f64) < 1e-10 { break; }
            num = 1.0 / (num - ai as f64);

        }


        // Since `p` and `r` are coprime so that we do not need reduction
        Self {
            numer: sign * p,
            denom: r
        }
    }
}


impl Neg for Rational {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        Rational::raw(-1 * self.numer, self.denom)
    }
}


macro_rules! impl_rational_inner {
    ($t1:ty, $t2:ty) => (
        impl Add<$t1> for $t2 {
            type Output = Rational;
            #[inline]
            fn add(self, rhs: $t1) -> Self::Output {
                let lhs: Rational = self.into();
                let rhs: Rational = rhs.into();
                let left  = lhs.numer * rhs.denom;
                let right = rhs.numer * lhs.denom;

                let numer = left + right;
                let denom = lhs.denom * rhs.denom;

                Rational::raw(numer, denom)
            }
        }


        impl Sub<$t1> for $t2 {
            type Output = Rational;
            #[inline]
            fn sub(self, rhs: $t1) -> Self::Output {
                let lhs: Rational = self.into();
                let rhs: Rational = rhs.into();
                let rhs = -rhs;

                lhs + rhs
            }
        }


        impl Mul<$t1> for $t2 {
            type Output = Rational;
            #[inline]
            fn mul(self, rhs: $t1) -> Self::Output {
                let lhs: Rational = self.into();
                let rhs: Rational = rhs.into();
                let numer = lhs.numer * rhs.numer;
                let denom = lhs.denom * rhs.denom;

                Rational::raw(numer, denom)
            }
        }


        impl Div<$t1> for $t2 {
            type Output = Rational;
            #[inline]
            fn div(self, rhs: $t1) -> Self::Output {
                let lhs: Rational = self.into();
                let rhs: Rational = Rational::from(rhs).inv();


                lhs * rhs
            }
        }
    )
}


impl_rational_inner! { Rational, Rational }

macro_rules! impl_rational_bin_ops {
    ($($t:ty)*) => ($(
        impl_rational_inner! { $t, Rational }
        impl_rational_inner! { Rational, $t }
    )*)
}
impl_rational_bin_ops! { u8 u16 u32 i8 i16 i32 f32 f64 }


macro_rules! impl_rational_op_assign {
    ($($t:ty)*) => ($(
        impl AddAssign<$t> for Rational {
            #[inline]
            fn add_assign(&mut self, rhs: $t) {
                let rhs: Rational = rhs.into();
                let left  = self.numer * rhs.denom;
                let right = rhs.numer * self.denom;


                let numer = left + right;
                let denom = self.denom * rhs.denom;
                *self = Rational::raw(numer, denom);
            }
        }


        impl SubAssign<$t> for Rational {
            #[inline]
            fn sub_assign(&mut self, rhs: $t) {
                let rhs = -Rational::from(rhs);

                *self += rhs;
            }
        }


        impl MulAssign<$t> for Rational {
            #[inline]
            fn mul_assign(&mut self, rhs: $t) {
                let rhs: Rational = rhs.into();
                let numer = self.numer * rhs.numer;
                let denom = self.denom * rhs.denom;
                *self = Rational::raw(numer, denom);
            }
        }


        impl DivAssign<$t> for Rational {
            #[inline]
            fn div_assign(&mut self, rhs: $t) {
                let rhs: Rational = Rational::from(rhs).inv();
                *self *= rhs;
            }
        }
    )*)
}
impl_rational_op_assign! { u8 u16 u32 i8 i16 i32 f32 f64 Rational }


impl Default for Rational {
    fn default() -> Self {
        Rational { numer: 0, denom: 1 }
    }
}


impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({numer:>3}/{denom:>2})",
            numer = self.numer, denom = self.denom
        )
    }
}


impl PartialOrd for Rational {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let lhs = self.numer * other.denom;
        let rhs = other.numer * self.denom;

        lhs.partial_cmp(&rhs)
    }
}


