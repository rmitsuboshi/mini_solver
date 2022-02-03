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
use std::convert::From;
use std::ops::{Add, Sub, Mul, Div, Neg};

/// The maximum numerical error allowed in conversion
/// from a real number to a rational number.
const MAXIMUM_DENOMINATOR: u64 = 1_000_000_000_000;

/// Defines the rational numbers.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Rational {
    sign:  i8,
    numer: u64,
    denom: u64,
}


impl Rational {
    /// Maximum value that the `Rational` can deal with.
    #[inline(always)]
    pub const fn infinity() -> Rational {
        Rational { sign: 1, numer: u64::MAX, denom: 1 }
    }

    /// Minimum value that the `Rational` can deal with.
    #[inline(always)]
    pub const fn neg_infinity() -> Rational {
        Rational { sign: -1, numer: u64::MAX, denom: 1 }
    }

    /// Construct the `Rational` from raw representation
    #[inline(always)]
    pub fn raw(sign: i8, numer: u64, denom: u64) -> Self {
        if sign == 0 {
            panic!("sign must be a non-zero integer");
        }

        let sign = sign.signum();

        if denom == 0 {
            panic!("Zero division occurred!");
        }
        let mut r = Rational { sign, numer, denom };
        r.reduction();
        r
    }

    /// Convert the rational number into an `f64`.
    pub fn as_f64(&self) -> f64 {
        self.sign as f64 * (self.numer as f64 / self.denom as f64)
    }


    /// Reduce the fraction of `self`
    #[inline(always)]
    fn reduction(&mut self) {
        let divisor = gcd(self.denom, self.numer);
        self.numer /= divisor;
        self.denom /= divisor;
    }


    #[inline(always)]
    pub(crate) fn inv(mut self) -> Self {
        std::mem::swap(&mut self.numer, &mut self.denom);
        self
    }
}


/// Finds the greatest commondivisor by the Euclidean algorithm
fn gcd(m: u64, n: u64) -> u64 {
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
        let mut sign: i8 = 1;

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
        let mut p: u64 = 1;
        let mut q: u64 = 0;
        let mut r: u64 = 0;
        let mut s: u64 = 1;

        let mut ai = num as u64;

        // LOOP until the denominator reaches some limit.
        // Note that this stopping criterion guarantees that
        // the resulting rational number is
        // `1.0 / MAXIMUM_DENOMINATOR.pow(2)`-close to the original num.
        // Check the number theory.
        while MAXIMUM_DENOMINATOR > r * ai + s {
            ai = num as u64;

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
        Rational::raw(sign, p, r)
    }
}


impl Neg for Rational {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self::Output {
        let sign = self.sign * -1;
        Rational::raw(sign, self.numer, self.denom)
    }
}


impl Add for Rational {
    type Output = Rational;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {

        let left  = self.sign as i64 * self.numer as i64 *  rhs.denom as i64;
        let right =  rhs.sign as i64 *  rhs.numer as i64 * self.denom as i64;

        let mut numer = left + right;
        let mut sign  = 1;
        if numer < 0 {
            sign  = -1;
            numer = numer.abs();
        }
        let numer = numer as u64;
        let denom = self.denom * rhs.denom;

        Self::raw(sign, numer, denom)
    }
}


impl Sub for Rational {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        let mut rhs = rhs;
        rhs.sign *= -1;

        self + rhs
    }
}


impl Mul for Rational {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        let sign  = self.sign  * rhs.sign;
        let numer = self.numer * rhs.numer;
        let denom = self.denom * rhs.denom;

        Rational::raw(sign, numer, denom)
    }
}


impl Div for Rational {
    type Output = Self;
    #[inline]
    fn div(self, rhs: Self) -> Self::Output {
        let sign  = self.sign  * rhs.sign;
        let numer = self.numer * rhs.denom;
        let denom = self.denom * rhs.numer;

        Rational::raw(sign, numer, denom)
    }
}


impl Default for Rational {
    fn default() -> Self {
        Rational::raw(1, 0, 1)
    }
}
