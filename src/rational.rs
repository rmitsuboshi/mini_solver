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
use std::ops::{Neg, Add, Sub, Mul, Div};

/// The maximum numerical error allowed in conversion
/// from a real number to a rational number.
const MAXIMUM_DENOMINATOR: u64 = 1_000_000_000;

/// Defines the rational numbers.
#[derive(Debug, PartialEq, Eq)]
pub struct Rational {
    sign:  i8,
    numer: u64,
    denom: u64,
}


impl Rational {
    /// Construct the `Rational` from row representation
    #[inline(always)]
    fn row(sign: i8, numer: u64, denom: u64) -> Self {
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

            if num == ai as f64 {
                break;
            }
            num = 1.0 / (num - ai as f64);
        }


        // Since `p` and `r` are coprime so that we do not need reduction
        Rational::row(sign, p, r)
    }
}


impl Neg for Rational {
    type Output = Self;
    fn neg(self) -> Self::Output {
        let sign = self.sign * -1;
        Rational::row(sign, self.numer, self.denom)
    }
}


impl Add for Rational {
    type Output = Rational;
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

        Self::row(sign, numer, denom)
    }
}


impl Sub for Rational {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        let mut rhs = rhs;
        rhs.sign *= -1;

        self + rhs
    }
}


impl Mul for Rational {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let sign  = self.sign  * rhs.sign;
        let numer = self.numer * rhs.numer;
        let denom = self.denom * rhs.denom;

        Rational::row(sign, numer, denom)
    }
}


impl Div for Rational {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let sign  = self.sign  * rhs.sign;
        let numer = self.numer * rhs.denom;
        let denom = self.denom * rhs.numer;

        Rational::row(sign, numer, denom)
    }
}


#[cfg(test)]
mod rational_conversion {
    use super::*;
    #[test]
    fn pi() {
        let num = std::f64::consts::PI;
        let rational = Rational::from(num);

        println!("{rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn e() {
        let num = std::f64::consts::E;
        let rational = Rational::from(num);

        println!("{rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn quater() {
        let num = 1.0 / 4.0;
        let rational = Rational::from(num);

        println!("{rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn three() {
        let num = 3.0;
        let rational = Rational::from(num);

        println!("{rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn i32_three() {
        let num = 3_i32;
        let rational = Rational::from(num);

        println!("{rational:?}");
        assert!((num as f64 - rational.as_f64()).abs() < 1e-9);
    }
}


#[cfg(test)]
mod others {
    use super::*;
    #[test]
    fn reduction() {
        let r1 = Rational::row(1, 2, 4);
        let r2 = Rational::from(0.5);

        assert!(r1 == r2);
    }


    #[test]
    fn add() {
        let r1 = Rational::row( 1, 2, 9);
        let r2 = Rational::row(-1, 3, 4);

        let r3 = r1 + r2;

        let expected = Rational::row(-1, 19, 36);

        assert_eq!(r3, expected);

        let r1 = Rational::from(1e9);
        let r2 = Rational::from(1e9);

        let r3 = r1 + r2;
        println!("Added: {r3:?}");
        assert!(true);
    }


    #[test]
    fn sub() {
        let r1 = Rational::row( 1, 2, 9);
        let r2 = Rational::row(-1, 3, 4);

        let r3 = r1 - r2;

        let expected = Rational::row(1, 35, 36);

        assert_eq!(r3, expected);

        let r1 = Rational::from(1e9);
        let r2 = Rational::from(1e9);

        let r3 = r1 - r2;
        println!("Subtracted: {r3:?}");
        assert!(true);
    }


    #[test]
    fn mul() {
        let r1 = Rational::row( 1, 2, 9);
        let r2 = Rational::row(-1, 3, 4);

        let r3 = r1 * r2;

        let expected = Rational::row(-1, 1, 6);

        assert_eq!(r3, expected);

        let r1 = Rational::from(1e9);
        let r2 = Rational::from(1e9);

        let r3 = r1 * r2;
        println!("Multiplied: {r3:?}");
        assert!(true);
    }


    #[test]
    fn div() {
        let r1 = Rational::row( 1, 2, 9);
        let r2 = Rational::row(-1, 4, 3);

        let r3 = r1 / r2;

        let expected = Rational::row(-1, 1, 6);

        assert_eq!(r3, expected);

        let r1 = Rational::from(1e9);
        let r2 = Rational::from(1e9);

        let r3 = r1 / r2;
        println!("Divided: {r3:?}");
        assert!(true);
    }
}

