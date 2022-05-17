#[cfg(test)]
mod rational_conversion {
    use crate::rational::Rational;
    #[test]
    fn pi() {
        let num = std::f64::consts::PI;
        let rational = Rational::from(num);

        // println!("{rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn e() {
        let num = std::f64::consts::E;
        let rational = Rational::from(num);

        // println!("{rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn quater() {
        let num = 1.0 / 4.0;
        let rational = Rational::from(num);

        // println!("{rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn three() {
        let num = 3.0;
        let rational = Rational::from(num);

        // println!("{rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn three_point_two() {
        let num = 3.2_f64;
        let rational = Rational::from(num);

        // println!("3.2 as rational: {rational:?}");
        assert!((num - rational.as_f64()).abs() < 1e-9);
    }


    #[test]
    fn i32_three() {
        let num = 3_i32;
        let rational = Rational::from(num);

        // println!("{rational:?}");
        assert!((num as f64 - rational.as_f64()).abs() < 1e-9);
    }
}


#[cfg(test)]
mod operators {
    use crate::rational::Rational;
    #[test]
    fn binary_operator() {
        let a = Rational::from(1.0 / 3.0);
        let b = Rational::from(1.0 / 2.0);
        let c = a + b;
        println!("a + b = {c}");
        assert_eq!(c, Rational::from(5.0 / 6.0));


        let a = Rational::from(1.0 / 3.0);
        let b = Rational::from(1.0 / 2.0);
        let c = a - b;
        println!("a - b = {c}");
        assert_eq!(c, Rational::from(-1.0 / 6.0));


        let a = Rational::from(1.0 / 3.0);
        let b = Rational::from(1.0 / 2.0);
        let c = a * b;
        println!("a * b = {c}");
        assert_eq!(c, Rational::from(1.0 / 6.0));


        let a = Rational::from(1.0 / 3.0);
        let b = Rational::from(1.0 / 2.0);
        let c = a / b;
        println!("a / b = {c}");
        assert_eq!(c, Rational::from(2.0 / 3.0));
    }




    #[test]
    fn assigns() {
        let mut a = Rational::from(1.0 / 3.0);
        let b = Rational::from(1.0 / 2.0);

        a += b;
        println!("a + b = {a}");
        assert_eq!(a, Rational::from(5.0 / 6.0));


        let mut a = Rational::from(1.0 / 3.0);
        let b = Rational::from(1.0 / 2.0);

        a -= b;
        println!("a - b = {a}");
        assert_eq!(a, Rational::from(-1.0 / 6.0));


        let mut a = Rational::from(1.0 / 3.0);
        let b = Rational::from(1.0 / 2.0);

        a *= b;
        println!("a * b = {a}");
        assert_eq!(a, Rational::from(1.0 / 6.0));


        let mut a = Rational::from(1.0 / 3.0);
        let b = Rational::from(1.0 / 2.0);

        a /= b;
        println!("a / b = {a}");
        assert_eq!(a, Rational::from(2.0 / 3.0));
    }
}


#[cfg(test)]
mod others {
    use crate::rational::Rational;
    #[test]
    fn reduction() {
        let r1 = Rational::raw(2, 4);
        let r2 = Rational::from(0.5);

        assert!(r1 == r2);
    }
}


