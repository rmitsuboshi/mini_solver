#[cfg(test)]
mod var_test {
    use crate::var::Var;
    use crate::rational::Rational;
    #[test]
    fn test_01() {
        let name = "x";

        let x = Var::new(&name, 0..);

        let expected = Var {
            name:  String::from("x"),
            range: Rational::from(0)..Rational::infinity()
        };

        assert_eq!(x, expected);
    }
}

