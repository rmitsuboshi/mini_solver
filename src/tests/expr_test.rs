#[cfg(test)]
mod multiply {
    use crate::rational::Rational;
    use crate::expr::Expr;
    use crate::var::Var;
    #[test]
    fn test_01() {
        let x = Var::new(&"x", 0.0_f64..);

        let expected = {
            let mut linear = std::collections::HashMap::new();
            linear.insert(&x, 3.2_f64.into());

            Expr {
                linear, constant: Rational::default()
            }
        };

        let expr = 3.2_f64 * &x;

        // println!("{expr:?}");


        assert_eq!(expr, expected);
    }
}
