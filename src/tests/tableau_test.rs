
use crate::tableau::Tableau;

#[test]
fn tableau_test() {
    let mut tableau = Tableau::new();
    tableau.set_objective(&[4.0, 1.0, 1.0]);

    tableau.extend(&[2.0, 1.0, 2.0], 4.0);
    tableau.extend(&[3.0, 3.0, 1.0], 3.0);


    println!("Init: \n{tableau}\n\n");


    tableau.solve();
    println!("Canonical: \n{tableau}\n");
    assert!(true);
}
