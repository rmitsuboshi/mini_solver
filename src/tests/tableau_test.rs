
use crate::tableau::Tableau;

#[test]
fn example_2() {
    let mut tableau = Tableau::new();
    tableau.set_objective(&[4.0, 1.0, 1.0]);

    tableau.extend(&[2.0, 1.0, 2.0], 4.0);
    tableau.extend(&[3.0, 3.0, 1.0], 3.0);


    println!("Init: \n{tableau}\n\n");


    tableau.solve().unwrap();
    println!("Solved: \n{tableau}\n");
    assert!(true);
}


// /// This test starts from the **Initial tableau** in page 52 of
// /// Linear and non-linear programming
// #[test]
// fn example_3() {
//     let mut tableau = Tableau::new();
//     tableau.set_objective(&[-2, 4, 7, 1, 5]);
// 
//     tableau.extend(&[-1, 1, 2, 1, 2], 7);
//     tableau.extend(&[-1, 2, 3, 1, 1], 6);
//     tableau.extend(&[-1, 1, 1, 2, 1], 4);
// 
// 
//     println!("Init: \n{tableau}\n\n");
// 
// 
//     tableau.solve().unwrap();
//     println!("Solved: \n{tableau}\n");
//     assert!(true);
// }
