//! Simple example demonstrating B-field computation from a cuboid magnet

use magrslib::sources::magnets::Cuboid;
use magrslib::types::Observers;
use magrslib::get_b;

fn main() {
    // Create a 1cm cube magnet at the origin with 1 Tesla polarization in z-direction
    let magnet = Cuboid::builder()
        .dimension([0.01, 0.01, 0.01])      // 1cm x 1cm x 1cm
        .polarization([0.0, 0.0, 1.0])       // 1 Tesla in +z direction
        .position([0.0, 0.0, 0.0])           // Centered at origin
        .build()
        .unwrap();

    // Define observer positions (in meters)
    let observers = Observers::from_vec3s(vec![
        [0.0, 0.0, 0.02],   // 2cm above magnet
        [0.02, 0.0, 0.0],   // 2cm to the right
        [0.0, 0.02, 0.0],   // 2cm in front
    ]);

    // Compute B-field
    let b_field = get_b(&[magnet.into()], &observers);

    // Display results
    println!("B-field computation for 1cm cuboid magnet with 1T z-polarization\n");
    println!("Observer positions and computed B-fields (in Tesla):");
    println!("{:-<70}", "");
    println!("{:>20} {:>20} {:>20}", "Position (m)", "B-field (T)", "Magnitude (T)");
    println!("{:-<70}", "");

    for i in 0..observers.len() {
        let obs_pos = observers.as_array();
        let bx = b_field[[i, 0]];
        let by = b_field[[i, 1]];
        let bz = b_field[[i, 2]];
        let b_mag = (bx * bx + by * by + bz * bz).sqrt();

        println!(
            "[{:6.3}, {:6.3}, {:6.3}]  [{:8.5}, {:8.5}, {:8.5}]  {:8.5}",
            obs_pos[[i, 0]], obs_pos[[i, 1]], obs_pos[[i, 2]],
            bx, by, bz,
            b_mag
        );
    }
    println!("{:-<70}", "");
}
