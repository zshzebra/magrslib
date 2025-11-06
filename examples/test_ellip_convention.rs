// Test to understand ellip crate's parameter convention

fn main() {
    println!("=== Testing ellip parameter convention ===\n");

    // Scipy uses parameter m = k²
    // K(m=0.25) in scipy = K(k²=0.25) = K(k=0.5) = 1.8540746773013719

    println!("--- ELLIPK Tests ---");
    println!("ellip::ellipk(0.25) = {:?}", ellip::ellipk(0.25).unwrap());
    println!("ellip::ellipk(0.5) = {:?}", ellip::ellipk(0.5).unwrap());
    println!("ellip::ellipk(sqrt(0.25)) = {:?}", ellip::ellipk((0.25_f64).sqrt()).unwrap());
    println!("Scipy K(m=0.25) = 1.8540746773013719");
    println!("→ ellipk uses modulus k, need to call with sqrt(m)!\n");

    println!("--- ELLIPE Tests ---");
    println!("ellip::ellipe(0.25) = {:?}", ellip::ellipe(0.25).unwrap());
    println!("Scipy E(m=0.25) = 1.4674622093394271");
    println!("→ ellipe already uses parameter m! ✓\n");

    println!("--- ELLIPPI Tests ---");
    println!("ellip::ellippi(0.2, 0.25) = {:?}", ellip::ellippi(0.2, 0.25).unwrap());
    println!("ellip::ellippi(0.2, sqrt(0.25)) = {:?}", ellip::ellippi(0.2, (0.25_f64).sqrt()).unwrap());
    println!("Scipy Π(n=0.2, m=0.25) = 1.9867935682987379");

    // Try different n values too
    println!("\nTrying different n parameter conversions:");
    let n: f64 = 0.2;
    let m: f64 = 0.25;
    println!("ellippi(n={}, m={}) variations:", n, m);
    println!("  ellippi({}, {}) = {:?}", n, m, ellip::ellippi(n, m).unwrap());
    println!("  ellippi({}, sqrt({})) = {:?}", n, m, ellip::ellippi(n, m.sqrt()).unwrap());

    // Maybe n also needs conversion?
    println!("  ellippi(sqrt({}), {}) = {:?}", n, m, ellip::ellippi(n.sqrt(), m).unwrap());
    println!("  ellippi(sqrt({}), sqrt({})) = {:?}", n, m, ellip::ellippi(n.sqrt(), m.sqrt()).unwrap());

    println!("\n--- Identity Check: Π(0, m) = K(m) ---");
    let pi_0_m = ellip::ellippi(0.0, m).unwrap();
    let k_k = ellip::ellipk(m).unwrap();
    let k_m = ellip::ellipk(m.sqrt()).unwrap();

    println!("Π(0, {}) = {}", m, pi_0_m);
    println!("K({}) = {} (if ellipk uses modulus k)", m, k_k);
    println!("K(sqrt({})) = {} (if ellipk uses param m)", m, k_m);

    if (pi_0_m - k_m).abs() < 1e-10 {
        println!("→ Π(0, m) = K(sqrt(m)), so ellippi uses parameter m, ellipk uses modulus k (INCONSISTENT)");
    } else if (pi_0_m - k_k).abs() < 1e-10 {
        println!("→ Π(0, m) = K(m), both use the SAME convention");
    } else {
        println!("→ No identity match - different convention entirely!");
    }

    println!("\n--- Testing n parameter transformation ---");
    let n_scipy: f64 = 0.2;
    let m_scipy: f64 = 0.25;
    let k = m_scipy.sqrt();  // k = sqrt(m)

    println!("Scipy: Π(n={}, m={}) = 1.9867935682987379", n_scipy, m_scipy);
    println!("\nTrying different n transformations with k=sqrt(m):");
    println!("  ellippi({}, {}) = {:?}", n_scipy, k, ellip::ellippi(n_scipy, k).unwrap());
    println!("  ellippi(1-{}, {}) = {:?}", n_scipy, k, ellip::ellippi(1.0 - n_scipy, k).unwrap());
    println!("  ellippi(-{}, {}) = {:?}", n_scipy, k, ellip::ellippi(-n_scipy, k).unwrap());
    println!("  ellippi({}/({}-1), {}) = {:?}", n_scipy, n_scipy, k, ellip::ellippi(n_scipy/(n_scipy-1.0), k).unwrap());
}
