use convolution::{intt, ntt, P1};
use std::time::Instant;
fn main() {
    for h in [20usize, 23] {
        let n = 1usize << h;
        let mut a: Vec<u32> = (0..n as u32).map(|i| i.wrapping_mul(2654435761) % P1).collect();
        let t = Instant::now();
        for _ in 0..3 {
            ntt::<P1>(&mut a);
            intt::<P1>(&mut a);
        }
        println!("2^{h}: {:?} per transform", t.elapsed() / 6);
    }
}
