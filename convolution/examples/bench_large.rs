use convolution::{convolution_large, intt, ntt, P1};
use std::time::Instant;
fn main() {
    let n = 1usize << 24;
    let a: Vec<u32> = (0..n as u32).map(|i| i.wrapping_mul(2654435761) % P1).collect();
    let t = Instant::now();
    let c = convolution_large(&a, &a);
    println!("convolution_large 2^24 x 2^24: {:?} (len {})", t.elapsed(), c.len());
    let mut f = vec![1u32; 1 << 23];
    let t = Instant::now();
    for _ in 0..4 {
        ntt::<P1>(&mut f);
        intt::<P1>(&mut f);
    }
    println!("8 transforms of 2^23: {:?}", t.elapsed());
}
