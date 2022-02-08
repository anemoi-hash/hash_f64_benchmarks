// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use criterion::{black_box, criterion_group, criterion_main, Criterion};

extern crate hash;
use cheetah::Fp;
use hash::rescue_64_8_4::{digest::RescueDigest, hasher::RescueHash};
use hash::traits::Hasher;
use rand::rngs::OsRng;
use rand::thread_rng;
use rand::RngCore;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("rescue-64-8-4 - merge", |bench| {
        let v: [RescueDigest; 2] = [RescueHash::hash(&[1u8]), RescueHash::hash(&[2u8])];
        bench.iter(|| RescueHash::merge(black_box(&v)))
    });

    c.bench_function("rescue-64-8-4 - hash 25 Fp elements", |bench| {
        let mut v = [Fp::zero(); 25];
        for e in v.iter_mut() {
            *e = Fp::random(OsRng);
        }
        bench.iter(|| RescueHash::digest(black_box(&v)))
    });

    c.bench_function("rescue-64-8-4 - hash 10KB", |bench| {
        let mut data = vec![0u8; 10 * 1024];
        let mut rng = thread_rng();
        rng.fill_bytes(&mut data);
        bench.iter(|| RescueHash::hash(black_box(&data)))
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = criterion_benchmark);
criterion_main!(benches);
