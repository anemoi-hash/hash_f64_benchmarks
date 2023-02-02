// Copyright (c) 2021-2023 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use cheetah::Fp;

/// Digest for Poseidon
mod digest;
/// Hasher for Poseidon
mod hasher;
/// MDS matrix for Poseidon
mod mds;
/// Round constants for Poseidon
mod round_constants;

pub use digest::PoseidonDigest;
pub use hasher::PoseidonHash;

// RESCUE CONSTANTS
// ================================================================================================

/// Function state is set to 8 field elements or 64 bytes;
/// 4 elements of the state are reserved for capacity
pub const STATE_WIDTH: usize = 8;
/// 8 elements of the state are reserved for rate
pub const RATE_WIDTH: usize = 4;

/// Seven elements (32-bytes) are returned as digest.
pub const DIGEST_SIZE: usize = 4;

/// The number of full rounds is set to 8 to provide 128-bit security level with security margin.
/// Computed using <https://github.com/Nashtare/goldilocks-hadeshash>. Note that considering
/// BBLP22 attack may reduce the actual security level.
pub const NUM_HALF_FULL_ROUNDS: usize = 4;

/// The number of partial rounds is set to 22 to provide 128-bit security level with security margin.
/// Computed using <https://github.com/Nashtare/goldilocks-hadeshash>. Note that considering
/// BBLP22 attack may reduce the actual security level.
pub const NUM_PARTIAL_ROUNDS: usize = 22;

// HELPER FUNCTIONS
// ================================================================================================

#[inline(always)]
/// Applies matrix-vector multiplication of the current
/// hash state with the Poseidon MDS matrix.
pub(crate) fn apply_mds(state: &mut [Fp; STATE_WIDTH]) {
    let mut result = [Fp::zero(); STATE_WIDTH];
    for (i, r) in result.iter_mut().enumerate() {
        for (j, s) in state.iter().enumerate() {
            *r += mds::MDS[i * STATE_WIDTH + j] * s;
        }
    }

    state.copy_from_slice(&result);
}

#[inline(always)]
pub(crate) fn apply_full_sbox(state: &mut [Fp; STATE_WIDTH]) {
    state.iter_mut().for_each(pow_7);
}

#[inline(always)]
/// Applies exponentiation of the current element by 7
pub(crate) fn pow_7(x: &mut Fp) {
    let t2 = x.square();
    let t4 = t2.square();
    *x *= t2 * t4;
}

// RESCUE PERMUTATION
// ================================================================================================

/// Applies Poseidon permutation to the provided state.
pub(crate) fn apply_permutation(state: &mut [Fp; STATE_WIDTH]) {
    for i in 0..NUM_HALF_FULL_ROUNDS {
        apply_full_round(state, i);
    }

    for i in 0..NUM_PARTIAL_ROUNDS {
        apply_partial_round(state, NUM_HALF_FULL_ROUNDS + i);
    }

    for i in 0..NUM_HALF_FULL_ROUNDS {
        apply_full_round(state, NUM_PARTIAL_ROUNDS + NUM_HALF_FULL_ROUNDS + i);
    }
}

/// Poseidon full round function;
#[inline(always)]
pub(crate) fn apply_full_round(state: &mut [Fp; STATE_WIDTH], step: usize) {
    // determine which round constants to use
    let ark = round_constants::ARK[step % (2 * NUM_HALF_FULL_ROUNDS + NUM_PARTIAL_ROUNDS)];

    for i in 0..STATE_WIDTH {
        state[i] += ark[i];
    }

    apply_full_sbox(state);
    apply_mds(state);
}

/// Poseidon full round function;
#[inline(always)]
pub(crate) fn apply_partial_round(state: &mut [Fp; STATE_WIDTH], step: usize) {
    // determine which round constants to use
    let ark = round_constants::ARK[step % (2 * NUM_HALF_FULL_ROUNDS + NUM_PARTIAL_ROUNDS)];

    for i in 0..STATE_WIDTH {
        state[i] += ark[i];
    }

    pow_7(&mut state[0]);
    apply_mds(state);
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::OsRng;

    const INV_MDS: [Fp; STATE_WIDTH * STATE_WIDTH] = [
        Fp::new(5232231990142370708),
        Fp::new(1456292618197229064),
        Fp::new(2585965691860279897),
        Fp::new(8596043838921749788),
        Fp::new(8501843150900594164),
        Fp::new(15661370495349411547),
        Fp::new(8335511691205600286),
        Fp::new(9472710036860008327),
        Fp::new(12909463115742621643),
        Fp::new(17812790881402916949),
        Fp::new(17084854063114100360),
        Fp::new(1743778736584506472),
        Fp::new(931032031734907883),
        Fp::new(1182029477279544872),
        Fp::new(4218744783891988473),
        Fp::new(16674624160942802487),
        Fp::new(2105794376044984331),
        Fp::new(11936208987889491064),
        Fp::new(3020825386030820825),
        Fp::new(16385199278367979228),
        Fp::new(15663647911850050672),
        Fp::new(13385591541441954316),
        Fp::new(10851630570597723463),
        Fp::new(4009661047516767761),
        Fp::new(3257879938515494937),
        Fp::new(17198486073375638065),
        Fp::new(5354045220144466484),
        Fp::new(1350784205834799709),
        Fp::new(70213105078461578),
        Fp::new(1303555312236393700),
        Fp::new(10775713498649348807),
        Fp::new(11941690491548710456),
        Fp::new(1362147347033862645),
        Fp::new(4465586076580454760),
        Fp::new(3496493523135751019),
        Fp::new(8152485191176261434),
        Fp::new(14326166154068252015),
        Fp::new(11352384956803611766),
        Fp::new(16036981034433906839),
        Fp::new(6095410813907601113),
        Fp::new(13638149741021140024),
        Fp::new(3727429434539678642),
        Fp::new(1214738526414198494),
        Fp::new(919217129542191474),
        Fp::new(10705274983492634622),
        Fp::new(11625725003743044981),
        Fp::new(9985537962900532294),
        Fp::new(3026827904242613038),
        Fp::new(15810453901481577735),
        Fp::new(6918287518874692323),
        Fp::new(4549388663112521530),
        Fp::new(13958999623219919370),
        Fp::new(12799213451866585641),
        Fp::new(7736627847927338554),
        Fp::new(6988424089742464285),
        Fp::new(16371861096259113405),
        Fp::new(15530649383689463387),
        Fp::new(12515141549291147991),
        Fp::new(1079698422632889951),
        Fp::new(17605257525871053374),
        Fp::new(2371697140819154646),
        Fp::new(16208462304098589476),
        Fp::new(15354444556517360040),
        Fp::new(11228143819241617025),
    ];

    /// Applies matrix-vector multiplication of the current
    /// hash state with the inverse Poseidon MDS matrix.
    fn apply_inv_mds(state: &mut [Fp; STATE_WIDTH]) {
        let mut result = [Fp::zero(); STATE_WIDTH];
        for (i, r) in result.iter_mut().enumerate() {
            for (j, s) in state.iter().enumerate() {
                *r += INV_MDS[i * STATE_WIDTH + j] * s;
            }
        }

        state.copy_from_slice(&result);
    }

    #[test]
    fn test_mds() {
        let mut state = [Fp::zero(); STATE_WIDTH];
        let mut rng = OsRng;

        for _ in 0..100 {
            for s in state.iter_mut() {
                *s = Fp::random(&mut rng);
            }

            let state_copy = state;
            apply_mds(&mut state);

            // Check that matrix multiplication was consistent
            apply_inv_mds(&mut state);
            assert_eq!(state, state_copy);
        }
    }
}
