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

/// Function state is set to 12 field elements or 96 bytes;
/// 4 elements of the state are reserved for capacity
pub const STATE_WIDTH: usize = 12;
/// 8 elements of the state are reserved for rate
pub const RATE_WIDTH: usize = 8;

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
        Fp::new(68622165223932372),
        Fp::new(6423413484163019822),
        Fp::new(11597469379118512612),
        Fp::new(6693206300350993384),
        Fp::new(2065922988851736597),
        Fp::new(13738169560427216462),
        Fp::new(3999864863072067615),
        Fp::new(6929688507379293701),
        Fp::new(10939726434240958765),
        Fp::new(17381238908510441624),
        Fp::new(17967444139409001592),
        Fp::new(4160929357686324685),
        Fp::new(962906224985451064),
        Fp::new(3412729377562427811),
        Fp::new(15381675735446757833),
        Fp::new(3355898117843322763),
        Fp::new(5653529456122806411),
        Fp::new(15408406709966119947),
        Fp::new(6213314615124871018),
        Fp::new(15018305290447027939),
        Fp::new(15681622795007326000),
        Fp::new(9161344146446309100),
        Fp::new(15017990541592107289),
        Fp::new(2694475991053984403),
        Fp::new(15358712594869477458),
        Fp::new(13550384111829059629),
        Fp::new(3013472385312798486),
        Fp::new(1783936659194570174),
        Fp::new(15240820427238996824),
        Fp::new(15376436872631782544),
        Fp::new(11478346532443214525),
        Fp::new(17932679705047642608),
        Fp::new(14385431538825757895),
        Fp::new(13907612839919404645),
        Fp::new(7553298447390248532),
        Fp::new(2498300816278705528),
        Fp::new(4176269975610862725),
        Fp::new(17313638398999194265),
        Fp::new(16431515285104634449),
        Fp::new(3159896988600578797),
        Fp::new(11564798120917478562),
        Fp::new(14777505433983289239),
        Fp::new(8895064120050592987),
        Fp::new(17728341733513937089),
        Fp::new(14125777345863740501),
        Fp::new(16699204621948978910),
        Fp::new(10083363504152939340),
        Fp::new(3277298255495454242),
        Fp::new(424451089698404918),
        Fp::new(13142429649367055605),
        Fp::new(7422006391945449305),
        Fp::new(13576825471923697616),
        Fp::new(11884450659505440692),
        Fp::new(16036524304026957320),
        Fp::new(12903042565238961434),
        Fp::new(1839784979566983995),
        Fp::new(10731956518745180741),
        Fp::new(9064634639324016657),
        Fp::new(104959021593271033),
        Fp::new(17894098489119110259),
        Fp::new(8125544303556132304),
        Fp::new(7461356050557721141),
        Fp::new(8010189432831732484),
        Fp::new(13051890357089566712),
        Fp::new(10974756077724511900),
        Fp::new(4849461727590433168),
        Fp::new(11544108603476792233),
        Fp::new(12277744068938583968),
        Fp::new(11774136866147683772),
        Fp::new(12413998154537224783),
        Fp::new(16679107944731793382),
        Fp::new(3158969135684587311),
        Fp::new(13486878740094718956),
        Fp::new(18029394588831285319),
        Fp::new(16226102210667248849),
        Fp::new(10118992141386478318),
        Fp::new(13042501354505552816),
        Fp::new(9237721818436805545),
        Fp::new(1623524271288047591),
        Fp::new(14589811346531505008),
        Fp::new(15617359384457502500),
        Fp::new(10534771303392773385),
        Fp::new(10730306957302970670),
        Fp::new(2750247152207505680),
        Fp::new(302827037583520315),
        Fp::new(6569066158181164366),
        Fp::new(17903309058277915166),
        Fp::new(11288934695231577930),
        Fp::new(13606628368490230374),
        Fp::new(14676429020030778210),
        Fp::new(7190670641540717266),
        Fp::new(8795656286887879175),
        Fp::new(17726648107018320332),
        Fp::new(2136872793132379667),
        Fp::new(9354103294418872262),
        Fp::new(16870240418942142605),
        Fp::new(17423606395803825692),
        Fp::new(17810152721973227038),
        Fp::new(5276616971899867946),
        Fp::new(14092158052521266715),
        Fp::new(11792370592662370154),
        Fp::new(10127728304017080154),
        Fp::new(13140262324230679697),
        Fp::new(13340973077862799362),
        Fp::new(4170094242677824358),
        Fp::new(2235912091038239392),
        Fp::new(13733348966500808022),
        Fp::new(17888903734393052461),
        Fp::new(15717305987983622902),
        Fp::new(5503791211243659890),
        Fp::new(12605988001826888844),
        Fp::new(10208287535816569072),
        Fp::new(9005567418461722042),
        Fp::new(11023164660361714158),
        Fp::new(6822632235366134315),
        Fp::new(7602167680112454408),
        Fp::new(7876854404428960942),
        Fp::new(16272556316998173434),
        Fp::new(1346376997412500931),
        Fp::new(14238348923043852337),
        Fp::new(10911507843421905466),
        Fp::new(3518798397517446944),
        Fp::new(10249993590248610388),
        Fp::new(16024363189778294874),
        Fp::new(18335874451370590186),
        Fp::new(1026053375930291079),
        Fp::new(11834921377926000913),
        Fp::new(6944670064989008921),
        Fp::new(6928775255898974859),
        Fp::new(3915568827579082059),
        Fp::new(7686494658467474318),
        Fp::new(1407506715057306219),
        Fp::new(8875238732439747549),
        Fp::new(18192514923978097968),
        Fp::new(9693453785298102915),
        Fp::new(9600704992970800406),
        Fp::new(11534525146993156014),
        Fp::new(13496523045778335117),
        Fp::new(8772713780970466358),
        Fp::new(7136352296718252935),
        Fp::new(11226407776063457820),
        Fp::new(3113793859449668510),
        Fp::new(11268892989914614068),
        Fp::new(16044230002153338636),
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
