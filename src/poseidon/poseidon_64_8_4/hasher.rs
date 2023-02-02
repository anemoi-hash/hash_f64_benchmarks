// Copyright (c) 2021-2023 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! Hasher trait implementation for Poseidon

use core::convert::TryInto;

use super::digest::PoseidonDigest;
use super::{apply_permutation, DIGEST_SIZE, RATE_WIDTH, STATE_WIDTH};
use crate::traits::Hasher;

use cheetah::Fp;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
/// A Poseidon Hash over Fp
pub struct PoseidonHash {
    state: [Fp; STATE_WIDTH],
    idx: usize,
}

impl Default for PoseidonHash {
    fn default() -> Self {
        Self {
            state: [Fp::zero(); STATE_WIDTH],
            idx: 0,
        }
    }
}

impl Hasher<Fp> for PoseidonHash {
    type Digest = PoseidonDigest;

    fn hash(bytes: &[Fp]) -> Self::Digest {
        // initialize state to all zeros
        let mut state = [Fp::zero(); STATE_WIDTH];

        let mut i = 0;
        for &element in bytes.iter() {
            state[i] += element;
            i += 1;
            if i % RATE_WIDTH == 0 {
                apply_permutation(&mut state);
                i = 0;
            }
        }

        if i > 0 {
            state[i] += Fp::one();
            i += 1;

            while i % RATE_WIDTH != 0 {
                state[i] = Fp::zero();
                i += 1;
            }

            apply_permutation(&mut state);
        }

        PoseidonDigest::new(state[..DIGEST_SIZE].try_into().unwrap())
    }

    fn merge(values: &[Self::Digest; 2]) -> Self::Digest {
        let mut state = [Fp::zero(); STATE_WIDTH];
        let digest1 = values[0].as_elements();
        let digest2 = values[1].as_elements();
        // Uses Jive compression to fill the whole state and perform a single permutation
        state[..RATE_WIDTH].copy_from_slice(digest1);
        state[RATE_WIDTH..STATE_WIDTH].copy_from_slice(digest2);
        apply_permutation(&mut state);

        let mut result = [Fp::zero(); DIGEST_SIZE];
        for (i, r) in result.iter_mut().enumerate() {
            *r = digest1[i] + digest2[i] + state[i] + state[i + STATE_WIDTH / 2];
        }

        PoseidonDigest::new(state[..DIGEST_SIZE].try_into().unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poseidon_hash() {
        // Hardcoded input / output list generated from the
        // Sagemath code at <https://github.com/Nashtare/goldilocks-hadeshash>

        let input_data = [
            [Fp::zero(); 4],
            [Fp::one(); 4],
            [
                Fp::new(6569531088581465466),
                Fp::new(16633559500314554860),
                Fp::new(5886252279469044386),
                Fp::new(7254493620959209513),
            ],
            [
                Fp::new(14451671431242580637),
                Fp::new(12408955884259542783),
                Fp::new(6699700836437690609),
                Fp::new(6930173452022628101),
            ],
            [
                Fp::new(18395072982230283623),
                Fp::new(17422138460915468648),
                Fp::new(302928547891428217),
                Fp::new(17610521183302992333),
            ],
            [
                Fp::new(18068273339104138164),
                Fp::new(1911447746536554313),
                Fp::new(15648741352653160565),
                Fp::new(6668708735458596065),
            ],
            [
                Fp::new(11143425136806953711),
                Fp::new(18102105924703184544),
                Fp::new(7182071222255043774),
                Fp::new(9055541327943176217),
            ],
            [
                Fp::new(8655442662468800475),
                Fp::new(16505130820128335063),
                Fp::new(6493352664926984222),
                Fp::new(14967895957577870330),
            ],
            [
                Fp::new(9323549209372186683),
                Fp::new(10793117907827765253),
                Fp::new(182614859344069638),
                Fp::new(5967861484968599356),
            ],
            [
                Fp::new(10467738873745304557),
                Fp::new(16593721817482730487),
                Fp::new(16833487465053662813),
                Fp::new(13645546544866037446),
            ],
            [
                Fp::new(5197638704396892579),
                Fp::new(2260999998602521401),
                Fp::new(14339582556914007502),
                Fp::new(2059336184840445928),
            ],
            [
                Fp::new(11474010417503429257),
                Fp::new(6562199940984629687),
                Fp::new(1595612647018183356),
                Fp::new(15125081316530224427),
            ],
        ];

        // Generated from <https://github.com/Nashtare/goldilocks-hadeshash>
        let output_data = [
            [
                Fp::new(11238356114731002609),
                Fp::new(1817990373503212709),
                Fp::new(16350570428038487437),
                Fp::new(6614850657780244420),
            ],
            [
                Fp::new(8041252807103335744),
                Fp::new(13257044945362841715),
                Fp::new(3670917474116293738),
                Fp::new(3404292384679666951),
            ],
            [
                Fp::new(13340802203704231334),
                Fp::new(934216092306344555),
                Fp::new(15476427833015873916),
                Fp::new(6797038095273097506),
            ],
            [
                Fp::new(6177086991679129976),
                Fp::new(3412748413220833157),
                Fp::new(5084951472248880373),
                Fp::new(11112442887702404566),
            ],
            [
                Fp::new(5854725137514659260),
                Fp::new(12129411317354654554),
                Fp::new(15660366587783925572),
                Fp::new(8235503794135806089),
            ],
            [
                Fp::new(13719489232811922375),
                Fp::new(7466384294629751394),
                Fp::new(6421785214497179280),
                Fp::new(5117420494878168921),
            ],
            [
                Fp::new(16286727075544267236),
                Fp::new(14947718250167055938),
                Fp::new(5945660771842955379),
                Fp::new(18278342428984400181),
            ],
            [
                Fp::new(15194862054182941774),
                Fp::new(5997581910190898133),
                Fp::new(13948082932178022451),
                Fp::new(1038658250246226378),
            ],
            [
                Fp::new(4066880073755180414),
                Fp::new(16966267686805338947),
                Fp::new(8487193859494998664),
                Fp::new(17381647074250870936),
            ],
            [
                Fp::new(15065645856954930412),
                Fp::new(2598005939680799392),
                Fp::new(15280880696196736397),
                Fp::new(13052727819461015270),
            ],
            [
                Fp::new(271779881367265232),
                Fp::new(13017842100463108287),
                Fp::new(9657458714975527415),
                Fp::new(17681722652525508162),
            ],
            [
                Fp::new(13507364143313959582),
                Fp::new(17781762901407462186),
                Fp::new(4757605172902052948),
                Fp::new(7599183616003930405),
            ],
        ];

        for (input, expected) in input_data.iter().zip(output_data) {
            assert_eq!(expected, PoseidonHash::hash(input).to_elements());
        }
    }
}
