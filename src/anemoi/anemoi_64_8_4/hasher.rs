// Copyright (c) 2021-2022 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! Hasher trait implementation for Anemoi

use core::convert::TryInto;

use super::digest::AnemoiDigest;
use super::{apply_permutation, DIGEST_SIZE, NUM_COLUMNS, RATE_WIDTH, STATE_WIDTH};
use crate::traits::Hasher;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

use cheetah::Fp;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
/// A Anemoi Hash over Fp
pub struct AnemoiHash {
    state: [Fp; STATE_WIDTH],
    idx: usize,
}

impl Default for AnemoiHash {
    fn default() -> Self {
        Self {
            state: [Fp::zero(); STATE_WIDTH],
            idx: 0,
        }
    }
}

impl Hasher<Fp> for AnemoiHash {
    type Digest = AnemoiDigest;

    fn hash(bytes: &[Fp]) -> Self::Digest {
        // initialize state to all zeros.
        let mut state = [Fp::zero(); STATE_WIDTH];

        let sigma = if bytes.len() % RATE_WIDTH == 0 {
            Fp::one()
        } else {
            Fp::zero()
        };

        let mut i = 0;
        for &element in bytes.iter() {
            state[i] += element;
            i += 1;
            if i % RATE_WIDTH == 0 {
                apply_permutation(&mut state);
                i = 0;
            }
        }

        // If the message length is not a multiple of RATE_WIDTH, we append 1 to the rate cell
        // next to the one where we previously appended the last message element. This is
        // guaranted to be in the rate registers (i.e. to not require an extra permutation before
        // adding this constant) if sigma is equal to zero.
        if sigma.is_zero().into() {
            state[i] += Fp::one();
            apply_permutation(&mut state);
        }

        // We then add sigma to the last capacity register of the capacity.
        state[STATE_WIDTH - 1] += sigma;

        AnemoiDigest::new(state[..DIGEST_SIZE].try_into().unwrap())
    }

    // This merge function uses the compression approach of Anemoi-Jive
    // to save one Anemoi permutation call , which would be necessary if
    // using the regular Anemoi-Sponge to absorb two digests, both of
    // size RATE_WIDTH.
    fn merge(values: &[Self::Digest; 2]) -> Self::Digest {
        let mut state = [Fp::zero(); STATE_WIDTH];
        let digest1 = values[0].as_elements();
        let digest2 = values[1].as_elements();
        state[..RATE_WIDTH].copy_from_slice(digest1);
        state[RATE_WIDTH..STATE_WIDTH].copy_from_slice(digest2);
        apply_permutation(&mut state);

        let mut result = [Fp::zero(); DIGEST_SIZE];
        for (i, r) in result.iter_mut().enumerate() {
            *r = digest1[i] + digest2[i] + state[i] + state[i + NUM_COLUMNS];
        }

        AnemoiDigest::new(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_anemoi_hash() {
        // Hardcoded input / output list generated from the
        // Sagemath code at https://github.com/anemoi/anemoi-hash/

        let input_data = [
            vec![Fp::zero(); 8],
            vec![Fp::one(); 8],
            vec![
                Fp::zero(),
                Fp::zero(),
                Fp::zero(),
                Fp::zero(),
                Fp::one(),
                Fp::one(),
                Fp::one(),
                Fp::one(),
            ],
            vec![
                Fp::one(),
                Fp::one(),
                Fp::one(),
                Fp::one(),
                Fp::zero(),
                Fp::zero(),
                Fp::zero(),
                Fp::zero(),
            ],
            vec![Fp::new(17034304680378990593)],
            vec![Fp::new(14717522243836148645), Fp::new(13083744320516212099)],
            vec![
                Fp::new(12636374470493764585),
                Fp::new(4633685103624134705),
                Fp::new(3387113220078373436),
            ],
            vec![
                Fp::new(8456646190515084682),
                Fp::new(13422660296625058046),
                Fp::new(6412264356237719015),
                Fp::new(15681324668660747245),
            ],
            vec![
                Fp::new(5024098825481674854),
                Fp::new(17074804226935063103),
                Fp::new(12706716815026291874),
                Fp::new(15111404609115389377),
                Fp::new(14660803971611056841),
            ],
            vec![
                Fp::new(12593193490016384806),
                Fp::new(18163695353480745896),
                Fp::new(12379842404544411425),
                Fp::new(15341182951106044393),
                Fp::new(1499922380158571885),
                Fp::new(4387632416457430195),
            ],
        ];

        let output_data = [
            [
                Fp::new(14545405583772908021),
                Fp::new(1213427980949533238),
                Fp::new(5108363964065086531),
                Fp::new(11576263589351466107),
            ],
            [
                Fp::new(9034041492535792788),
                Fp::new(1880373787890293705),
                Fp::new(16508665783675182648),
                Fp::new(15004799909697058781),
            ],
            [
                Fp::new(4158861503480679584),
                Fp::new(15525323270749352759),
                Fp::new(3322238964359227464),
                Fp::new(429391678459631464),
            ],
            [
                Fp::new(11670455842144808209),
                Fp::new(6155050407516467312),
                Fp::new(16834017798801702367),
                Fp::new(5352738999477411269),
            ],
            [
                Fp::new(2120754002113216041),
                Fp::new(4176373977615675990),
                Fp::new(5380916603234155760),
                Fp::new(8459213842511311406),
            ],
            [
                Fp::new(3544328427077449214),
                Fp::new(15836917233876247314),
                Fp::new(7333626051256456686),
                Fp::new(17611887051951947693),
            ],
            [
                Fp::new(18402038807138490723),
                Fp::new(10713759395658533517),
                Fp::new(15564080858015170609),
                Fp::new(8529908362848541681),
            ],
            [
                Fp::new(17302804114778865555),
                Fp::new(849947101445175894),
                Fp::new(17782524703869970480),
                Fp::new(13986605488954541106),
            ],
            [
                Fp::new(10995727710267236590),
                Fp::new(622816764469086576),
                Fp::new(6998511466743063111),
                Fp::new(10059021300276080542),
            ],
            [
                Fp::new(3566384133654955307),
                Fp::new(2093744914411683144),
                Fp::new(12462545152074913887),
                Fp::new(4526022743638213689),
            ],
        ];

        for (input, expected) in input_data.iter().zip(output_data) {
            assert_eq!(expected, AnemoiHash::hash(input).to_elements());
        }
    }
}
