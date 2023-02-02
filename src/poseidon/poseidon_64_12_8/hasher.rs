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
        state[..DIGEST_SIZE].copy_from_slice(values[0].as_elements());
        state[DIGEST_SIZE..RATE_WIDTH].copy_from_slice(values[1].as_elements());
        apply_permutation(&mut state);

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
            [Fp::zero(); 8],
            [Fp::one(); 8],
            [
                Fp::new(1854186296009606004),
                Fp::new(11095952686320463486),
                Fp::new(1070369769687682777),
                Fp::new(5343455209103938928),
                Fp::new(968038157567727892),
                Fp::new(1911442725744079089),
                Fp::new(3511989225473916843),
                Fp::new(6616889790702343087),
            ],
            [
                Fp::new(3323504858425716648),
                Fp::new(7109147163868326989),
                Fp::new(14497056920727745515),
                Fp::new(9735581944069477148),
                Fp::new(5059436216248265674),
                Fp::new(8687938807992652586),
                Fp::new(1488381457275351091),
                Fp::new(8510592317513165148),
            ],
            [
                Fp::new(12647595006684713791),
                Fp::new(16038075251281990240),
                Fp::new(7438545350394189372),
                Fp::new(15802953074880975807),
                Fp::new(6346189954103175228),
                Fp::new(5773357976661773899),
                Fp::new(5484471395266653033),
                Fp::new(14137053600772666816),
            ],
            [
                Fp::new(18363177278156677388),
                Fp::new(10049421326749834952),
                Fp::new(16125014385895703365),
                Fp::new(6554523500998570287),
                Fp::new(2895147016209460605),
                Fp::new(11144136859546019030),
                Fp::new(15119850878804967273),
                Fp::new(11158823799820779550),
            ],
            [
                Fp::new(2518793593577017764),
                Fp::new(6216798898306003941),
                Fp::new(2825204806749464984),
                Fp::new(9726080653174839564),
                Fp::new(2242725448663943001),
                Fp::new(3577540141786750843),
                Fp::new(3981544363677331272),
                Fp::new(3755058646578539437),
            ],
            [
                Fp::new(748671726625133085),
                Fp::new(16457973663592160784),
                Fp::new(10335281425178509852),
                Fp::new(3513269257673292609),
                Fp::new(17127892871817407645),
                Fp::new(11069817215846519480),
                Fp::new(4229244207371962035),
                Fp::new(13521307887360970542),
            ],
            [
                Fp::new(6661956584137326678),
                Fp::new(14782478975879682753),
                Fp::new(6944879797312322601),
                Fp::new(4579128815758068637),
                Fp::new(17186672198037898028),
                Fp::new(3213310086807780296),
                Fp::new(2624424946607387331),
                Fp::new(9976624920510702649),
            ],
            [
                Fp::new(12837232005968094574),
                Fp::new(17893295118307922862),
                Fp::new(5534459983427118936),
                Fp::new(13298264085937057286),
                Fp::new(14747519306540139279),
                Fp::new(3311609252193461521),
                Fp::new(16260381367032123171),
                Fp::new(7244632472809117001),
            ],
            [
                Fp::new(10121226257530989283),
                Fp::new(17367044537479253),
                Fp::new(12749234296868924688),
                Fp::new(9137015959255606132),
                Fp::new(835775949203352656),
                Fp::new(14610618494036701430),
                Fp::new(5221240419910088067),
                Fp::new(8924356331170882496),
            ],
            [
                Fp::new(13848057127539975177),
                Fp::new(13398126226246728166),
                Fp::new(5903652180275854393),
                Fp::new(6108248817856190712),
                Fp::new(9675303959472926462),
                Fp::new(14630442237764846362),
                Fp::new(969289219878513206),
                Fp::new(2432080838343066343),
            ],
        ];

        // Generated from <https://github.com/Nashtare/goldilocks-hadeshash>
        let output_data = [
            [
                Fp::new(15200053835380386475),
                Fp::new(7591143832479391787),
                Fp::new(13656430951822477850),
                Fp::new(15075858837303362416),
            ],
            [
                Fp::new(3215110189550297925),
                Fp::new(7516979907590577539),
                Fp::new(17567539633201555290),
                Fp::new(4107813960770512499),
            ],
            [
                Fp::new(15914729853967880741),
                Fp::new(16572743799122799482),
                Fp::new(5697070153484269464),
                Fp::new(4664357757829673898),
            ],
            [
                Fp::new(9735476015564042468),
                Fp::new(4122294128506713140),
                Fp::new(12898425918614055059),
                Fp::new(10823450933787025625),
            ],
            [
                Fp::new(7415215638649750072),
                Fp::new(18420005941885010640),
                Fp::new(13647481739944295852),
                Fp::new(475855118340298999),
            ],
            [
                Fp::new(10504139805571550621),
                Fp::new(12247695678219710578),
                Fp::new(11620044851717693425),
                Fp::new(15799850271774583902),
            ],
            [
                Fp::new(7882923087848433944),
                Fp::new(14080275475924169620),
                Fp::new(1287233461951528305),
                Fp::new(13766298192192916055),
            ],
            [
                Fp::new(11403604029004570029),
                Fp::new(17692898850746996029),
                Fp::new(17489435212744446315),
                Fp::new(14010279029037125313),
            ],
            [
                Fp::new(17644851234394584663),
                Fp::new(11169548518106916847),
                Fp::new(17932632661388521827),
                Fp::new(13652373624339841158),
            ],
            [
                Fp::new(10867544499815531068),
                Fp::new(10478642228242492296),
                Fp::new(3724676220942414780),
                Fp::new(12405640553565815424),
            ],
            [
                Fp::new(9845876827594337898),
                Fp::new(8082649445755307351),
                Fp::new(15835539655927413717),
                Fp::new(9288782363442385984),
            ],
            [
                Fp::new(6802076088657719444),
                Fp::new(13373393578746327142),
                Fp::new(14703479682310980618),
                Fp::new(4292034528789329986),
            ],
        ];

        for (input, expected) in input_data.iter().zip(output_data) {
            assert_eq!(expected, PoseidonHash::hash(input).to_elements());
        }
    }
}
