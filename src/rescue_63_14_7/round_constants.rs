// Copyright (c) 2021 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::{NUM_HASH_ROUNDS, STATE_WIDTH};
use cheetah::Fp;

/// Additive Round Keys constants for Rescue,
/// computed using algorithm 5 from <https://eprint.iacr.org/2020/1143.pdf>
pub const ARK: [[Fp; STATE_WIDTH * 2]; NUM_HASH_ROUNDS] = [
    [
        Fp::new(0x1f0a24c4ed42b7df),
        Fp::new(0x23966eb7b343720e),
        Fp::new(0x14bbfa44ff5b743f),
        Fp::new(0xe664c9986cb8a9e),
        Fp::new(0x4119c0c05c7ecd7e),
        Fp::new(0x32ce8901c4293486),
        Fp::new(0x3e68c5c98d4b4cb8),
        Fp::new(0x2a63cb703b3572a0),
        Fp::new(0x370da12ca562d56d),
        Fp::new(0x1da6d3d90c15b05d),
        Fp::new(0x2eaf791c2a38d572),
        Fp::new(0x2bb3461b78a1f224),
        Fp::new(0x397fea4351111fe6),
        Fp::new(0x1fe11370e8a410d8),
        Fp::new(0x287cc57b73b216c4),
        Fp::new(0x31d43141acff6960),
        Fp::new(0x24a060674a8713ea),
        Fp::new(0x41181e510c8dbc78),
        Fp::new(0x28eea3b98c6b9ee7),
        Fp::new(0x3ce13e44655b3186),
        Fp::new(0xd825b0db466b46d),
        Fp::new(0x4d55c6b88df6972),
        Fp::new(0x11847585b3e06d1e),
        Fp::new(0x2686f84c862f4896),
        Fp::new(0x3faec01f47b5a468),
        Fp::new(0x32010b89ce5a5c16),
        Fp::new(0x3a8a353735812e88),
        Fp::new(0x19acb2c8c419d69),
    ],
    [
        Fp::new(0x2a67f2def9434b18),
        Fp::new(0x5938a7b8a911856),
        Fp::new(0x345ae70ea4fdf960),
        Fp::new(0x383b69f66d65b559),
        Fp::new(0x1eea20fb14fcc9cf),
        Fp::new(0x40d6bd565c2cf37),
        Fp::new(0x368944bc3e1ae57d),
        Fp::new(0x3449b6bb664d184a),
        Fp::new(0x416c90ce460c7258),
        Fp::new(0x1e270c06d813795d),
        Fp::new(0xf9e18710b4874f2),
        Fp::new(0x2c13bdf5b184c1b2),
        Fp::new(0xce723b3c4e32ff6),
        Fp::new(0x3b4f9580c0a03588),
        Fp::new(0x309fa4dadff08e09),
        Fp::new(0xdb312001f5d8e61),
        Fp::new(0x3553e9ffa77bc9ee),
        Fp::new(0x177ddfc84dcab572),
        Fp::new(0x3a2e9ce68b1a5115),
        Fp::new(0x7858b9979f77e46),
        Fp::new(0x29f254b2d69334fc),
        Fp::new(0x2104a7ca8fb2d70f),
        Fp::new(0x394054b0791650e9),
        Fp::new(0x246e4e5b18f07e54),
        Fp::new(0x6fedf25cfedde0c),
        Fp::new(0x31c2caf0082ccb62),
        Fp::new(0x236548e19637b41a),
        Fp::new(0x29f6a61610faf9c1),
    ],
    [
        Fp::new(0x819707ec9a67813),
        Fp::new(0x20c2f6a293cc0a87),
        Fp::new(0x1d3709b68192a421),
        Fp::new(0x645fe7901df574e),
        Fp::new(0x21f889c67d7e3ba3),
        Fp::new(0x2c460161b3914236),
        Fp::new(0x1bac0ef8e49616b0),
        Fp::new(0x70aff1238d34e11),
        Fp::new(0x38c168bc2f68832b),
        Fp::new(0x412a168e21bf2b53),
        Fp::new(0x287ba21a54154ce6),
        Fp::new(0x21ff3f2653cdd1eb),
        Fp::new(0x2f173cca8668ff8f),
        Fp::new(0x3a8696d71835f516),
        Fp::new(0x6d60271f19bdec5),
        Fp::new(0xb3f07f7039df8bf),
        Fp::new(0x345c46a0cc5fb5ed),
        Fp::new(0x1f284a385da803d2),
        Fp::new(0x31272f6ad3863843),
        Fp::new(0x1ce856afd3537362),
        Fp::new(0x1b008de8c1c3ca3a),
        Fp::new(0x3acbde69cfc423a3),
        Fp::new(0x1fb8d8f1e44dfd37),
        Fp::new(0x25bcbadef12e3474),
        Fp::new(0x3ce171702963c13d),
        Fp::new(0x276fd7aed3f312a2),
        Fp::new(0x40a4e27e824c2d),
        Fp::new(0x3e9c674a8002ca62),
    ],
    [
        Fp::new(0x1239d2a50c98ee11),
        Fp::new(0x128aa086b005e82c),
        Fp::new(0x2a7e981d4efe8b00),
        Fp::new(0xaaec7ccfe7f2324),
        Fp::new(0x15e97e0d0e1b7358),
        Fp::new(0x447697cb53e2335),
        Fp::new(0x353490eeef4f707c),
        Fp::new(0xa844d78c57c82ae),
        Fp::new(0x37209c0bb193f4a),
        Fp::new(0x39ed12078e2206da),
        Fp::new(0x3f1f3091852b09cc),
        Fp::new(0x208c0a8b88fc9e3e),
        Fp::new(0x1444d6073161e6e3),
        Fp::new(0x1393abe3ac44a731),
        Fp::new(0x954901d34f08c2f),
        Fp::new(0x434d68beff8bc3c),
        Fp::new(0x289b878613113d7b),
        Fp::new(0x11571f4113f74aea),
        Fp::new(0x295d7a74aecbd738),
        Fp::new(0x3fc9cb8bc9e5ce6b),
        Fp::new(0xdbd33109a6a49f7),
        Fp::new(0x1322f7c31be4be9f),
        Fp::new(0x1ce0bb10c065e5d3),
        Fp::new(0xb952ab8628cb682),
        Fp::new(0x40f814133438cbdf),
        Fp::new(0x25722ec1766cd448),
        Fp::new(0x5d49fd46561472d),
        Fp::new(0x991bb35cb7052ca),
    ],
    [
        Fp::new(0x2823a1b0d2646a2c),
        Fp::new(0x3a0dc712d799b107),
        Fp::new(0x8c6e77050f662e4),
        Fp::new(0x2fbcf9c0fe368312),
        Fp::new(0x399d5795595c979b),
        Fp::new(0x3ddbaac6e5cab794),
        Fp::new(0x3e3abc3c104634a5),
        Fp::new(0x58218618c424b24),
        Fp::new(0x28d45bdc3c867372),
        Fp::new(0x1f04d7b485f02826),
        Fp::new(0x12b38b2b8757364d),
        Fp::new(0x8faf4eef692d005),
        Fp::new(0x1d9175f53e6c64a1),
        Fp::new(0x30cd988a0ce61ca3),
        Fp::new(0x1dd0bdfacfc9ff80),
        Fp::new(0x22245428977637c7),
        Fp::new(0x1ce88dba021e6543),
        Fp::new(0x40293474d9e4eb72),
        Fp::new(0xf6618f49fa7a229),
        Fp::new(0x12dad1e5ae9d67),
        Fp::new(0x36923eef059e5918),
        Fp::new(0xa5f8accc1bd2c6e),
        Fp::new(0xf7cb307d0f31bd5),
        Fp::new(0x26522ba1cc28828c),
        Fp::new(0x1090b6e701b628b6),
        Fp::new(0x3a97813f9a5a82eb),
        Fp::new(0x8fe5b3cb78acf),
        Fp::new(0x17d261078e8b32c3),
    ],
    [
        Fp::new(0x1e4acd6ff3382eef),
        Fp::new(0x3d17ca86a7651d49),
        Fp::new(0x2d804138338b7f72),
        Fp::new(0x152788e7fc018214),
        Fp::new(0x22bbf35179db337),
        Fp::new(0xeaae2acc8190a60),
        Fp::new(0x20196ecc727e035b),
        Fp::new(0x2698f7a3485ea605),
        Fp::new(0x19c8deaacaf65443),
        Fp::new(0xde9eb1d8981506a),
        Fp::new(0x3935397a2890ca7a),
        Fp::new(0xedf58cc48004974),
        Fp::new(0x136d5c0e55f1170b),
        Fp::new(0x2b2aa453d1bdb322),
        Fp::new(0x219c52e273d5a977),
        Fp::new(0x16dfe0dc6eb7456f),
        Fp::new(0x188f3d51ce9efd25),
        Fp::new(0x4605ce20f8e3da4),
        Fp::new(0x380547e70a777777),
        Fp::new(0x2b5c71584d414c99),
        Fp::new(0x37507fe409d8339d),
        Fp::new(0x2d40ddcb229e22ac),
        Fp::new(0x11ca5ec22bea5bf4),
        Fp::new(0xd7af7899ff5344b),
        Fp::new(0xb5a2470994c53da),
        Fp::new(0xa3c737f8f73b866),
        Fp::new(0x67bb69f329faed4),
        Fp::new(0xbace861a72434a),
    ],
    [
        Fp::new(0x3f000ab3f6cd3732),
        Fp::new(0xa30620299cb250),
        Fp::new(0xfe651ea8b4878a7),
        Fp::new(0x2c65270c6cc5551d),
        Fp::new(0x40a7e4c790eefa0c),
        Fp::new(0x1ed6c721331db81b),
        Fp::new(0x268cc27cc0f0dc74),
        Fp::new(0x34684cc51ff6b99f),
        Fp::new(0x24ded79a44672f43),
        Fp::new(0x1cdccb24cf696ed5),
        Fp::new(0x11ea5613d92dade1),
        Fp::new(0x1ab67d02d67b0c37),
        Fp::new(0xe09d70f3d80af1c),
        Fp::new(0x333b98b0e36e8bd8),
        Fp::new(0x3b4a7dec3394c686),
        Fp::new(0x19378054e5a32c32),
        Fp::new(0x40bc3d05a339f10),
        Fp::new(0x1f1d63e484ef0021),
        Fp::new(0x15ce5213567419f),
        Fp::new(0x1b08d33c33000502),
        Fp::new(0x3d6176e2dbb0bb17),
        Fp::new(0x1cdde9da8a6083c),
        Fp::new(0x2f552238068a6303),
        Fp::new(0x16794069461cc8dc),
        Fp::new(0x3ac5564dcac156b3),
        Fp::new(0x124db611cbae1828),
        Fp::new(0x128ca6911c1f9f69),
        Fp::new(0xba762a29f69c886),
    ],
];
