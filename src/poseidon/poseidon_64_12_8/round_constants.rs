// Copyright (c) 2021-2023 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::{NUM_HALF_FULL_ROUNDS, NUM_PARTIAL_ROUNDS, STATE_WIDTH};
use cheetah::Fp;

/// Additive Round Keys constants for Poseidon.
pub(crate) const ARK: [[Fp; STATE_WIDTH]; 2 * NUM_HALF_FULL_ROUNDS + NUM_PARTIAL_ROUNDS] = [
    [
        Fp::new(3342624911463171251),
        Fp::new(6781356195391537436),
        Fp::new(4697929572322733707),
        Fp::new(4179687232228901671),
        Fp::new(17841073646522133059),
        Fp::new(18340176721233187897),
        Fp::new(13152929999122219197),
        Fp::new(6306257051437840427),
        Fp::new(4974451914008050921),
        Fp::new(11258703678970285201),
        Fp::new(581736081259960204),
        Fp::new(18323286026903235604),
    ],
    [
        Fp::new(10250026231324330997),
        Fp::new(13321947507807660157),
        Fp::new(13020725208899496943),
        Fp::new(11416990495425192684),
        Fp::new(7221795794796219413),
        Fp::new(2607917872900632985),
        Fp::new(2591896057192169329),
        Fp::new(10485489452304998145),
        Fp::new(9480186048908910015),
        Fp::new(2645141845409940474),
        Fp::new(16242299839765162610),
        Fp::new(12203738590896308135),
    ],
    [
        Fp::new(5395176197344543510),
        Fp::new(17941136338888340715),
        Fp::new(7559392505546762987),
        Fp::new(549633128904721280),
        Fp::new(15658455328409267684),
        Fp::new(10078371877170729592),
        Fp::new(2349868247408080783),
        Fp::new(13105911261634181239),
        Fp::new(12868653202234053626),
        Fp::new(9471330315555975806),
        Fp::new(4580289636625406680),
        Fp::new(13222733136951421572),
    ],
    [
        Fp::new(4555032575628627551),
        Fp::new(7619130111929922899),
        Fp::new(4547848507246491777),
        Fp::new(5662043532568004632),
        Fp::new(15723873049665279492),
        Fp::new(13585630674756818185),
        Fp::new(6990417929677264473),
        Fp::new(6373257983538884779),
        Fp::new(1005856792729125863),
        Fp::new(17850970025369572891),
        Fp::new(14306783492963476045),
        Fp::new(12653264875831356889),
    ],
    [
        Fp::new(10887434669785806501),
        Fp::new(5860438919639852469),
        Fp::new(13562758051366197240),
        Fp::new(11281423163572576941),
        Fp::new(11010364684624673712),
        Fp::new(3609894717198373475),
        Fp::new(2857078765033295647),
        Fp::new(10613907436992354031),
        Fp::new(2127516485693786130),
        Fp::new(5179216700721206500),
        Fp::new(667064457668387208),
        Fp::new(3269565065785113276),
    ],
    [
        Fp::new(10623812734732325494),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(3482313124719033742),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(16573951386448822347),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(11207203987960290077),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(15797451987878016994),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(8582621896649111280),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(17314968117136177458),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(662276588715934376),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(17440734083735621885),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(11157378678473576256),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(13379216091107483417),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(14817179324865344437),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(1846395309883290947),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(4870909675229475951),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(5661125543739831339),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(11127221369240904114),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(4250810027222144733),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(11760434366235360549),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(3330921671784295973),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(7122707995732409864),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(17522374822599693986),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::new(12723275943377720767),
        Fp::new(14785736031955679545),
        Fp::new(15257683393549924851),
        Fp::new(14586462537439744229),
        Fp::new(13109892360729616102),
        Fp::new(18054952537889795742),
        Fp::new(12589969976105374274),
        Fp::new(1436163932748701916),
        Fp::new(14879322534176465619),
        Fp::new(17580838042056220468),
        Fp::new(17970300042937392952),
        Fp::new(1420156878331078790),
    ],
    [
        Fp::new(17310902395782251544),
        Fp::new(9021117459098865178),
        Fp::new(9956374953785489337),
        Fp::new(9283926179170577664),
        Fp::new(2866744588122882663),
        Fp::new(12613310502798528952),
        Fp::new(48642999969593367),
        Fp::new(5069344854700671784),
        Fp::new(17704314310866354161),
        Fp::new(15988800480645163458),
        Fp::new(5818851986787837003),
        Fp::new(2578102338873304736),
    ],
    [
        Fp::new(13779049950938903889),
        Fp::new(9431620200784999174),
        Fp::new(16739974985861420869),
        Fp::new(15052263707393285803),
        Fp::new(11794044470831060964),
        Fp::new(12407310182639291233),
        Fp::new(16772957376150280651),
        Fp::new(7712582098709140624),
        Fp::new(1838118114767177406),
        Fp::new(14993893147640414015),
        Fp::new(1508554399796017280),
        Fp::new(14168503072706820248),
    ],
    [
        Fp::new(5465099181299765844),
        Fp::new(5157031041920728309),
        Fp::new(6420136306126768660),
        Fp::new(16331322219387907554),
        Fp::new(14945672812425851725),
        Fp::new(10584746352376021094),
        Fp::new(17049935967678990406),
        Fp::new(15218583273832305844),
        Fp::new(8057411706461231029),
        Fp::new(13457234708307314114),
        Fp::new(12805943922448093993),
        Fp::new(4260119674133082996),
    ],
];
