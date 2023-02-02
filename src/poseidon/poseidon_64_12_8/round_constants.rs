// Copyright (c) 2021-2023 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::{NUM_HALF_FULL_ROUNDS, NUM_PARTIAL_ROUNDS, STATE_WIDTH};
use cheetah::Fp;

/// Additive Round Keys constants for Poseidon,
/// computed using algorithm 5 from <https://eprint.iacr.org/2020/1143.pdf>
pub(crate) const ARK: [[Fp; STATE_WIDTH]; 2 * NUM_HALF_FULL_ROUNDS + NUM_PARTIAL_ROUNDS] = [
    [
        Fp::new(1431286215153372998),
        Fp::new(3509349009260703107),
        Fp::new(2289575380984896342),
        Fp::new(10625215922958251110),
        Fp::new(17137022507167291684),
        Fp::new(17143426961497010024),
        Fp::new(9589775313463224365),
        Fp::new(7736066733515538648),
        Fp::new(2217569167061322248),
        Fp::new(10394930802584583083),
        Fp::new(4612393375016695705),
        Fp::new(5332470884919453534),
    ],
    [
        Fp::new(8724526834049581439),
        Fp::new(17673787971454860688),
        Fp::new(2519987773101056005),
        Fp::new(7999687124137420323),
        Fp::new(18312454652563306701),
        Fp::new(15136091233824155669),
        Fp::new(1257110570403430003),
        Fp::new(5665449074466664773),
        Fp::new(16178737609685266571),
        Fp::new(52855143527893348),
        Fp::new(8084454992943870230),
        Fp::new(2597062441266647183),
    ],
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
        Fp::new(7221072982690633460),
        Fp::new(9953585853856674407),
        Fp::new(13497620366078753434),
        Fp::new(18140292631504202243),
        Fp::new(17311934738088402529),
        Fp::new(6686302214424395771),
        Fp::new(11193071888943695519),
        Fp::new(10233795775801758543),
        Fp::new(3362219552562939863),
        Fp::new(8595401306696186761),
        Fp::new(7753411262943026561),
    ],
    [
        Fp::new(12415218859476220947),
        Fp::new(12517451587026875834),
        Fp::new(3257008032900598499),
        Fp::new(2187469039578904770),
        Fp::new(657675168296710415),
        Fp::new(8659969869470208989),
        Fp::new(12526098871288378639),
        Fp::new(12525853395769009329),
        Fp::new(15388161689979551704),
        Fp::new(7880966905416338909),
        Fp::new(2911694411222711481),
        Fp::new(6420652251792580406),
    ],
    [
        Fp::new(323544930728360053),
        Fp::new(11718666476052241225),
        Fp::new(2449132068789045592),
        Fp::new(17993014181992530560),
        Fp::new(15161788952257357966),
        Fp::new(3788504801066818367),
        Fp::new(1282111773460545571),
        Fp::new(8849495164481705550),
        Fp::new(8380852402060721190),
        Fp::new(2161980224591127360),
        Fp::new(2440151485689245146),
        Fp::new(17521895002090134367),
    ],
    [
        Fp::new(13821005335130766955),
        Fp::new(17513705631114265826),
        Fp::new(17068447856797239529),
        Fp::new(17964439003977043993),
        Fp::new(5685000919538239429),
        Fp::new(11615940660682589106),
        Fp::new(2522854885180605258),
        Fp::new(12584118968072796115),
        Fp::new(17841258728624635591),
        Fp::new(10821564568873127316),
        Fp::new(12929526205313074951),
        Fp::new(15240209309138869842),
    ],
    [
        Fp::new(8112988184280322821),
        Fp::new(10264318651796760217),
        Fp::new(11567563749053508498),
        Fp::new(10342172001635729828),
        Fp::new(8518076871621000645),
        Fp::new(9443305710168864155),
        Fp::new(12258139284331692775),
        Fp::new(11225713976478342221),
        Fp::new(1083829959428202152),
        Fp::new(13295679221277307734),
        Fp::new(8702942527907868190),
        Fp::new(3447159893350309030),
    ],
    [
        Fp::new(16331987863400672412),
        Fp::new(17004721198375099349),
        Fp::new(14568842036851006853),
        Fp::new(14031093640500276073),
        Fp::new(8047796853787800360),
        Fp::new(18176470296573070531),
        Fp::new(1733280390763076136),
        Fp::new(15280460251950617888),
        Fp::new(5319165528697198957),
        Fp::new(4130010739946422935),
        Fp::new(4862639442103099490),
        Fp::new(11947225653897253435),
    ],
    [
        Fp::new(16093634485870170562),
        Fp::new(466101267687143357),
        Fp::new(5269775209624779324),
        Fp::new(12661180512164132421),
        Fp::new(8527855600080265358),
        Fp::new(3509637282341164493),
        Fp::new(14524011473168972347),
        Fp::new(9558935312509120777),
        Fp::new(8282858737521047195),
        Fp::new(10171277103718892682),
        Fp::new(12294317531079789416),
        Fp::new(7182028925080765556),
    ],
    [
        Fp::new(2038954051047328382),
        Fp::new(1572125904757759485),
        Fp::new(6023737508444785880),
        Fp::new(8798428950960158590),
        Fp::new(1968909394335647758),
        Fp::new(16968160382228211614),
        Fp::new(32551027029362334),
        Fp::new(3205180815856999908),
        Fp::new(10740246361676213188),
        Fp::new(10169158339754762156),
        Fp::new(15226715702476100867),
        Fp::new(8966100427867584251),
    ],
    [
        Fp::new(17917233579925756683),
        Fp::new(7959268962897120034),
        Fp::new(532408456989891872),
        Fp::new(9851667167813963284),
        Fp::new(13448506932345489306),
        Fp::new(16135486720253939622),
        Fp::new(8458050899770540390),
        Fp::new(6021254166081897382),
        Fp::new(15552837092683737625),
        Fp::new(15440505484365682848),
        Fp::new(16088056409693275462),
        Fp::new(6169635475476966421),
    ],
    [
        Fp::new(5480704578777097169),
        Fp::new(7516526247262867111),
        Fp::new(3438140470099985472),
        Fp::new(13048600081642942971),
        Fp::new(9829255629799717904),
        Fp::new(17311489510949436164),
        Fp::new(15254947846872712175),
        Fp::new(5825939868327872570),
        Fp::new(850656437239379199),
        Fp::new(12619934071925039179),
        Fp::new(15233049780346247641),
        Fp::new(9298309061465962971),
    ],
    [
        Fp::new(741424706267005090),
        Fp::new(17203483336096778815),
        Fp::new(6919908349347460635),
        Fp::new(863377837517698584),
        Fp::new(11632281421519826770),
        Fp::new(17750153240261395489),
        Fp::new(14753366294352507072),
        Fp::new(12793355793496405427),
        Fp::new(16289545878058120229),
        Fp::new(6368259120071113126),
        Fp::new(4057875983396832839),
        Fp::new(13847225916600191037),
    ],
    [
        Fp::new(7872218736019578342),
        Fp::new(5426064199624116028),
        Fp::new(9479822711840773905),
        Fp::new(10634838597871962689),
        Fp::new(7081809782259040995),
        Fp::new(1440626909472018594),
        Fp::new(6603963598898808862),
        Fp::new(12662045888242770199),
        Fp::new(18036285107641934643),
        Fp::new(15828843208411476617),
        Fp::new(14102670999874605825),
        Fp::new(15585654191999307703),
    ],
    [
        Fp::new(940187017142450256),
        Fp::new(8747386241522630712),
        Fp::new(6750641561540124748),
        Fp::new(7440998025584530008),
        Fp::new(6136358134615751537),
        Fp::new(12413576830284969612),
        Fp::new(11675438539028694710),
        Fp::new(17580553691069642927),
        Fp::new(892707462476851332),
        Fp::new(15167485180850043745),
        Fp::new(9924997173903409412),
        Fp::new(9613966396549972013),
    ],
    [
        Fp::new(3242363036477934858),
        Fp::new(8529581814542674199),
        Fp::new(1460135031320476117),
        Fp::new(15230276901939640657),
        Fp::new(3034222759280296577),
        Fp::new(2536834233629877234),
        Fp::new(12229748406346543211),
        Fp::new(13166855996952940567),
        Fp::new(16039201196582061794),
        Fp::new(14239610657545203244),
        Fp::new(4079052969819075917),
        Fp::new(2550303736432259954),
    ],
    [
        Fp::new(15415646525902701306),
        Fp::new(16984207496990988313),
        Fp::new(6195489392633771043),
        Fp::new(15696991486732177869),
        Fp::new(17238905290121258980),
        Fp::new(16082743896956175460),
        Fp::new(2607127875797716838),
        Fp::new(4163972359010584653),
        Fp::new(2369705041192477687),
        Fp::new(12936899802672086396),
        Fp::new(17399492193998111961),
        Fp::new(14701188996710188063),
    ],
    [
        Fp::new(10673647621461954174),
        Fp::new(10187656820932330866),
        Fp::new(14253604578356758004),
        Fp::new(10632764261170436503),
        Fp::new(2575456097595068268),
        Fp::new(14486510292332525540),
        Fp::new(857634655205127854),
        Fp::new(11539936742927634064),
        Fp::new(3025473245387650600),
        Fp::new(3072205393568168823),
        Fp::new(16220766505279212230),
        Fp::new(13095270286885528495),
    ],
    [
        Fp::new(10043771903993878423),
        Fp::new(4580450255883541632),
        Fp::new(5546821308061729354),
        Fp::new(17932404490144193348),
        Fp::new(4055843989895157237),
        Fp::new(506731346742428544),
        Fp::new(1750774988219982266),
        Fp::new(13647783723546009630),
        Fp::new(17180411145007510672),
        Fp::new(7092939346849547588),
        Fp::new(2004811345434270086),
        Fp::new(3930380885080085231),
    ],
    [
        Fp::new(5731056810399963425),
        Fp::new(16339249658689415041),
        Fp::new(10896947625319492019),
        Fp::new(58048537304546191),
        Fp::new(12301681553475871944),
        Fp::new(15410898306178483444),
        Fp::new(5248513067045859782),
        Fp::new(11268429244640014487),
        Fp::new(3785322258417388297),
        Fp::new(12573604913857968925),
        Fp::new(10088460126056383905),
        Fp::new(9505879368173225761),
    ],
    [
        Fp::new(12331335364636844807),
        Fp::new(15800425329127532993),
        Fp::new(17233569579365152217),
        Fp::new(6580598753390726049),
        Fp::new(6332388716747236070),
        Fp::new(14837976254465985338),
        Fp::new(1387653002144476724),
        Fp::new(15556347971769261667),
        Fp::new(7571094906243962853),
        Fp::new(14097015672565897063),
        Fp::new(1689918468007574312),
        Fp::new(16247594734699408053),
    ],
    [
        Fp::new(6376995477333092352),
        Fp::new(962981388472387485),
        Fp::new(2846128944153513179),
        Fp::new(11832408739941285626),
        Fp::new(16892791912968591653),
        Fp::new(14660122210495197643),
        Fp::new(16446079849332856874),
        Fp::new(7976724875926637635),
        Fp::new(13842280498640749771),
        Fp::new(15375657835094741734),
        Fp::new(8871752519026737048),
        Fp::new(6979293996243387512),
    ],
    [
        Fp::new(10552448846206288151),
        Fp::new(14987673924494666433),
        Fp::new(18035303280469462414),
        Fp::new(16595113834715919465),
        Fp::new(15208661533916677630),
        Fp::new(4170608138187333497),
        Fp::new(16304084357983152470),
        Fp::new(2331503858766652994),
        Fp::new(8776079357547932587),
        Fp::new(18299646478835171989),
        Fp::new(3681263166902989193),
        Fp::new(12612029705709390274),
    ],
    [
        Fp::new(12014669431902405777),
        Fp::new(11319504285297576766),
        Fp::new(5234999940078631477),
        Fp::new(1125448944938006422),
        Fp::new(2164405204907480972),
        Fp::new(6168495504522907053),
        Fp::new(6250236942243891229),
        Fp::new(18269902991411124149),
        Fp::new(9426885685329917236),
        Fp::new(4521800374915508165),
        Fp::new(2213719649464492152),
        Fp::new(9422759956003735939),
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
];
