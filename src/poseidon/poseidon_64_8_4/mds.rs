// Copyright (c) 2021-2023 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use super::{NUM_PARTIAL_ROUNDS, STATE_WIDTH};
use cheetah::Fp;

/// Maximum Distance Separable matrix for Poseidon.
pub(crate) const MDS: [Fp; STATE_WIDTH * STATE_WIDTH] = [
    Fp::new(15665739677678329448),
    Fp::new(14761037431676564634),
    Fp::new(15232655799733888386),
    Fp::new(7219896008038694144),
    Fp::new(16833143609560954478),
    Fp::new(13371594715003915405),
    Fp::new(11025193874933309488),
    Fp::new(2115462038151941103),
    Fp::new(14714234413772989334),
    Fp::new(15709239302486602610),
    Fp::new(11400395424283092088),
    Fp::new(3307336641928453331),
    Fp::new(2083650131756605006),
    Fp::new(1900217460120204977),
    Fp::new(16949700886266222431),
    Fp::new(11823490811171925303),
    Fp::new(15295442582497673479),
    Fp::new(17396583668478712252),
    Fp::new(1343987031226794836),
    Fp::new(2339730835464162424),
    Fp::new(7413931974161992290),
    Fp::new(13534370939837294384),
    Fp::new(9832848613173085858),
    Fp::new(7844357158177642373),
    Fp::new(2696925559007867278),
    Fp::new(1557049008711126064),
    Fp::new(7360817417939431019),
    Fp::new(11078027825891222598),
    Fp::new(6607339708200885002),
    Fp::new(4684539808108992398),
    Fp::new(4717760937828897088),
    Fp::new(15964374282815256859),
    Fp::new(17965905919598481955),
    Fp::new(15407633716183131496),
    Fp::new(28955141318220746),
    Fp::new(13479003427736828777),
    Fp::new(4160715655911201686),
    Fp::new(11201498084541442593),
    Fp::new(8917207575607047436),
    Fp::new(7116052015563793652),
    Fp::new(9295364978535393065),
    Fp::new(4817503758083861761),
    Fp::new(16447695158761396496),
    Fp::new(4753576845484507428),
    Fp::new(7235543102394329024),
    Fp::new(7498238539141902111),
    Fp::new(9331742461781196439),
    Fp::new(12949701259399224346),
    Fp::new(15256992577437847297),
    Fp::new(6652925989773138817),
    Fp::new(5079847034286052344),
    Fp::new(4550811009404289617),
    Fp::new(1458397351619239493),
    Fp::new(17792684201245919421),
    Fp::new(6150173047518599418),
    Fp::new(1994521646738237037),
    Fp::new(3158367583368769003),
    Fp::new(14738065501122616700),
    Fp::new(301916385472689059),
    Fp::new(9357604073739055570),
    Fp::new(12675749649111885170),
    Fp::new(12491452512285176950),
    Fp::new(4730418877145812395),
    Fp::new(2086640369749051342),
];

pub(crate) const M_I: [[Fp; STATE_WIDTH]; STATE_WIDTH] = [
    [
        Fp::one(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::zero(),
        Fp::new(11417175008092614480),
        Fp::new(539339684927131214),
        Fp::new(9370096670683349538),
        Fp::new(17735590064681481672),
        Fp::new(15860523797886301402),
        Fp::new(331904639706794286),
        Fp::new(14847318717282070426),
    ],
    [
        Fp::zero(),
        Fp::new(8073066592479890593),
        Fp::new(4958763416319432662),
        Fp::new(4346299147281121438),
        Fp::new(7156152544070263116),
        Fp::new(6441351369888798847),
        Fp::new(12302195703376707011),
        Fp::new(13775848424915007996),
    ],
    [
        Fp::zero(),
        Fp::new(1521336658734357176),
        Fp::new(9609797837462109879),
        Fp::new(8696838521930320976),
        Fp::new(16546384001894101779),
        Fp::new(11905134903355379932),
        Fp::new(2746812526880489727),
        Fp::new(8792827167878067852),
    ],
    [
        Fp::zero(),
        Fp::new(7688692089893502263),
        Fp::new(1506157853856276434),
        Fp::new(8605622456439300958),
        Fp::new(2185874837149726550),
        Fp::new(2505224718279382174),
        Fp::new(5732863445774430136),
        Fp::new(2956799298581338507),
    ],
    [
        Fp::zero(),
        Fp::new(123809535153350209),
        Fp::new(5881371791636834400),
        Fp::new(12096906975260248744),
        Fp::new(9083397940013258929),
        Fp::new(5888207532787990280),
        Fp::new(653300070183740841),
        Fp::new(8374231230061881601),
    ],
    [
        Fp::zero(),
        Fp::new(16893340412411273511),
        Fp::new(15694659029239593053),
        Fp::new(5600689532757754456),
        Fp::new(14985967037687339922),
        Fp::new(7477095283177689305),
        Fp::new(2098213422797841468),
        Fp::new(6884501437450402973),
    ],
    [
        Fp::zero(),
        Fp::new(9165715784062500506),
        Fp::new(10356878312485286851),
        Fp::new(13160721654300935699),
        Fp::new(4814244272938750345),
        Fp::new(12718611369423834193),
        Fp::new(13484733700056817242),
        Fp::new(17942021740832824),
    ],
];

pub(crate) const V_COL: [[Fp; STATE_WIDTH]; NUM_PARTIAL_ROUNDS] = [
    [
        Fp::zero(),
        Fp::new(14714234413772989334),
        Fp::new(15295442582497673479),
        Fp::new(2696925559007867278),
        Fp::new(17965905919598481955),
        Fp::new(9295364978535393065),
        Fp::new(15256992577437847297),
        Fp::new(3158367583368769003),
    ],
    [
        Fp::zero(),
        Fp::new(17434648884499039172),
        Fp::new(15395171557141269351),
        Fp::new(17165130218314876351),
        Fp::new(2318989611215810392),
        Fp::new(17983354292234754098),
        Fp::new(3329479811376337598),
        Fp::new(7215420945515093878),
    ],
    [
        Fp::zero(),
        Fp::new(17396816530289318624),
        Fp::new(10818042536989618729),
        Fp::new(3286812612465662602),
        Fp::new(10551392076678634864),
        Fp::new(8239914318844963640),
        Fp::new(3708392512710007349),
        Fp::new(4211188505236435417),
    ],
    [
        Fp::zero(),
        Fp::new(13920731545234664628),
        Fp::new(13202538296811479702),
        Fp::new(4442397242288509728),
        Fp::new(12044996773503518821),
        Fp::new(10667432157756740067),
        Fp::new(7861826730005481550),
        Fp::new(1708356356398162853),
    ],
    [
        Fp::zero(),
        Fp::new(7953196718822256178),
        Fp::new(11664155898882213900),
        Fp::new(7885959015860783051),
        Fp::new(5192182123754012307),
        Fp::new(4009807460571720306),
        Fp::new(7815984065605510937),
        Fp::new(8757968947205464174),
    ],
    [
        Fp::zero(),
        Fp::new(9287716530172836684),
        Fp::new(9864232958758783932),
        Fp::new(12630587719892330711),
        Fp::new(13491890605628657297),
        Fp::new(3671628873646733547),
        Fp::new(3958394425879049514),
        Fp::new(13189122458756487660),
    ],
    [
        Fp::zero(),
        Fp::new(7878426467880568926),
        Fp::new(2458143552811094413),
        Fp::new(1864114828393993565),
        Fp::new(13854150286571093044),
        Fp::new(16922577471080933344),
        Fp::new(289713472824945579),
        Fp::new(3660794549807575902),
    ],
    [
        Fp::zero(),
        Fp::new(8011156089730631401),
        Fp::new(12192090446642039270),
        Fp::new(13531436176528423855),
        Fp::new(18279402914315908450),
        Fp::new(12122263267233787223),
        Fp::new(4649039760877044545),
        Fp::new(774117055211209308),
    ],
    [
        Fp::zero(),
        Fp::new(10187729073612718515),
        Fp::new(9831215990813509174),
        Fp::new(15894327751972047367),
        Fp::new(17678208387410779102),
        Fp::new(12370724615153314357),
        Fp::new(5612049686543101843),
        Fp::new(12941189046191301012),
    ],
    [
        Fp::zero(),
        Fp::new(259529520163092294),
        Fp::new(6678737766211136590),
        Fp::new(971968362182839516),
        Fp::new(5758641078030592869),
        Fp::new(6353776433509552482),
        Fp::new(4650461526464300370),
        Fp::new(3137875695092797832),
    ],
    [
        Fp::zero(),
        Fp::new(1373803839828072035),
        Fp::new(8722324625935338590),
        Fp::new(15183481750700167211),
        Fp::new(952948248534488297),
        Fp::new(17156331445851101253),
        Fp::new(15437010288959432158),
        Fp::new(3685347770567811932),
    ],
    [
        Fp::zero(),
        Fp::new(13153581878748812938),
        Fp::new(4435858640749877085),
        Fp::new(5536862541576354420),
        Fp::new(2697347949762671413),
        Fp::new(10272913315251174916),
        Fp::new(11305085310113155374),
        Fp::new(14265925088526361110),
    ],
    [
        Fp::zero(),
        Fp::new(4620840999367217461),
        Fp::new(9134293570234229924),
        Fp::new(17166279391362831041),
        Fp::new(2494769171254721370),
        Fp::new(16664561206228220272),
        Fp::new(5284489402954674742),
        Fp::new(12312420844101303761),
    ],
    [
        Fp::zero(),
        Fp::new(1654629273931315435),
        Fp::new(18254965379032004659),
        Fp::new(8880869464815497145),
        Fp::new(231324553483287115),
        Fp::new(6068831571199242158),
        Fp::new(14676454573891149722),
        Fp::new(8173461891987916228),
    ],
    [
        Fp::zero(),
        Fp::new(3315127955573276263),
        Fp::new(1051471900491075124),
        Fp::new(16657921117706684514),
        Fp::new(15588164500884865995),
        Fp::new(9432351979821966671),
        Fp::new(146172509861691911),
        Fp::new(14381341306358726870),
    ],
    [
        Fp::zero(),
        Fp::new(17052097107742048020),
        Fp::new(7438923376739246340),
        Fp::new(8135064217904846460),
        Fp::new(18282063799536637630),
        Fp::new(1629375766535839618),
        Fp::new(14187311274670463421),
        Fp::new(9014555022063543799),
    ],
    [
        Fp::zero(),
        Fp::new(15559204037314418191),
        Fp::new(2890334797296737012),
        Fp::new(4191660376196357349),
        Fp::new(16018664954818099696),
        Fp::new(7983751776208009484),
        Fp::new(12651490725524204962),
        Fp::new(155897827789606440),
    ],
    [
        Fp::zero(),
        Fp::new(3367423674765160956),
        Fp::new(12928436102312095142),
        Fp::new(7456338901332129971),
        Fp::new(17361978922988111419),
        Fp::new(4559459950508901974),
        Fp::new(1917836068719068979),
        Fp::new(1576114184934801301),
    ],
    [
        Fp::zero(),
        Fp::new(3475085957796945797),
        Fp::new(3295849077756365516),
        Fp::new(9438697043878875672),
        Fp::new(1503978507612151321),
        Fp::new(3312356440497309665),
        Fp::new(31143325588064814),
        Fp::new(11788246494843619827),
    ],
    [
        Fp::zero(),
        Fp::new(4838430657107732788),
        Fp::new(9355951818074364835),
        Fp::new(15403830109779901164),
        Fp::new(3884072746597139146),
        Fp::new(8134550246837891083),
        Fp::new(10206642302609506686),
        Fp::new(5647842928246683832),
    ],
    [
        Fp::zero(),
        Fp::new(17022142438444346336),
        Fp::new(17152505894966579525),
        Fp::new(2482349221949635510),
        Fp::new(11060929051930900139),
        Fp::new(17070005525500714889),
        Fp::new(15690262577663963466),
        Fp::new(7569882159656954547),
    ],
    [
        Fp::zero(),
        Fp::new(14435535831228617659),
        Fp::new(1796720195831182232),
        Fp::new(1442360747831927003),
        Fp::new(3368328654499589437),
        Fp::new(1606227434279650234),
        Fp::new(9548960172372227806),
        Fp::new(9298681339224169477),
    ],
];

pub(crate) const W_HAT: [[Fp; STATE_WIDTH]; NUM_PARTIAL_ROUNDS] = [
    [
        Fp::new(15665739677678329448),
        Fp::new(17232883515396795863),
        Fp::new(14236601672255005261),
        Fp::new(6833728164056015567),
        Fp::new(16744272745294593339),
        Fp::new(977502436806863166),
        Fp::new(7736136323470716757),
        Fp::new(6219350691475573577),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(4661801157882984172),
        Fp::new(6580924646000666874),
        Fp::new(1491456614677402743),
        Fp::new(17620039427400252496),
        Fp::new(9577938345332203423),
        Fp::new(8553962580890007161),
        Fp::new(14407132943251198789),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(314947281145154056),
        Fp::new(2639276743140882414),
        Fp::new(18359617890805587632),
        Fp::new(17720135834359893189),
        Fp::new(3717248477080227792),
        Fp::new(2271164900914575158),
        Fp::new(3645734884209347175),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(608875068351247300),
        Fp::new(10513948249862711204),
        Fp::new(16042957565259377754),
        Fp::new(4738996224031044393),
        Fp::new(9498844674137415389),
        Fp::new(14895916716216094503),
        Fp::new(4883015024678323663),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(2671808849846657829),
        Fp::new(18099494010101798878),
        Fp::new(9332824505110791188),
        Fp::new(5224322033233971607),
        Fp::new(12150193460642903704),
        Fp::new(8765357771888501148),
        Fp::new(7884162349068708932),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(18055148734663628248),
        Fp::new(9779014104285707040),
        Fp::new(10076064710657143718),
        Fp::new(18130396764255554272),
        Fp::new(10535916727797006992),
        Fp::new(13245845421496886075),
        Fp::new(6946854582558382830),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(8636118467569073221),
        Fp::new(9019105853311912199),
        Fp::new(9683845632922380100),
        Fp::new(16544164614968831326),
        Fp::new(9852843708879439779),
        Fp::new(6632939486909968000),
        Fp::new(13979627164078645610),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(10208810653387952974),
        Fp::new(1822699096543926793),
        Fp::new(3059924690378257357),
        Fp::new(8038098565083286089),
        Fp::new(4495069545760351456),
        Fp::new(3327203617447245207),
        Fp::new(10822129188448966135),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(3074903898833815658),
        Fp::new(5836424161120531724),
        Fp::new(16579797455531270547),
        Fp::new(3142347461449201457),
        Fp::new(16989524593149446515),
        Fp::new(13617954913757031979),
        Fp::new(2740237490072006335),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(12770176443797132592),
        Fp::new(8846000752338083651),
        Fp::new(10672378758952145557),
        Fp::new(2934048404118528784),
        Fp::new(5987294971168543875),
        Fp::new(18283816555037449308),
        Fp::new(17521646125737263736),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(10477798739121162779),
        Fp::new(6298727276747768874),
        Fp::new(11013708236169365910),
        Fp::new(11210191079330499899),
        Fp::new(507909967100413077),
        Fp::new(2171610989259822371),
        Fp::new(5232943736465415166),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(6783768545294738226),
        Fp::new(14678703841381030384),
        Fp::new(9983741824202831481),
        Fp::new(14856679798235258963),
        Fp::new(14689973280222977585),
        Fp::new(10409316311995350006),
        Fp::new(14067258760709817042),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(13448256280298132724),
        Fp::new(10276774410130051289),
        Fp::new(13925975751618149941),
        Fp::new(15772917982255413389),
        Fp::new(3842441561710749847),
        Fp::new(2879362729761045254),
        Fp::new(3249793397834684295),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(7410641574912520668),
        Fp::new(15320816675489988652),
        Fp::new(5657763507863891187),
        Fp::new(9119867548485666918),
        Fp::new(6314430726896704725),
        Fp::new(5711452648741414528),
        Fp::new(17602658866987624023),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(10994047073907868483),
        Fp::new(8216641614210664126),
        Fp::new(3503858270552599384),
        Fp::new(10011046567191683738),
        Fp::new(14932774705604500283),
        Fp::new(2316253901020668110),
        Fp::new(18149748225829126241),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(8375970984005269224),
        Fp::new(3435980036161356512),
        Fp::new(5767070434603677324),
        Fp::new(1452235317312856512),
        Fp::new(16136258063284525039),
        Fp::new(14733492472701417223),
        Fp::new(7183095870196989032),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(17557010166753490623),
        Fp::new(14028280845694952237),
        Fp::new(12914216367085308653),
        Fp::new(5822048435689966810),
        Fp::new(18146230864389865958),
        Fp::new(9849996067876521215),
        Fp::new(17831160909507072350),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(14343470352978796062),
        Fp::new(15918016519014119368),
        Fp::new(1517569989932517265),
        Fp::new(14183095349894257225),
        Fp::new(7981405543668413997),
        Fp::new(18236278366866574291),
        Fp::new(12383051004075138289),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(16814247010191196818),
        Fp::new(8247690514607887002),
        Fp::new(5331770784558034365),
        Fp::new(1305149629581187109),
        Fp::new(81613922796915882),
        Fp::new(17750390909352117992),
        Fp::new(13690790489809218463),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(18425006711959819293),
        Fp::new(2234547165584397360),
        Fp::new(187374063863358105),
        Fp::new(10249579720640525061),
        Fp::new(8676144487086415955),
        Fp::new(2263417139807844923),
        Fp::new(16320520364905964003),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(11071864438782463275),
        Fp::new(6496150729843188361),
        Fp::new(86943579747116357),
        Fp::new(10467289479100461068),
        Fp::new(7560708523540007019),
        Fp::new(9729076366842668834),
        Fp::new(219373675070972461),
    ],
    [
        Fp::new(15665739677678329448),
        Fp::new(10829072901710157659),
        Fp::new(11409850929191849698),
        Fp::new(12797502451450543895),
        Fp::new(9513101583707246677),
        Fp::new(11873701286632384758),
        Fp::new(2717809327157260447),
        Fp::new(1852895931045625524),
    ],
];
