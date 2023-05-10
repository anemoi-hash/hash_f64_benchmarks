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
    Fp::new(8062354434719430557),
    Fp::new(14961898876701168261),
    Fp::new(14072440978656289413),
    Fp::new(3813898698085822056),
    Fp::new(2661992107161659435),
    Fp::new(3094630301859670042),
    Fp::new(3977106272425873496),
    Fp::new(13933779115385017346),
    Fp::new(14189793454789259666),
    Fp::new(4589564893138867870),
    Fp::new(10656609797473742311),
    Fp::new(4351795636554906703),
    Fp::new(3629225075051654566),
    Fp::new(5505041507055232084),
    Fp::new(10613583062517235368),
    Fp::new(1959763024107998211),
    Fp::new(7567922212435485310),
    Fp::new(2733130895475987737),
    Fp::new(4343102577131897859),
    Fp::new(17270065735636854691),
    Fp::new(3752691483242407468),
    Fp::new(5081949864520279640),
    Fp::new(8073294779762788586),
    Fp::new(13956275920678724480),
    Fp::new(6376192859023102565),
    Fp::new(12742498195983677857),
    Fp::new(14773641347982789813),
    Fp::new(14837508958922443859),
    Fp::new(16763356155844416005),
    Fp::new(10598503663924302579),
    Fp::new(16184038293031944329),
    Fp::new(13621822751103418285),
    Fp::new(4305242994748865404),
    Fp::new(7879428323599719493),
    Fp::new(14047355934964505240),
    Fp::new(2705439982723035579),
    Fp::new(15291615095972091554),
    Fp::new(13266027256011357758),
    Fp::new(7833042914132636033),
    Fp::new(2067639338136971208),
    Fp::new(9895526792122960223),
    Fp::new(10209551385997036050),
    Fp::new(8773753981941606356),
    Fp::new(13205642198100979132),
    Fp::new(7837802265059468963),
    Fp::new(6408204717473862567),
    Fp::new(15830614506969495927),
    Fp::new(15013501726694741751),
    Fp::new(4231867082619939226),
    Fp::new(9487199445522374556),
    Fp::new(10193853044559582482),
    Fp::new(13922939181059841659),
    Fp::new(17666289061339892708),
    Fp::new(10897134177118973490),
    Fp::new(16386144691484835842),
    Fp::new(1072548541978236852),
    Fp::new(6069893080856189971),
    Fp::new(8881735842565329257),
    Fp::new(5850024874013740724),
    Fp::new(6976954589483661694),
    Fp::new(4184087375896651807),
    Fp::new(3638918904749292497),
    Fp::new(17239574026122101423),
    Fp::new(12786639009708118542),
    Fp::new(6254862572315590282),
    Fp::new(1629590728702708772),
    Fp::new(6673060922912979208),
    Fp::new(6516195972713186877),
    Fp::new(16915033719725765719),
    Fp::new(4450092469307778592),
    Fp::new(6146395711041360497),
    Fp::new(12203405532558077535),
    Fp::new(3447947061680092303),
    Fp::new(14698586745078822058),
    Fp::new(7417801860138218805),
    Fp::new(7688875089822095986),
    Fp::new(2037829315856691427),
    Fp::new(18012203483038413286),
    Fp::new(10134743567522000194),
    Fp::new(17690756460959553601),
    Fp::new(18192098666773533285),
    Fp::new(6311149141618075588),
    Fp::new(14728574619499559417),
    Fp::new(11983667693583105106),
    Fp::new(2569967209706010234),
    Fp::new(7580063482011404161),
    Fp::new(1037264167681961852),
    Fp::new(8606342303052326554),
    Fp::new(14674322088923028330),
    Fp::new(17504536570532066534),
    Fp::new(8146281342193668769),
    Fp::new(9557835506920111885),
    Fp::new(850433035699737413),
    Fp::new(11898265417469071053),
    Fp::new(7779157783667474031),
    Fp::new(8995535601867580306),
    Fp::new(7581636945872566763),
    Fp::new(14601060572455277959),
    Fp::new(3016935651704872944),
    Fp::new(7026634763238041693),
    Fp::new(196656516477314644),
    Fp::new(12243723558789076687),
    Fp::new(6216628555349727933),
    Fp::new(2072871153194036346),
    Fp::new(13605384252814469191),
    Fp::new(11815698918607316697),
    Fp::new(9467288609154657287),
    Fp::new(15812850586328277755),
    Fp::new(4316582198064733503),
    Fp::new(4791006047222050110),
    Fp::new(1145519502995968566),
    Fp::new(4285137197096394109),
    Fp::new(9117066297380536120),
    Fp::new(9054380764130329792),
    Fp::new(17228068624121884666),
    Fp::new(11082693959710948890),
    Fp::new(16616264979856204026),
    Fp::new(12374626663989955442),
    Fp::new(9731001592189262557),
    Fp::new(2396964613666436532),
    Fp::new(16201557689612905826),
    Fp::new(125361537117473690),
    Fp::new(15830772518552943808),
    Fp::new(10511623029446172327),
    Fp::new(10570957911041485345),
    Fp::new(4962838094393792127),
    Fp::new(9076005812400608041),
    Fp::new(12254095523339714984),
    Fp::new(12890969368040571960),
    Fp::new(17094410138259985087),
    Fp::new(6696785618317865133),
    Fp::new(15368358043273625229),
    Fp::new(5140309616425735503),
    Fp::new(12788654981167798303),
    Fp::new(16891596484841155792),
    Fp::new(12991174693007313597),
    Fp::new(8473418439344353796),
    Fp::new(5691475968038033816),
    Fp::new(16396843216027926601),
    Fp::new(14015446524163305690),
    Fp::new(1321923216714856411),
    Fp::new(12535014903554896842),
    Fp::new(3217272198733985269),
    Fp::new(15190344371695839041),
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
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
        Fp::zero(),
    ],
    [
        Fp::zero(),
        Fp::new(2158060718382483437),
        Fp::new(11754347615332992199),
        Fp::new(3486218395374702936),
        Fp::new(4653298637352379119),
        Fp::new(1953963121003461671),
        Fp::new(15548184796780012808),
        Fp::new(4543539838695305424),
        Fp::new(5784489872872352520),
        Fp::new(3415021772692410798),
        Fp::new(805230818572575460),
        Fp::new(6774082512037881392),
    ],
    [
        Fp::zero(),
        Fp::new(12633424942224798421),
        Fp::new(15221410160123011287),
        Fp::new(7048571406563264703),
        Fp::new(16666615639444866792),
        Fp::new(10851184604049541589),
        Fp::new(2902262303638539229),
        Fp::new(14593589866004278817),
        Fp::new(15471886733357500598),
        Fp::new(3489621560158593854),
        Fp::new(12267358463494859140),
        Fp::new(16819193817668948100),
    ],
    [
        Fp::zero(),
        Fp::new(4777031255670923347),
        Fp::new(9576416882454739319),
        Fp::new(3390695810356745466),
        Fp::new(13610041889433331621),
        Fp::new(4260292174967232937),
        Fp::new(15694280630349127258),
        Fp::new(14832434754214995386),
        Fp::new(9951348002575169362),
        Fp::new(11056398025294547380),
        Fp::new(3981890029904108424),
        Fp::new(9931598532857996089),
    ],
    [
        Fp::zero(),
        Fp::new(13810382489330343170),
        Fp::new(4934201052230350224),
        Fp::new(17815710998065970247),
        Fp::new(16508747845641758791),
        Fp::new(2093595995981523053),
        Fp::new(11933460275235303455),
        Fp::new(10187121468219971046),
        Fp::new(10886010806916911815),
        Fp::new(1376074721579680681),
        Fp::new(11268169567552075960),
        Fp::new(17733071428694541110),
    ],
    [
        Fp::zero(),
        Fp::new(5980937736834687786),
        Fp::new(18131261854779456237),
        Fp::new(5683127066971111493),
        Fp::new(12940031247090567099),
        Fp::new(16923628384061636138),
        Fp::new(16770585064217543619),
        Fp::new(512996911318723683),
        Fp::new(10573764975764958572),
        Fp::new(13420511277149628271),
        Fp::new(16174409121062643931),
        Fp::new(12072613736336094882),
    ],
    [
        Fp::zero(),
        Fp::new(17615625872697810809),
        Fp::new(3696372200709561186),
        Fp::new(16534123397177176779),
        Fp::new(7956713201180548613),
        Fp::new(5404224727068436824),
        Fp::new(16744644483687293548),
        Fp::new(10576287842559590775),
        Fp::new(3016791753848188),
        Fp::new(11653065982005406757),
        Fp::new(5222141619784919975),
        Fp::new(1286242358904612343),
    ],
    [
        Fp::zero(),
        Fp::new(8519081369454728786),
        Fp::new(1830975090686683544),
        Fp::new(4818144576565365187),
        Fp::new(5036568925971673080),
        Fp::new(12141104289011305014),
        Fp::new(13695162428859011932),
        Fp::new(13488838327375509042),
        Fp::new(5601639491406028618),
        Fp::new(16133371653697095553),
        Fp::new(7960809792654116594),
        Fp::new(6502657508317070596),
    ],
    [
        Fp::zero(),
        Fp::new(17983666362169740948),
        Fp::new(7402203711856260174),
        Fp::new(11045624522517011798),
        Fp::new(10058280257916306652),
        Fp::new(17570862541139196820),
        Fp::new(7152331555942404225),
        Fp::new(698534714526749786),
        Fp::new(8881011192109198104),
        Fp::new(12302470161068436372),
        Fp::new(10869283568618512040),
        Fp::new(9221307689218580481),
    ],
    [
        Fp::zero(),
        Fp::new(352717809370599872),
        Fp::new(15993806483102991232),
        Fp::new(11845781958513489049),
        Fp::new(18285750637631418779),
        Fp::new(8627836555235183197),
        Fp::new(2970007195368798264),
        Fp::new(2987548373835153756),
        Fp::new(9182460993923536403),
        Fp::new(17780027919379006106),
        Fp::new(13254853361609060086),
        Fp::new(10017444878623389708),
    ],
    [
        Fp::zero(),
        Fp::new(14596717606170371766),
        Fp::new(2194427887710972645),
        Fp::new(10972555231289128451),
        Fp::new(13191954341739574054),
        Fp::new(15970493132971852340),
        Fp::new(17943264901704182513),
        Fp::new(6907532238036774524),
        Fp::new(1462585920642101189),
        Fp::new(2126860165578362362),
        Fp::new(14340036452199296972),
        Fp::new(3614239723492323729),
    ],
    [
        Fp::zero(),
        Fp::new(11878471173434419322),
        Fp::new(6116900358833301660),
        Fp::new(9618413108795116640),
        Fp::new(16034138042725216380),
        Fp::new(294786144583825969),
        Fp::new(1767322893549682670),
        Fp::new(18009420070725193469),
        Fp::new(2317779622249459788),
        Fp::new(5241564959198791370),
        Fp::new(10579988325410622082),
        Fp::new(16274346933542694497),
    ],
];

pub(crate) const V_COL: [[Fp; STATE_WIDTH - 1]; NUM_PARTIAL_ROUNDS] = [
    [
        Fp::new(3629225075051654566),
        Fp::new(6376192859023102565),
        Fp::new(15291615095972091554),
        Fp::new(4231867082619939226),
        Fp::new(4184087375896651807),
        Fp::new(3447947061680092303),
        Fp::new(2569967209706010234),
        Fp::new(7581636945872566763),
        Fp::new(4316582198064733503),
        Fp::new(16201557689612905826),
        Fp::new(5140309616425735503),
    ],
    [
        Fp::new(3174769152765695593),
        Fp::new(15398814493500280447),
        Fp::new(7102572464868618193),
        Fp::new(16994932202887483160),
        Fp::new(7595997084663333207),
        Fp::new(17194652927860183796),
        Fp::new(11871246013281161847),
        Fp::new(12655654669163645933),
        Fp::new(5597619770687835499),
        Fp::new(9527594309024257404),
        Fp::new(5774890393891431239),
    ],
    [
        Fp::new(9900823259986584006),
        Fp::new(7897512944977505901),
        Fp::new(6359896515385500498),
        Fp::new(17179422913435930974),
        Fp::new(4154948033283679269),
        Fp::new(5671013427478139457),
        Fp::new(3528074043479024517),
        Fp::new(18324059333920206712),
        Fp::new(1668112970542290857),
        Fp::new(3625261527297370832),
        Fp::new(13539537832154938145),
    ],
    [
        Fp::new(15950971046952987094),
        Fp::new(16757976457980953654),
        Fp::new(13162381019784039009),
        Fp::new(18313989005710292068),
        Fp::new(268453962807870821),
        Fp::new(6098660200059101590),
        Fp::new(17270759517455511357),
        Fp::new(1808547483243666174),
        Fp::new(10561631511550265701),
        Fp::new(5441948435397805079),
        Fp::new(13085113480963835196),
    ],
    [
        Fp::new(2065177704099479340),
        Fp::new(11627079021215808113),
        Fp::new(14903643111665525020),
        Fp::new(16369376864207675730),
        Fp::new(9060376496652906879),
        Fp::new(754550844407590986),
        Fp::new(14069983996522666276),
        Fp::new(12832955135057147400),
        Fp::new(6436779036950727791),
        Fp::new(3975477911503929212),
        Fp::new(8987970270225818452),
    ],
    [
        Fp::new(13139392767217129631),
        Fp::new(7911906441256040991),
        Fp::new(11274432502984851096),
        Fp::new(11211209709553037295),
        Fp::new(1372859574983125504),
        Fp::new(12397366047374959673),
        Fp::new(17772165147772940859),
        Fp::new(5084288824275860183),
        Fp::new(1738018633395983276),
        Fp::new(11843853300175599524),
        Fp::new(15472393889590690640),
    ],
    [
        Fp::new(3904941661289751125),
        Fp::new(16460825649328888370),
        Fp::new(10008761613534697254),
        Fp::new(17502493990320430191),
        Fp::new(4277916158331180526),
        Fp::new(884123356170170309),
        Fp::new(5523035328343569958),
        Fp::new(13794118152417037822),
        Fp::new(600111217372614462),
        Fp::new(7805094262647227503),
        Fp::new(17062759630687559556),
    ],
    [
        Fp::new(6461764041951615326),
        Fp::new(4656265152802223966),
        Fp::new(16046553876337806210),
        Fp::new(4514234704320432401),
        Fp::new(11091741279640103475),
        Fp::new(15062759977368613582),
        Fp::new(3868179120167499416),
        Fp::new(7041075066228536609),
        Fp::new(15622226894059770996),
        Fp::new(16883668369804892660),
        Fp::new(3675545581701288072),
    ],
    [
        Fp::new(3114161665082839122),
        Fp::new(12753443089443893759),
        Fp::new(688046305758890303),
        Fp::new(9724916664503296832),
        Fp::new(9949803546411180303),
        Fp::new(3234289782707954702),
        Fp::new(7206596811010183353),
        Fp::new(12458667016285968957),
        Fp::new(8041024608611587254),
        Fp::new(6085249416353960285),
        Fp::new(5631414445789734936),
    ],
    [
        Fp::new(14328022145509946308),
        Fp::new(10281273059164539845),
        Fp::new(10953870767391011328),
        Fp::new(9215264849880260265),
        Fp::new(739903065261847740),
        Fp::new(14217092720568361632),
        Fp::new(4476233894250140047),
        Fp::new(1078026271423843718),
        Fp::new(11738310925408701982),
        Fp::new(13229892590184962812),
        Fp::new(5491076495390020912),
    ],
    [
        Fp::new(5881986625426797506),
        Fp::new(6938819412898842),
        Fp::new(3192263132204025985),
        Fp::new(564404582539101496),
        Fp::new(9598341869352134181),
        Fp::new(6322739125080706197),
        Fp::new(278711029585087508),
        Fp::new(3560662468400468052),
        Fp::new(10094244102135883495),
        Fp::new(10368866570483793634),
        Fp::new(18173424469866095558),
    ],
    [
        Fp::new(4340233673262743195),
        Fp::new(5412520373127784619),
        Fp::new(15087018169083941570),
        Fp::new(9556794505221451448),
        Fp::new(4364021002675748317),
        Fp::new(2104067839451465804),
        Fp::new(11074071509739717784),
        Fp::new(9928047471811309533),
        Fp::new(1820748710299608652),
        Fp::new(11197934435654809429),
        Fp::new(5892729950681645780),
    ],
    [
        Fp::new(1976645583306580231),
        Fp::new(10879448862362332422),
        Fp::new(3680486528272222032),
        Fp::new(1414763658050872129),
        Fp::new(15475980094393448281),
        Fp::new(13723475442576965524),
        Fp::new(4271432360104067112),
        Fp::new(3806075920426632338),
        Fp::new(2820403172525422582),
        Fp::new(2383048935092489737),
        Fp::new(11786889667190877716),
    ],
    [
        Fp::new(103100088301150079),
        Fp::new(5320097144863780227),
        Fp::new(987335612830152350),
        Fp::new(1552585381046976392),
        Fp::new(13882596209699906017),
        Fp::new(18093675644281219501),
        Fp::new(12596149597955621142),
        Fp::new(2070197435224727453),
        Fp::new(3680603003014812465),
        Fp::new(14112546989100657610),
        Fp::new(12494359290403337778),
    ],
    [
        Fp::new(14151934581786944899),
        Fp::new(13887699195350627988),
        Fp::new(8795178045181095299),
        Fp::new(4399117974048598099),
        Fp::new(10999089880759249935),
        Fp::new(11924725603586529716),
        Fp::new(8710459237669954437),
        Fp::new(12651238610945171587),
        Fp::new(4641159618944760377),
        Fp::new(6169141624008518477),
        Fp::new(4881059707347026494),
    ],
    [
        Fp::new(4342599146291414314),
        Fp::new(13580352481588642514),
        Fp::new(11795467523314822754),
        Fp::new(5886070139709419559),
        Fp::new(2508866944054020039),
        Fp::new(1630754566261052977),
        Fp::new(16468765799918809378),
        Fp::new(13413684424144493660),
        Fp::new(129989508857472649),
        Fp::new(9252185312076567640),
        Fp::new(1613341563578894196),
    ],
    [
        Fp::new(5881161376590569367),
        Fp::new(18238408201523176044),
        Fp::new(4814836150892889140),
        Fp::new(4613463776836364995),
        Fp::new(2090919146578941871),
        Fp::new(8502911456828685773),
        Fp::new(12506467312026599905),
        Fp::new(16445802965524940471),
        Fp::new(249437767226346575),
        Fp::new(1825205014635838860),
        Fp::new(7163429820956120073),
    ],
    [
        Fp::new(12793466016333952410),
        Fp::new(11431540265716897966),
        Fp::new(1775225793782912443),
        Fp::new(18281609648435992199),
        Fp::new(16227569315938784950),
        Fp::new(17049530573936359840),
        Fp::new(6147457812712741548),
        Fp::new(12686937074915559964),
        Fp::new(4925115237899128480),
        Fp::new(2579719555369486540),
        Fp::new(7182060210949979232),
    ],
    [
        Fp::new(12454228053357765563),
        Fp::new(2629692629805585783),
        Fp::new(1093612361041406535),
        Fp::new(6936414575345642244),
        Fp::new(8896009754285000024),
        Fp::new(7744378096686790110),
        Fp::new(15661543027297758554),
        Fp::new(825969605583424032),
        Fp::new(3675939254525501555),
        Fp::new(13832550918826371654),
        Fp::new(7477046261159793357),
    ],
    [
        Fp::new(13004316706408310304),
        Fp::new(747875737245729263),
        Fp::new(17925664271468088999),
        Fp::new(7833034556437495734),
        Fp::new(17267948274759122624),
        Fp::new(13456101099020500377),
        Fp::new(432486153647843301),
        Fp::new(8346576719598541686),
        Fp::new(466546354564520542),
        Fp::new(560654515341796946),
        Fp::new(10216086068922803667),
    ],
    [
        Fp::new(12495389817112187194),
        Fp::new(4997848048970644713),
        Fp::new(8820153561103573931),
        Fp::new(1135816623197994199),
        Fp::new(14182178139782796191),
        Fp::new(12918211244185358863),
        Fp::new(9975242396411442972),
        Fp::new(8469371067231849782),
        Fp::new(16784703395658719491),
        Fp::new(2358923861924523686),
        Fp::new(1856315545073838847),
    ],
    [
        Fp::new(10886701715421322959),
        Fp::new(1160652906737679474),
        Fp::new(11568002309872769775),
        Fp::new(6391424287152370448),
        Fp::new(2586034218968342471),
        Fp::new(15723735635326701099),
        Fp::new(9220105308603778506),
        Fp::new(4976751516332797122),
        Fp::new(7240687926687836112),
        Fp::new(5404232401386714988),
        Fp::new(12029411450680498427),
    ],
];

pub(crate) const W_HAT: [[Fp; STATE_WIDTH - 1]; NUM_PARTIAL_ROUNDS] = [
    [
        Fp::new(1867177648889215775),
        Fp::new(14047986573916550888),
        Fp::new(13678177067582171850),
        Fp::new(12863463961918205696),
        Fp::new(2816458212870361606),
        Fp::new(6488135874391666051),
        Fp::new(3437116424102173320),
        Fp::new(8842687689188467876),
        Fp::new(2826064837365678581),
        Fp::new(658840792183071885),
        Fp::new(15721493965206350769),
    ],
    [
        Fp::new(16597666069766672862),
        Fp::new(18115039219464537806),
        Fp::new(15402418414365514836),
        Fp::new(15881049221707682123),
        Fp::new(44038082423104488),
        Fp::new(493477868698632600),
        Fp::new(16593990050305803718),
        Fp::new(7282461366904484602),
        Fp::new(14013131729240129536),
        Fp::new(16633394452620714095),
        Fp::new(12311998408006928744),
    ],
    [
        Fp::new(2907512777859410844),
        Fp::new(9208503536965406982),
        Fp::new(16216374438449666965),
        Fp::new(14844938336004632567),
        Fp::new(5789190283744873881),
        Fp::new(15342158928834710541),
        Fp::new(662750964240076544),
        Fp::new(2425807139631389546),
        Fp::new(1662654065826901278),
        Fp::new(3289334310811262971),
        Fp::new(8248327568107419799),
    ],
    [
        Fp::new(3111377728668235695),
        Fp::new(6742868611572023696),
        Fp::new(14104621223398144042),
        Fp::new(8154768093785310338),
        Fp::new(11546672114354462090),
        Fp::new(12907946023713321434),
        Fp::new(244172551034307786),
        Fp::new(11243434373698781387),
        Fp::new(16370042466624472008),
        Fp::new(7745439498859874252),
        Fp::new(3416622226921401933),
    ],
    [
        Fp::new(6262731050311674691),
        Fp::new(3626814920263147703),
        Fp::new(5492674111065556195),
        Fp::new(14458343832821467249),
        Fp::new(16253146333818890399),
        Fp::new(6451961014696003330),
        Fp::new(14751938965208982692),
        Fp::new(6861030538469410506),
        Fp::new(17351496915228422105),
        Fp::new(3600152177012496450),
        Fp::new(10396985109687211952),
    ],
    [
        Fp::new(10090613519320443604),
        Fp::new(14488742663342256606),
        Fp::new(6951290240380398652),
        Fp::new(4146252529551147990),
        Fp::new(13133229638559488566),
        Fp::new(6809328597228316785),
        Fp::new(5378433248986465834),
        Fp::new(4748552632287917540),
        Fp::new(9791220523410774835),
        Fp::new(13971468217297075439),
        Fp::new(11090756600293661206),
    ],
    [
        Fp::new(13837602407884432356),
        Fp::new(7834482161321069499),
        Fp::new(9357723951755562414),
        Fp::new(8902058118751616812),
        Fp::new(9972196027117396602),
        Fp::new(2906545461378228693),
        Fp::new(18188156607311196991),
        Fp::new(10705143547155513456),
        Fp::new(11061387832956172260),
        Fp::new(10634538347498592605),
        Fp::new(1493830521795341373),
    ],
    [
        Fp::new(16624695864080499990),
        Fp::new(13791686741552249569),
        Fp::new(9823492037892770377),
        Fp::new(13482883000361194719),
        Fp::new(17659478796588298971),
        Fp::new(15087807004491322280),
        Fp::new(1958900311927974762),
        Fp::new(18221749351506921616),
        Fp::new(15971448179451753804),
        Fp::new(8579073481743428462),
        Fp::new(3490664135870704446),
    ],
    [
        Fp::new(8870679208657801571),
        Fp::new(13291085073230109121),
        Fp::new(9445946220599352859),
        Fp::new(4068613378454613166),
        Fp::new(7219976475861648115),
        Fp::new(9425784200244693910),
        Fp::new(5608999798069142587),
        Fp::new(17658277199332239983),
        Fp::new(177468111938977865),
        Fp::new(5248759849419656971),
        Fp::new(9388947517913306928),
    ],
    [
        Fp::new(10863737689744328425),
        Fp::new(651616585337249751),
        Fp::new(1987616664402086318),
        Fp::new(11285689554159888461),
        Fp::new(14551926402005037496),
        Fp::new(10776802328952256806),
        Fp::new(3368501198988478445),
        Fp::new(3140142022946087091),
        Fp::new(1985580634208993405),
        Fp::new(6966742083457508697),
        Fp::new(4699718752057263567),
    ],
    [
        Fp::new(15256972163927412934),
        Fp::new(14660379870084621829),
        Fp::new(17904624720016941378),
        Fp::new(17454967393183206167),
        Fp::new(12178846075159167958),
        Fp::new(9032980769023847223),
        Fp::new(3920437407641323089),
        Fp::new(15097990686465748669),
        Fp::new(4616484942963685939),
        Fp::new(18117188739823719914),
        Fp::new(1422873100358434209),
    ],
    [
        Fp::new(17366439487366385578),
        Fp::new(1760120500336612355),
        Fp::new(10136299924184357869),
        Fp::new(4638856500017385281),
        Fp::new(42387822155267507),
        Fp::new(2369106805807710920),
        Fp::new(18359150410469708264),
        Fp::new(3117745207992482838),
        Fp::new(8735029656671880111),
        Fp::new(219126756713711933),
        Fp::new(15927989835793909982),
    ],
    [
        Fp::new(3016687123234865260),
        Fp::new(7168959499650507543),
        Fp::new(3621666702397166038),
        Fp::new(6790995893349037069),
        Fp::new(763872415753733329),
        Fp::new(982995134154506940),
        Fp::new(4345554329803862621),
        Fp::new(15097519625812440149),
        Fp::new(2977331789057996380),
        Fp::new(7403401748101264905),
        Fp::new(8642092307169954224),
    ],
    [
        Fp::new(13554361889612894269),
        Fp::new(12653635137853489769),
        Fp::new(4208693237683046800),
        Fp::new(5048957914766516287),
        Fp::new(3603232588024916254),
        Fp::new(8043443424472390486),
        Fp::new(2793772692762899036),
        Fp::new(2812395623847006123),
        Fp::new(17814903877314473260),
        Fp::new(4617387519742390220),
        Fp::new(17271454912209369477),
    ],
    [
        Fp::new(501154067648155634),
        Fp::new(776413635959717188),
        Fp::new(7246212854545673760),
        Fp::new(17574855500146756211),
        Fp::new(7986525476453271668),
        Fp::new(13649171140377007879),
        Fp::new(7797254645544706695),
        Fp::new(18433244849619790138),
        Fp::new(1234303719231864833),
        Fp::new(14298385539756159739),
        Fp::new(14403164449988626923),
    ],
    [
        Fp::new(3380811521665421538),
        Fp::new(8727321330908146361),
        Fp::new(3955409799579559961),
        Fp::new(3339333888336303924),
        Fp::new(759323967050526397),
        Fp::new(5205294308108865305),
        Fp::new(5900300433153818643),
        Fp::new(13851869672402414027),
        Fp::new(5386682578087687601),
        Fp::new(12086883693046019951),
        Fp::new(6423748754294986832),
    ],
    [
        Fp::new(6166063930501477577),
        Fp::new(13705550032600702069),
        Fp::new(719934703848329223),
        Fp::new(1077623097608858594),
        Fp::new(14564835855045911933),
        Fp::new(13525602396347738632),
        Fp::new(16400857888223360782),
        Fp::new(16714146392101839677),
        Fp::new(6910924023381260731),
        Fp::new(184041681713360185),
        Fp::new(1715564452880218307),
    ],
    [
        Fp::new(15279455465556128221),
        Fp::new(2280314165251043834),
        Fp::new(20408913243633763),
        Fp::new(15242886575581237147),
        Fp::new(72284874393716348),
        Fp::new(12509080551577613041),
        Fp::new(516172342737654039),
        Fp::new(7754926961164928811),
        Fp::new(15625746105591169890),
        Fp::new(7509898374806321037),
        Fp::new(12420852404554686392),
    ],
    [
        Fp::new(11038503205720364921),
        Fp::new(14844915958505263777),
        Fp::new(6865108060410275137),
        Fp::new(9371749032461006049),
        Fp::new(5443932381017742197),
        Fp::new(8007797513539847331),
        Fp::new(12562284025768602160),
        Fp::new(1787272432393771459),
        Fp::new(3823083059709586122),
        Fp::new(7473177234877659716),
        Fp::new(11287353350233770681),
    ],
    [
        Fp::new(6660350873769463660),
        Fp::new(6282689462668041393),
        Fp::new(3035644164662993689),
        Fp::new(3340647035578798023),
        Fp::new(10442111131997969928),
        Fp::new(16396862477987587126),
        Fp::new(10681243678428642138),
        Fp::new(14616530772654872360),
        Fp::new(12480949792847319034),
        Fp::new(8413217378708148047),
        Fp::new(7723266752541471048),
    ],
    [
        Fp::new(6912706702384783280),
        Fp::new(4424306136217508506),
        Fp::new(1803381219602908141),
        Fp::new(2660010183792201644),
        Fp::new(2606947164709508210),
        Fp::new(6822158739120852807),
        Fp::new(10732006137166454252),
        Fp::new(17053668181508829905),
        Fp::new(9789728551021227611),
        Fp::new(11305393825637579091),
        Fp::new(4738316080822123878),
    ],
    [
        Fp::new(1156555263915682148),
        Fp::new(6884817270917034867),
        Fp::new(2840238617231832669),
        Fp::new(6887637099744311427),
        Fp::new(16382821351275608019),
        Fp::new(13392830502066914253),
        Fp::new(7069105398725090623),
        Fp::new(18075455748201948790),
        Fp::new(15753004024578308038),
        Fp::new(15532599079399910129),
        Fp::new(16631295750232733355),
    ],
];

pub(crate) const M_00: Fp = Fp::new(8062354434719430557);
