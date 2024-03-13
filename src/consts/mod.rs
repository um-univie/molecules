use crate::molecular_formula::Isotope;
use lazy_static::lazy_static;
use phf::phf_map;

// This speeds up by checking only within a certain range
pub const BOND_SEARCH_THRESHOLD: f64 = 2.0;
// This was chosen to be the same as the threshold used in QCxMS
pub const BOND_TOLERANCE: f64 = 1.5;


lazy_static! {
    pub static ref OXIDATION_STATES: [Option<Vec<i8>>;119] = [
    None, // Padding
    Some(vec![-1, 0, 1]), // H
    Some(vec![0]), // HE
    Some(vec![0, 1]), // LI
    Some(vec![0, 1, 2]), // BE
    Some(vec![-5, -1, 0, 1, 2, 3]), // B
    Some(vec![-4, -3, -2, -1, 0, 1, 2, 3, 4]), // C
    Some(vec![-3, -2, -1, 0, 1, 2, 3, 4, 5]), // N
    Some(vec![-2, -1, 0, 1, 2]), // O
    Some(vec![-1, 0]), // F
    Some(vec![0]), // NE
    Some(vec![-1, 0, 1]), // NA
    Some(vec![0, 1, 2]), // MG
    Some(vec![-2, -1, 0, 1, 2, 3]), // AL
    Some(vec![-4, -3, -2, -1, 0, 1, 2, 3, 4]), // SI
    Some(vec![-3, -2, -1, 0, 1, 2, 3, 4, 5]), // P
    Some(vec![-2, -1, 0, 1, 2, 3, 4, 5, 6]), // S
    Some(vec![-1, 0, 1, 2, 3, 4, 5, 6, 7]), // CL
    Some(vec![0]), // AR
    Some(vec![-1, 1]), // K
    Some(vec![1, 2]), // CA
    Some(vec![0, 1, 2, 3]), // SC
    Some(vec![-2, -1, 0, 1, 2, 3, 4]), // TI
    Some(vec![-3, -1, 0, 1, 2, 3, 4, 5]), // V
    Some(vec![-4, -2, -1, 0, 1, 2, 3, 4, 5, 6]), // CR
    Some(vec![-3, -1, 0, 1, 2, 3, 4, 5, 6, 7]), // MN
    Some(vec![-4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]), // FE
    Some(vec![-3, -1, 0, 1, 2, 3, 4, 5]), // CO
    Some(vec![-2, -1, 0, 1, 2, 3, 4]), // NI
    Some(vec![-2, 0, 1, 2, 3, 4]), // CU
    Some(vec![-2, 0, 1, 2]), // ZN
    Some(vec![-5, -4, -3, -2, -1, 0, 1, 2, 3]), // GA
    Some(vec![-4, -3, -2, -1, 0, 1, 2, 3, 4]), // GE
    Some(vec![-3, -2, -1, 0, 1, 2, 3, 4, 5]), // AS
    Some(vec![-2, -1, 0, 1, 2, 3, 4, 5, 6]), // SE
    Some(vec![-1, 0, 1, 2, 3, 4, 5, 7]), // BR
    Some(vec![0, 1, 2]), // KR
    Some(vec![-1, 1]), // RB
    Some(vec![1, 2]), // SR
    Some(vec![0, 1, 2, 3]), // Y
    Some(vec![-2, 0, 1, 2, 3, 4]), // ZR
    Some(vec![-3, -1, 0, 1, 2, 3, 4, 5]), // NB
    Some(vec![-4, -2, -1, 0, 1, 2, 3, 4, 5, 6]), // MO
    Some(vec![-1, 0, 1, 2, 3, 4, 5, 6, 7]), // TC
    Some(vec![-4, -2, 0, 1, 2, 3, 4, 5, 6, 7, 8]), // RU
    Some(vec![-3, -1, 0, 1, 2, 3, 4, 5, 6, 7]), // RH
    Some(vec![0, 1, 2, 3, 4, 5]), // PD
    Some(vec![-2, -1, 0, 1, 2, 3]), // AG
    Some(vec![-2, 1, 2]), // CD
    Some(vec![-5, -2, -1, 0, 1, 2, 3]), // IN
    Some(vec![-4, -3, -2, -1, 0, 1, 2, 3, 4]), // SN
    Some(vec![-3, -2, -1, 0, 1, 2, 3, 4, 5]), // SB
    Some(vec![-2, -1, 0, 1, 2, 3, 4, 5, 6]), // TE
    Some(vec![-1, 0, 1, 2, 3, 4, 5, 6, 7]), // I
    Some(vec![0, 2, 4, 6, 8]), // XE
    Some(vec![-1, 1]), // CS
    Some(vec![1, 2]), // BA
    Some(vec![0, 1, 2, 3]), // LA
    Some(vec![2, 3, 4]), // CE
    Some(vec![0, 1, 2, 3, 4, 5]), // PR
    Some(vec![0, 2, 3, 4]), // ND
    Some(vec![2, 3]), // PM
    Some(vec![0, 1, 2, 3]), // SM
    Some(vec![0, 2, 3]), // EU
    Some(vec![0, 1, 2, 3]), // GD
    Some(vec![0, 1, 2, 3, 4]), // TB
    Some(vec![0, 2, 3, 4]), // DY
    Some(vec![0, 2, 3]), // HO
    Some(vec![0, 2, 3]), // ER
    Some(vec![0, 1, 2, 3]), // TM
    Some(vec![0, 1, 2, 3]), // YB
    Some(vec![0, 2, 3]), // LU
    Some(vec![-2, 0, 1, 2, 3, 4]), // HF
    Some(vec![-3, -1, 0, 1, 2, 3, 4, 5]), // TA
    Some(vec![-4, -2, -1, 0, 1, 2, 3, 4, 5, 6]), // W
    Some(vec![-3, -1, 0, 1, 2, 3, 4, 5, 6, 7]), // RE
    Some(vec![-4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8]), // OS
    Some(vec![-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]), // IR
    Some(vec![-3, -2, -1, 0, 1, 2, 3, 4, 5, 6]), // PT
    Some(vec![-3, -2, -1, 0, 1, 2, 3, 5]), // AU
    Some(vec![-2, 1, 2]), // HG
    Some(vec![-5, -2, -1, 1, 2, 3]), // TL
    Some(vec![-4, -2, -1, 0, 1, 2, 3, 4]), // PB
    Some(vec![-3, -2, -1, 0, 1, 2, 3, 4, 5]), // BI
    Some(vec![-2, 2, 4, 5, 6]), // PO
    Some(vec![-1, 1, 3, 5, 7]), // AT
    Some(vec![2, 6]), // RN
    Some(vec![1]), // FR
    Some(vec![2]), // RA
    Some(vec![3]), // AC
    Some(vec![-1, 1, 2, 3, 4]), // TH
    Some(vec![2, 3, 4, 5]), // PA
    Some(vec![-1, 1, 2, 3, 4, 5, 6]), // U
    Some(vec![2, 3, 4, 5, 6, 7]), // NP
    Some(vec![2, 3, 4, 5, 6, 7, 8]), // PU
    Some(vec![2, 3, 4, 5, 6, 7]), // AM
    Some(vec![3, 4, 5, 6]), // CM
    Some(vec![2, 3, 4, 5]), // BK
    Some(vec![2, 3, 4, 5]), // CF
    Some(vec![2, 3, 4]), // ES
    Some(vec![2, 3]), // FM
    Some(vec![2, 3]), // MD
    Some(vec![2, 3]), // NO
    Some(vec![3]), // LR
    Some(vec![4]), // RF
    Some(vec![5]), // DB
    Some(vec![0, 6]), // SG
    Some(vec![7]), // BH
    Some(vec![8]), // HS
    None, // MT
    None, // DS
    None, // RG
    Some(vec![2]), // CN
    None, // NH
    None, // FL
    None, // MC
    None, // LV
    None, // TS
    None, // OG
    ];
}
pub const PRIMES: [u16; 1000] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
    809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929,
    937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039,
    1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153,
    1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279,
    1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409,
    1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
    1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613,
    1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741,
    1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873,
    1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999,
    2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113,
    2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251,
    2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371,
    2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477,
    2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647,
    2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731,
    2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857,
    2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001,
    3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163,
    3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299,
    3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407,
    3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539,
    3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659,
    3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793,
    3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919,
    3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051,
    4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201,
    4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327,
    4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463,
    4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603,
    4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733,
    4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903,
    4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009,
    5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153,
    5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303,
    5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441,
    5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569,
    5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701,
    5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843,
    5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987,
    6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131,
    6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269,
    6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373,
    6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553,
    6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691,
    6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829,
    6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967,
    6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109,
    7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247,
    7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451,
    7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559,
    7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687,
    7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841,
    7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
];

// Array based lookup with dummy value at position 0, based on the NIST dataset: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&isotype=some
// Only up to Xenon
pub const STANDARD_ATOMIC_WEIGHTS: [f64; 55] = [
    0.0, 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990,
    24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942,
    51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904,
    83.798, 85.4678, 87.62, 88.90584, 91.224, 92.90637, 95.95, 98.0, 101.07, 102.90549, 106.42,
    107.8682, 112.414, 114.818, 118.710, 121.760, 127.60, 126.90447, 131.293,
];

// Only single bonds are considered, and the values are taken from: https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.200800987
//pub static ATOMIC_RADII: [f64; 55] = [
//    0.0, 0.32, 0.46, 1.33, 1.02, 0.85, 0.75, 0.75, 0.73, 0.72, 0.58, 1.60, 1.39, 1.26, 1.16, 1.11,
//    1.03, 0.99, 0.97, 2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16, 1.15, 1.17, 1.25, 1.26,
//    1.22, 1.21, 1.16, 1.14, 1.12, 2.16, 1.91, 1.62, 1.45, 1.34, 1.29, 1.27, 1.25, 1.25, 1.20, 1.39,
//    1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40, // XE
//];
//
// These values are adapted from QCxMS and are based on the values from: https://pubmed.ncbi.nlm.nih.gov/19058281/ for all elements after Plutonium, the radii are reduced by 10% for metals
pub const ATOMIC_RADII: [f64; 119] = [
    0.0, 0.32, 0.37, 1.30, 0.99, 0.84, 0.75, 0.71, 0.64, 0.60, 0.62, 1.60, 1.40, 1.24, 1.14, 1.09,
    1.04, 1.00, 1.01, 2.00, 1.74, 1.59, 1.48, 1.44, 1.30, 1.29, 1.24, 1.18, 1.17, 1.22, 1.20, 1.23,
    1.20, 1.20, 1.18, 1.17, 1.16, 2.15, 1.90, 1.76, 1.64, 1.56, 1.46, 1.38, 1.36, 1.34, 1.30, 1.36,
    1.40, 1.42, 1.40, 1.40, 1.37, 1.36, 1.36, 2.38, 2.06, 1.94, 1.84, 1.90, 1.88, 1.86, 1.85, 1.83,
    1.82, 1.81, 1.80, 1.79, 1.77, 1.77, 1.78, 1.74, 1.64, 1.58, 1.50, 1.41, 1.36, 1.32, 1.30, 1.30,
    1.32, 1.44, 1.45, 1.50, 1.42, 1.48, 1.46, 2.42, 2.11, 2.01, 1.90, 1.84, 1.83, 1.80, 1.80, 1.49,
    1.49, 1.51, 1.51, 1.48, 1.50, 1.56, 1.58, 1.45, 1.41, 1.34, 1.29, 1.27, 1.21, 1.16, 1.15, 1.09,
    1.22, 1.36, 1.43, 1.46, 1.58, 1.48, 1.57,
];
// This is a heuristic based on experience, and is not (officially) based on any scientific data yet
// We will add an initial check soon that will try to find the initial valence of each atom
pub const VALENCIES: [i8; 33] = [
    0, // Dummy value
    1, // H
    0, // HE
    1, // LI
    2, // ...
    3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4,
];

// Taken from: https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
// The first element is not used, but is included for completeness
pub const MONOISOTOPIC_MASSES: [f64; 119] = [
    0.00000000000,
    1.00782503223,
    3.0160293201,
    6.0151228874,
    9.012183065,
    10.01293695,
    12.0000000,
    14.00307400443,
    15.99491461957,
    18.99840316273,
    19.9924401762,
    22.989769282,
    23.985041697,
    26.98153853,
    27.97692653465,
    30.97376199842,
    31.9720711744,
    34.968852682,
    35.967545105,
    38.9637064864,
    39.962590863,
    44.95590828,
    45.95262772,
    49.94715601,
    49.94604183,
    54.93804391,
    53.93960899,
    58.93319429,
    57.93534241,
    62.92959772,
    63.92914201,
    68.9255735,
    69.92424875,
    74.92159457,
    73.922475934,
    78.9183376,
    77.92036494,
    84.9117897379,
    83.9134191,
    88.9058403,
    89.9046977,
    92.906373,
    91.90680796,
    96.9063667,
    95.90759025,
    102.905498,
    101.9056022,
    106.9050916,
    105.9064599,
    112.90406184,
    111.90482387,
    120.903812,
    119.9040593,
    126.9044719,
    123.905892,
    132.905451961,
    129.9063207,
    137.9071149,
    135.90712921,
    140.9076576,
    141.907729,
    144.9127559,
    143.9120065,
    150.9198578,
    151.9197995,
    158.9253547,
    155.9242847,
    164.9303288,
    161.9287884,
    168.9342179,
    167.9338896,
    174.9407752,
    173.9400461,
    179.9474648,
    179.9467108,
    184.9529545,
    183.9524885,
    190.9605893,
    189.9599297,
    196.96656879,
    195.9658326,
    202.9723446,
    203.973044,
    208.9803991,
    208.9824308,
    209.9871479,
    210.9906011,
    223.019736,
    223.0185023,
    227.0277523,
    230.0331341,
    231.0358842,
    233.0396355,
    236.04657,
    238.0495601,
    241.0568293,
    243.0613893,
    247.0703073,
    249.0748539,
    252.08298,
    257.0951061,
    258.0984315,
    259.10103,
    262.10961,
    267.12179,
    268.12567,
    271.13393,
    272.13826,
    270.13429,
    276.15159,
    281.16451,
    280.16514,
    285.17712,
    284.17873,
    289.19042,
    288.19274,
    293.20449,
    292.20746,
    294.21392,
];

pub const ATOMIC_SYMBOLS: [&str; 119] = [
    "None", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

pub const ATOMIC_NUMBERS: phf::Map<&'static str, u8> = phf_map! {
"H" => 1,
"HE" => 2,
"LI" => 3,
"BE" => 4,
"B" => 5,
"C" => 6,
"N" => 7,
"O" => 8,
"F" => 9,
"NE" => 10,
"NA" => 11,
"MG" => 12,
"AL" => 13,
"SI" => 14,
"P" => 15,
"S" => 16,
"CL" => 17,
"AR" => 18,
"K" => 19,
"CA" => 20,
"SC" => 21,
"TI" => 22,
"V" => 23,
"CR" => 24,
"MN" => 25,
"FE" => 26,
"NI" => 27,
"CO" => 28,
"CU" => 29,
"ZN" => 30,
"GA" => 31,
"GE" => 32,
"AS" => 33,
"SE" => 34,
"BR" => 35,
"KR" => 36,
"RB" => 37,
"SR" => 38,
"Y" => 39,
"ZR" => 40,
"NB" => 41,
"MO" => 42,
"TC" => 43,
"RU" => 44,
"RH" => 45,
"PD" => 46,
"AG" => 47,
"CD" => 48,
"IN" => 49,
"SN" => 50,
"SB" => 51,
"TE" => 52,
"I" => 53,
"XE" => 54,
"CS" => 55,
"BA" => 56,
"LA" => 57,
"CE" => 58,
"PR" => 59,
"ND" => 60,
"PM" => 61,
"SM" => 62,
"EU" => 63,
"GD" => 64,
"TB" => 65,
"DY" => 66,
"HO" => 67,
"ER" => 68,
"TM" => 69,
"YB" => 70,
"LU" => 71,
"HF" => 72,
"TA" => 73,
"W" => 74,
"RE" => 75,
"OS" => 76,
"IR" => 77,
"PT" => 78,
"AU" => 79,
"HG" => 80,
"TL" => 81,
"PB" => 82,
"BI" => 83,
"TH" => 90,
"PA" => 91,
"U" => 92,
"NP" => 93,
"PU" => 94,
"AM" => 95,
"CM" => 96,
"BK" => 97,
"CF" => 98,
"ES" => 99,
"FM" => 100,
"MD" => 101,
"NO" => 102,
"LR" => 103,
"RF" => 104,
"DB" => 105,
"SG" => 106,
"BH" => 107,
"HS" => 108,
"MT" => 109,
"DS" => 110,
"RG" => 111,
"CN" => 112,
"NH" => 113,
"FL" => 114,
"MC" => 115,
"LV" => 116,
"TS" => 117,
"OG" => 118,
};

// Electronegativities based on the Pauling scale, scaled by a factor of 100 to avoid floating point errors
pub const ELECTRONEGATIVITIES: [u32; 103] = [
    0, 220, 0, 98, 157, 204, 255, 304, 350, 398, 0, 93, 131, 161, 190, 219, 258, 316, 0, 82, 100,
    136, 154, 163, 166, 155, 183, 188, 191, 190, 165, 181, 201, 218, 255, 296, 300, 82, 95, 122,
    133, 160, 216, 190, 220, 228, 220, 193, 169, 178, 196, 205, 210, 266, 260, 79, 89, 110, 112,
    113, 114, 0, 117, 0, 120, 0, 122, 123, 124, 125, 0, 127, 130, 150, 236, 190, 220, 220, 228,
    254, 200, 162, 233, 202, 200, 220, 0, 0, 90, 110, 130, 150, 138, 136, 128, 130, 130, 130, 130,
    130, 130, 130, 130,
];

// Taken from: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii&isotype=some
// Zero padding is used to make the array index match the atomic number
pub const ISOTOPES: [[Option<Isotope>; 4]; 119] = [
    [None, None, None, None],
    [
        Some(Isotope {
            mass: 1.00782503223,
            abundance: 0.999885,
        }),
        Some(Isotope {
            mass: 2.01410177812,
            abundance: 0.000115,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 4.00260325413,
            abundance: 0.99999866,
        }),
        Some(Isotope {
            mass: 3.0160293201,
            abundance: 1.34e-06,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 7.0160034366,
            abundance: 0.9241,
        }),
        Some(Isotope {
            mass: 6.0151228874,
            abundance: 0.0759,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 9.012183065,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 11.00930536,
            abundance: 0.801,
        }),
        Some(Isotope {
            mass: 10.01293695,
            abundance: 0.199,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 12.0,
            abundance: 0.9893,
        }),
        Some(Isotope {
            mass: 13.00335483507,
            abundance: 0.0107,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 14.00307400443,
            abundance: 0.99636,
        }),
        Some(Isotope {
            mass: 15.00010889888,
            abundance: 0.00364,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 15.99491461957,
            abundance: 0.99757,
        }),
        Some(Isotope {
            mass: 17.99915961286,
            abundance: 0.00205,
        }),
        Some(Isotope {
            mass: 16.9991317565,
            abundance: 0.00038,
        }),
        None,
    ],
    [
        Some(Isotope {
            mass: 18.99840316273,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 19.9924401762,
            abundance: 0.9048,
        }),
        Some(Isotope {
            mass: 21.991385114,
            abundance: 0.0925,
        }),
        Some(Isotope {
            mass: 20.993846685,
            abundance: 0.0027,
        }),
        None,
    ],
    [
        Some(Isotope {
            mass: 22.989769282,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 23.985041697,
            abundance: 0.7899,
        }),
        Some(Isotope {
            mass: 25.982592968,
            abundance: 0.1101,
        }),
        Some(Isotope {
            mass: 24.985836976,
            abundance: 0.1,
        }),
        None,
    ],
    [
        Some(Isotope {
            mass: 26.98153853,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 27.97692653465,
            abundance: 0.92223,
        }),
        Some(Isotope {
            mass: 28.9764946649,
            abundance: 0.04685,
        }),
        Some(Isotope {
            mass: 29.973770136,
            abundance: 0.03092,
        }),
        None,
    ],
    [
        Some(Isotope {
            mass: 30.97376199842,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 31.9720711744,
            abundance: 0.9499,
        }),
        Some(Isotope {
            mass: 33.967867004,
            abundance: 0.0425,
        }),
        Some(Isotope {
            mass: 32.9714589098,
            abundance: 0.0075,
        }),
        Some(Isotope {
            mass: 35.96708071,
            abundance: 0.0001,
        }),
    ],
    [
        Some(Isotope {
            mass: 34.968852682,
            abundance: 0.7576,
        }),
        Some(Isotope {
            mass: 36.965902602,
            abundance: 0.2424,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 39.9623831237,
            abundance: 0.996035,
        }),
        Some(Isotope {
            mass: 35.967545105,
            abundance: 0.003336,
        }),
        Some(Isotope {
            mass: 37.96273211,
            abundance: 0.000629,
        }),
        None,
    ],
    [
        Some(Isotope {
            mass: 38.9637064864,
            abundance: 0.932581,
        }),
        Some(Isotope {
            mass: 40.9618252579,
            abundance: 0.067302,
        }),
        Some(Isotope {
            mass: 39.963998166,
            abundance: 0.000117,
        }),
        None,
    ],
    [
        Some(Isotope {
            mass: 39.962590863,
            abundance: 0.96941,
        }),
        Some(Isotope {
            mass: 43.95548156,
            abundance: 0.02086,
        }),
        Some(Isotope {
            mass: 41.95861783,
            abundance: 0.00647,
        }),
        Some(Isotope {
            mass: 47.95252276,
            abundance: 0.00187,
        }),
    ],
    [
        Some(Isotope {
            mass: 44.95590828,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 47.94794198,
            abundance: 0.7372,
        }),
        Some(Isotope {
            mass: 45.95262772,
            abundance: 0.0825,
        }),
        Some(Isotope {
            mass: 46.95175879,
            abundance: 0.0744,
        }),
        Some(Isotope {
            mass: 48.94786568,
            abundance: 0.0541,
        }),
    ],
    [
        Some(Isotope {
            mass: 50.94395704,
            abundance: 0.9975,
        }),
        Some(Isotope {
            mass: 49.94715601,
            abundance: 0.0025,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 51.94050623,
            abundance: 0.83789,
        }),
        Some(Isotope {
            mass: 52.94064815,
            abundance: 0.09501,
        }),
        Some(Isotope {
            mass: 49.94604183,
            abundance: 0.04345,
        }),
        Some(Isotope {
            mass: 53.93887916,
            abundance: 0.02365,
        }),
    ],
    [
        Some(Isotope {
            mass: 54.93804391,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 55.93493633,
            abundance: 0.91754,
        }),
        Some(Isotope {
            mass: 53.93960899,
            abundance: 0.05845,
        }),
        Some(Isotope {
            mass: 56.93539284,
            abundance: 0.02119,
        }),
        Some(Isotope {
            mass: 57.93327443,
            abundance: 0.00282,
        }),
    ],
    [
        Some(Isotope {
            mass: 58.93319429,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 57.93534241,
            abundance: 0.68077,
        }),
        Some(Isotope {
            mass: 59.93078588,
            abundance: 0.26223,
        }),
        Some(Isotope {
            mass: 61.92834537,
            abundance: 0.036346,
        }),
        Some(Isotope {
            mass: 60.93105557,
            abundance: 0.011399,
        }),
    ],
    [
        Some(Isotope {
            mass: 62.92959772,
            abundance: 0.6915,
        }),
        Some(Isotope {
            mass: 64.9277897,
            abundance: 0.3085,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 63.92914201,
            abundance: 0.4917,
        }),
        Some(Isotope {
            mass: 65.92603381,
            abundance: 0.2773,
        }),
        Some(Isotope {
            mass: 67.92484455,
            abundance: 0.1845,
        }),
        Some(Isotope {
            mass: 66.92712775,
            abundance: 0.0404,
        }),
    ],
    [
        Some(Isotope {
            mass: 68.9255735,
            abundance: 0.60108,
        }),
        Some(Isotope {
            mass: 70.92470258,
            abundance: 0.39892,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 73.921177761,
            abundance: 0.365,
        }),
        Some(Isotope {
            mass: 71.922075826,
            abundance: 0.2745,
        }),
        Some(Isotope {
            mass: 69.92424875,
            abundance: 0.2057,
        }),
        Some(Isotope {
            mass: 72.923458956,
            abundance: 0.0775,
        }),
    ],
    [
        Some(Isotope {
            mass: 74.92159457,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 79.9165218,
            abundance: 0.4961,
        }),
        Some(Isotope {
            mass: 77.91730928,
            abundance: 0.2377,
        }),
        Some(Isotope {
            mass: 75.919213704,
            abundance: 0.0937,
        }),
        Some(Isotope {
            mass: 81.9166995,
            abundance: 0.0873,
        }),
    ],
    [
        Some(Isotope {
            mass: 78.9183376,
            abundance: 0.5069,
        }),
        Some(Isotope {
            mass: 80.9162897,
            abundance: 0.4931,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 83.9114977282,
            abundance: 0.56987,
        }),
        Some(Isotope {
            mass: 85.9106106269,
            abundance: 0.17279,
        }),
        Some(Isotope {
            mass: 81.91348273,
            abundance: 0.11593,
        }),
        Some(Isotope {
            mass: 82.91412716,
            abundance: 0.115,
        }),
    ],
    [
        Some(Isotope {
            mass: 84.9117897379,
            abundance: 0.7217,
        }),
        Some(Isotope {
            mass: 86.909180531,
            abundance: 0.2783,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 87.9056125,
            abundance: 0.8258,
        }),
        Some(Isotope {
            mass: 85.9092606,
            abundance: 0.0986,
        }),
        Some(Isotope {
            mass: 86.9088775,
            abundance: 0.07,
        }),
        Some(Isotope {
            mass: 83.9134191,
            abundance: 0.0056,
        }),
    ],
    [
        Some(Isotope {
            mass: 88.9058403,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 89.9046977,
            abundance: 0.5145,
        }),
        Some(Isotope {
            mass: 93.9063108,
            abundance: 0.1738,
        }),
        Some(Isotope {
            mass: 91.9050347,
            abundance: 0.1715,
        }),
        Some(Isotope {
            mass: 90.9056396,
            abundance: 0.1122,
        }),
    ],
    [
        Some(Isotope {
            mass: 92.906373,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 97.90540482,
            abundance: 0.2439,
        }),
        Some(Isotope {
            mass: 95.90467612,
            abundance: 0.1667,
        }),
        Some(Isotope {
            mass: 94.90583877,
            abundance: 0.1584,
        }),
        Some(Isotope {
            mass: 91.90680796,
            abundance: 0.1453,
        }),
    ],
    [None, None, None, None],
    [
        Some(Isotope {
            mass: 101.9043441,
            abundance: 0.3155,
        }),
        Some(Isotope {
            mass: 103.9054275,
            abundance: 0.1862,
        }),
        Some(Isotope {
            mass: 100.9055769,
            abundance: 0.1706,
        }),
        Some(Isotope {
            mass: 98.9059341,
            abundance: 0.1276,
        }),
    ],
    [
        Some(Isotope {
            mass: 102.905498,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 105.9034804,
            abundance: 0.2733,
        }),
        Some(Isotope {
            mass: 107.9038916,
            abundance: 0.2646,
        }),
        Some(Isotope {
            mass: 104.9050796,
            abundance: 0.2233,
        }),
        Some(Isotope {
            mass: 109.9051722,
            abundance: 0.1172,
        }),
    ],
    [
        Some(Isotope {
            mass: 106.9050916,
            abundance: 0.51839,
        }),
        Some(Isotope {
            mass: 108.9047553,
            abundance: 0.48161,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 113.90336509,
            abundance: 0.2873,
        }),
        Some(Isotope {
            mass: 111.90276287,
            abundance: 0.2413,
        }),
        Some(Isotope {
            mass: 110.90418287,
            abundance: 0.128,
        }),
        Some(Isotope {
            mass: 109.90300661,
            abundance: 0.1249,
        }),
    ],
    [
        Some(Isotope {
            mass: 114.903878776,
            abundance: 0.9571,
        }),
        Some(Isotope {
            mass: 112.90406184,
            abundance: 0.0429,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 119.90220163,
            abundance: 0.3258,
        }),
        Some(Isotope {
            mass: 117.90160657,
            abundance: 0.2422,
        }),
        Some(Isotope {
            mass: 115.9017428,
            abundance: 0.1454,
        }),
        Some(Isotope {
            mass: 118.90331117,
            abundance: 0.0859,
        }),
    ],
    [
        Some(Isotope {
            mass: 120.903812,
            abundance: 0.5721,
        }),
        Some(Isotope {
            mass: 122.9042132,
            abundance: 0.4279,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 129.906222748,
            abundance: 0.3408,
        }),
        Some(Isotope {
            mass: 127.90446128,
            abundance: 0.3174,
        }),
        Some(Isotope {
            mass: 125.9033109,
            abundance: 0.1884,
        }),
        Some(Isotope {
            mass: 124.9044299,
            abundance: 0.0707,
        }),
    ],
    [
        Some(Isotope {
            mass: 126.9044719,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 131.9041550856,
            abundance: 0.269086,
        }),
        Some(Isotope {
            mass: 128.9047808611,
            abundance: 0.264006,
        }),
        Some(Isotope {
            mass: 130.90508406,
            abundance: 0.212324,
        }),
        Some(Isotope {
            mass: 133.90539466,
            abundance: 0.104357,
        }),
    ],
    [
        Some(Isotope {
            mass: 132.905451961,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 137.905247,
            abundance: 0.71698,
        }),
        Some(Isotope {
            mass: 136.90582714,
            abundance: 0.11232,
        }),
        Some(Isotope {
            mass: 135.90457573,
            abundance: 0.07854,
        }),
        Some(Isotope {
            mass: 134.90568838,
            abundance: 0.06592,
        }),
    ],
    [
        Some(Isotope {
            mass: 138.9063563,
            abundance: 0.9991119,
        }),
        Some(Isotope {
            mass: 137.9071149,
            abundance: 0.0008881,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 139.9054431,
            abundance: 0.8845,
        }),
        Some(Isotope {
            mass: 141.9092504,
            abundance: 0.11114,
        }),
        Some(Isotope {
            mass: 137.905991,
            abundance: 0.00251,
        }),
        Some(Isotope {
            mass: 135.90712921,
            abundance: 0.00185,
        }),
    ],
    [
        Some(Isotope {
            mass: 140.9076576,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 141.907729,
            abundance: 0.27152,
        }),
        Some(Isotope {
            mass: 143.910093,
            abundance: 0.23798,
        }),
        Some(Isotope {
            mass: 145.9131226,
            abundance: 0.17189,
        }),
        Some(Isotope {
            mass: 142.90982,
            abundance: 0.12174,
        }),
    ],
    [None, None, None, None],
    [
        Some(Isotope {
            mass: 151.9197397,
            abundance: 0.2675,
        }),
        Some(Isotope {
            mass: 153.9222169,
            abundance: 0.2275,
        }),
        Some(Isotope {
            mass: 146.9149044,
            abundance: 0.1499,
        }),
        Some(Isotope {
            mass: 148.9171921,
            abundance: 0.1382,
        }),
    ],
    [
        Some(Isotope {
            mass: 152.921238,
            abundance: 0.5219,
        }),
        Some(Isotope {
            mass: 150.9198578,
            abundance: 0.4781,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 157.9241123,
            abundance: 0.2484,
        }),
        Some(Isotope {
            mass: 159.9270624,
            abundance: 0.2186,
        }),
        Some(Isotope {
            mass: 155.9221312,
            abundance: 0.2047,
        }),
        Some(Isotope {
            mass: 156.9239686,
            abundance: 0.1565,
        }),
    ],
    [
        Some(Isotope {
            mass: 158.9253547,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 163.9291819,
            abundance: 0.2826,
        }),
        Some(Isotope {
            mass: 161.9268056,
            abundance: 0.25475,
        }),
        Some(Isotope {
            mass: 162.9287383,
            abundance: 0.24896,
        }),
        Some(Isotope {
            mass: 160.9269405,
            abundance: 0.18889,
        }),
    ],
    [
        Some(Isotope {
            mass: 164.9303288,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 165.9302995,
            abundance: 0.33503,
        }),
        Some(Isotope {
            mass: 167.9323767,
            abundance: 0.26978,
        }),
        Some(Isotope {
            mass: 166.9320546,
            abundance: 0.22869,
        }),
        Some(Isotope {
            mass: 169.9354702,
            abundance: 0.1491,
        }),
    ],
    [
        Some(Isotope {
            mass: 168.9342179,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 173.9388664,
            abundance: 0.32026,
        }),
        Some(Isotope {
            mass: 171.9363859,
            abundance: 0.2168,
        }),
        Some(Isotope {
            mass: 172.9382151,
            abundance: 0.16103,
        }),
        Some(Isotope {
            mass: 170.9363302,
            abundance: 0.1409,
        }),
    ],
    [
        Some(Isotope {
            mass: 174.9407752,
            abundance: 0.97401,
        }),
        Some(Isotope {
            mass: 175.9426897,
            abundance: 0.02599,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 179.946557,
            abundance: 0.3508,
        }),
        Some(Isotope {
            mass: 177.9437058,
            abundance: 0.2728,
        }),
        Some(Isotope {
            mass: 176.9432277,
            abundance: 0.186,
        }),
        Some(Isotope {
            mass: 178.9458232,
            abundance: 0.1362,
        }),
    ],
    [
        Some(Isotope {
            mass: 180.9479958,
            abundance: 0.9998799,
        }),
        Some(Isotope {
            mass: 179.9474648,
            abundance: 0.0001201,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 183.95093092,
            abundance: 0.3064,
        }),
        Some(Isotope {
            mass: 185.9543628,
            abundance: 0.2843,
        }),
        Some(Isotope {
            mass: 181.94820394,
            abundance: 0.265,
        }),
        Some(Isotope {
            mass: 182.95022275,
            abundance: 0.1431,
        }),
    ],
    [
        Some(Isotope {
            mass: 186.9557501,
            abundance: 0.626,
        }),
        Some(Isotope {
            mass: 184.9529545,
            abundance: 0.374,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 191.961477,
            abundance: 0.4078,
        }),
        Some(Isotope {
            mass: 189.9584437,
            abundance: 0.2626,
        }),
        Some(Isotope {
            mass: 188.9581442,
            abundance: 0.1615,
        }),
        Some(Isotope {
            mass: 187.9558352,
            abundance: 0.1324,
        }),
    ],
    [
        Some(Isotope {
            mass: 192.9629216,
            abundance: 0.627,
        }),
        Some(Isotope {
            mass: 190.9605893,
            abundance: 0.373,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 194.9647917,
            abundance: 0.3378,
        }),
        Some(Isotope {
            mass: 193.9626809,
            abundance: 0.3286,
        }),
        Some(Isotope {
            mass: 195.96495209,
            abundance: 0.2521,
        }),
        Some(Isotope {
            mass: 197.9678949,
            abundance: 0.07356,
        }),
    ],
    [
        Some(Isotope {
            mass: 196.96656879,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 201.9706434,
            abundance: 0.2986,
        }),
        Some(Isotope {
            mass: 199.96832659,
            abundance: 0.231,
        }),
        Some(Isotope {
            mass: 198.96828064,
            abundance: 0.1687,
        }),
        Some(Isotope {
            mass: 200.97030284,
            abundance: 0.1318,
        }),
    ],
    [
        Some(Isotope {
            mass: 204.9744278,
            abundance: 0.7048,
        }),
        Some(Isotope {
            mass: 202.9723446,
            abundance: 0.2952,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 207.9766525,
            abundance: 0.524,
        }),
        Some(Isotope {
            mass: 205.9744657,
            abundance: 0.241,
        }),
        Some(Isotope {
            mass: 206.9758973,
            abundance: 0.221,
        }),
        Some(Isotope {
            mass: 203.973044,
            abundance: 0.014,
        }),
    ],
    [
        Some(Isotope {
            mass: 208.9803991,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [
        Some(Isotope {
            mass: 230.0331341,
            abundance: 232.0377,
        }),
        Some(Isotope {
            mass: 232.0380558,
            abundance: 1.0,
        }),
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 231.0358842,
            abundance: 1.0,
        }),
        None,
        None,
        None,
    ],
    [
        Some(Isotope {
            mass: 233.0396355,
            abundance: 238.02891,
        }),
        Some(Isotope {
            mass: 238.0507884,
            abundance: 0.992742,
        }),
        Some(Isotope {
            mass: 235.0439301,
            abundance: 0.007204,
        }),
        Some(Isotope {
            mass: 234.0409523,
            abundance: 5.4e-05,
        }),
    ],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
    [None, None, None, None],
];
