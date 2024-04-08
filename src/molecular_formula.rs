use crate::consts::{
    ATOMIC_NUMBERS, ATOMIC_SYMBOLS, ISOTOPES, MONOISOTOPIC_MASSES, STANDARD_ATOMIC_WEIGHTS,
};
use core::{
    fmt::{Display, Formatter},
    str::FromStr,
};
use ndarray::{s, Array1};

use nohash_hasher::IntMap;

// TODO: make this a parameter
const PRECISION: f64 = 1000.0;

#[derive(Clone, Debug, PartialEq)]
pub struct MolecularFormula {
    elements: IntMap<u8, usize>,
}

impl MolecularFormula {
    pub fn new() -> Self {
        MolecularFormula {
            elements: IntMap::default(),
        }
    }
    pub fn elements(&self) -> &IntMap<u8, usize> {
        &self.elements
    }
    /// Merges two molecular formulas.
    ///
    /// # Examples
    /// ```
    /// use molecules::molecular_formula::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// let mut water_methane = water.clone();
    /// water_methane.merge(&methane);
    /// assert_eq!(water_methane, "CH4H2O".parse::<MolecularFormula>().unwrap());
    /// assert_eq!(water, "H2O".parse::<MolecularFormula>().unwrap());
    /// assert_eq!(methane, "CH4".parse::<MolecularFormula>().unwrap());
    /// assert_eq!(water_methane, "CH4H2O".parse::<MolecularFormula>().unwrap());
    ///
    /// ```

    pub fn merge(&mut self, other: &MolecularFormula) {
        for (atom, count) in &other.elements {
            *self.elements.entry(*atom).or_insert(0) += count;
        }
    }
    /// Combines two molecular formulas.
    /// # Examples
    /// ```
    /// use molecules::molecular_formula::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// let water_methane = water.combine(&methane);
    /// assert_eq!(water_methane, "CH4H2O".parse::<MolecularFormula>().unwrap());
    /// ```
    pub fn combine(&self, other: &MolecularFormula) -> Self {
        let mut combined = MolecularFormula::new();
        for (atom, count) in &self.elements {
            *combined.elements.entry(*atom).or_insert(0) += count;
        }
        for (atom, count) in &other.elements {
            *combined.elements.entry(*atom).or_insert(0) += count;
        }
        combined
    }
    /// Calculates the monoisotopic mass of the molecular formula.
    ///
    /// # Examples
    /// ```
    /// use molecules::molecular_formula::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// println!("{:?}", water);
    /// assert_eq!(water.monoisotopic_mass(), 18.01056468403);
    /// assert_eq!(methane.monoisotopic_mass(), 16.03130012892);
    ///
    /// ```
    pub fn monoisotopic_mass(&self) -> f64 {
        self.elements.iter().fold(0.0, |acc, (atom, count)| {
            acc + MONOISOTOPIC_MASSES[*atom as usize] * *count as f64
        })
    }

    /// Calculates the molecular mass of the molecular formula.
    /// # Examples
    /// ```
    /// use molecules::molecular_formula::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// assert_eq!(water.molecular_mass(), 18.015);
    /// assert_eq!(methane.molecular_mass(), 16.043);
    pub fn molecular_mass(&self) -> f64 {
        self.elements.iter().fold(0.0, |acc, (atom, count)| {
            acc + STANDARD_ATOMIC_WEIGHTS[*atom as usize] * *count as f64
        })
    }
    //This is commented out due to a saftey issue with the fft crate

    // Calculates the isotopic pattern of the molecular formula.
    //
    // # Arguments
    // * `max_mass_difference` - The maximum mass difference between the monoisotopic peak and the
    // last peak in the isotopic pattern to be calculated.
    // # Examples
    // ```
    // use molecules::molecular_formula::MolecularFormula;
    // let water = "H2O".parse::<MolecularFormula>().unwrap();
    // let isotopic_pattern_h2o = water.isotopic_pattern(3,1);
    //
    // println!("{:?}", isotopic_pattern_h2o);
    // println!("{:?}", water);
    // assert!(isotopic_pattern_h2o.len() == 3);
    //
    // let methane = "CH4".parse::<MolecularFormula>().unwrap();
    // let isotopic_pattern_methane = methane.isotopic_pattern(3,1);
    // println!("{:?}", isotopic_pattern_methane);
    // assert!(isotopic_pattern_methane.len() == 3);
    // println!("{:?}", water);
    // ```
    //  pub fn isotopic_pattern(
    //      &self,
    //      max_mass_difference: usize,
    //      resolution: usize,
    //  ) -> Vec<(f64, f64)> {
    //      let min_mass: usize = self.monoisotopic_mass().round() as usize;
    //      let max_mass = min_mass + max_mass_difference;
    //      let mut total_distribution = Array1::<f64>::zeros(max_mass * resolution + 1);
    //      total_distribution[0] = 1.0;

    //      let mut planner = FftPlanner::new();
    //      let fft = planner.plan_fft_forward(total_distribution.len());
    //      let ifft = planner.plan_fft_inverse(total_distribution.len());

    //      let mut atom_distribution = Array1::<f64>::zeros(max_mass * resolution + 1);
    //      for (atomic_number, count) in self.elements.iter() {
    //          if let Some(atom) = ISOTOPES.get(*atomic_number as usize) {
    //              for _ in 0..*count {
    //                  atom_distribution.fill(0.0);
    //                  for isotope in atom.iter().flatten() {
    //                      atom_distribution[(isotope.mass * resolution as f64).round() as usize] =
    //                          isotope.abundance;
    //                  }

    //                  // Preparing input for FFT
    //                  let mut distribution_fft = total_distribution
    //                      .iter()
    //                      .map(|&v| Complex::new(v, 0.0))
    //                      .collect::<Vec<Complex<f64>>>();
    //                  let mut atom_distribution_fft = atom_distribution
    //                      .iter()
    //                      .map(|&v| Complex::new(v, 0.0))
    //                      .collect::<Vec<Complex<f64>>>();
    //                  // Apply FFT
    //                  fft.process(&mut distribution_fft);
    //                  fft.process(&mut atom_distribution_fft);

    //                  // Convolution in the frequency domain
    //                  let mut convoluted = distribution_fft
    //                      .iter()
    //                      .zip(atom_distribution_fft.iter())
    //                      .map(|(a, b)| a * b)
    //                      .collect::<Vec<Complex<f64>>>();

    //                  // Apply inverse FFT
    //                  ifft.process(&mut convoluted);
    //                  total_distribution =
    //                      Array1::from(convoluted.iter().map(|v| v.re).collect::<Vec<f64>>());
    //              }
    //          }
    //      }

    //      // Normalization
    //      let sum = total_distribution.sum();
    //      total_distribution
    //          .slice(s![min_mass * resolution..max_mass * resolution])
    //          .mapv(|v| v / sum)
    //          .into_iter()
    //          .zip(
    //              (min_mass * resolution..max_mass * resolution)
    //                  .map(|v| v as f64 / resolution as f64),
    //          )
    //          .collect::<Vec<(f64, f64)>>()
    //  }
    // pub fn isotopic_pattern(&self, max_mass_difference: usize, resolution: usize) -> Vec<f64> {
    //     let min_mass: usize = self.monoisotopic_mass().round() as usize;
    //     let max_mass = min_mass + max_mass_difference;
    //     let mut total_distribution = Array1::<f64>::zeros(max_mass * resolution + 1);
    //     total_distribution[0] = 1.0;
    //  for (atomic_number, count) in self.elements.iter() {
    //         if let Some(atom) = ISOTOPES.get(*atomic_number as usize) {
    //             for _ in 0..*count {
    //                 let mut atom_distribution = Array1::<f64>::zeros(max_mass * resolution + 1);
    //                 for isotope in atom {
    //                     let Some(isotope) = isotope else {
    //                         break
    //                     };
    //                     atom_distribution[(isotope.mass * resolution as f64).round() as usize] =
    //                         isotope.abundance;
    //                 }

    //                 let mut planner = FftPlanner::new();
    //                 let fft = planner.plan_fft_forward(total_distribution.len());
    //                 let ifft = planner.plan_fft_inverse(total_distribution.len());

    //                 let mut input1: Vec<Complex<f64>> = total_distribution
    //                     .iter()
    //                     .map(|&v| Complex::new(v, 0.0))
    //                     .collect();
    //                 let mut input2: Vec<Complex<f64>> = atom_distribution
    //                     .iter()
    //                     .map(|&v| Complex::new(v, 0.0))
    //                     .collect();

    //                 fft.process(&mut input1);
    //                 fft.process(&mut input2);

    //                 let mut output: Vec<Complex<f64>> = input1
    //                     .iter()
    //                     .zip(input2.iter())
    //                     .map(|(a, b)| a * b)
    //                     .collect();
    //                 ifft.process(&mut output);
    //                 total_distribution =
    //                     Array1::from(output.iter().map(|v| v.re).collect::<Vec<f64>>());
    //             }
    //         }
    //     }
    //     total_distribution
    //         .slice(s![min_mass * resolution..max_mass * resolution])
    //         .mapv(|v| v / total_distribution.sum())
    //         .to_vec()
    // }
}

impl Default for MolecularFormula {
    fn default() -> Self {
        Self::new()
    }
}

impl Display for MolecularFormula {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut elements = self
            .elements
            .iter()
            .map(|(&atom, &count)| (ATOMIC_SYMBOLS[atom as usize], count))
            .collect::<Vec<_>>();
        elements.sort_by_key(|&(atom, _)| atom);
        for (atom, count) in elements {
            write!(f, "{}{}", atom, count)?;
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct ParseFormulaError;

impl std::fmt::Display for ParseFormulaError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "could not parse the provided string as a molecular formula"
        )
    }
}

impl std::error::Error for ParseFormulaError {}

impl FromStr for MolecularFormula {
    type Err = ParseFormulaError;

    /// Parses a molecular formula from a string.
    ///
    /// # Panics
    ///
    /// This function does not panic. But it will ignore everything that is neither a number nor a
    /// letter
    ///
    /// # Examples
    /// ```
    /// use molecules::molecular_formula::MolecularFormula;
    /// let water = "H2O".parse::<MolecularFormula>().unwrap();
    /// let methane = "CH4".parse::<MolecularFormula>().unwrap();
    /// assert_eq!(water.monoisotopic_mass(), 18.01056468403);
    /// assert_eq!(methane.monoisotopic_mass(), 16.03130012892);
    ///
    /// ```

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Ok(MolecularFormula::new());
        }
        let mut number_buffer = String::new();
        let mut element_buffer = String::new();
        let mut elements: IntMap<u8, usize> = IntMap::default();
        s.chars().for_each(|c| {
            if c.is_alphabetic() {
                if !element_buffer.is_empty() && c.is_uppercase() {
                    let atomic_number = ATOMIC_NUMBERS[&element_buffer];
                    let count = if !number_buffer.is_empty() {
                        number_buffer.parse::<usize>().unwrap_or(1)
                    } else {
                        1
                    };
                    *elements.entry(atomic_number).or_insert(0) += count;
                    number_buffer.clear();
                    element_buffer.clear();
                }
                element_buffer.push(c);
            } else if c.is_numeric() {
                number_buffer.push(c);
            }
        });
        if !element_buffer.is_empty() {
            let atomic_number = ATOMIC_NUMBERS[&element_buffer];
            let count = if !number_buffer.is_empty() {
                number_buffer.parse::<usize>().unwrap_or(1)
            } else {
                1
            };
            *elements.entry(atomic_number).or_insert(0) += count;
        }
        Ok(MolecularFormula { elements })
    }
}

#[derive(Debug, Clone, Default)]
pub struct IsotopicPattern {
    pub buckets: IntMap<usize, (f64, f64)>,
}

impl IsotopicPattern {
    pub fn new(abundances: Vec<f64>, masses: Vec<f64>) -> Self {
        let mut buckets = IntMap::default();
        for (&intensity, &mass) in abundances.iter().zip(masses.iter()) {
            let mass_index = (PRECISION * mass).round() as usize;
            match buckets.entry(mass_index) {
                std::collections::hash_map::Entry::Occupied(mut entry) => {
                    let bucket: &mut (f64, f64) = entry.get_mut();
                    let mass_sum = mass * intensity + bucket.0 * bucket.1;
                    let intensity_sum = intensity + bucket.1;
                    bucket.0 = mass_sum / intensity_sum;
                    bucket.1 += intensity;
                }
                std::collections::hash_map::Entry::Vacant(entry) => {
                    entry.insert((mass, intensity));
                }
            }
        }
        IsotopicPattern { buckets }
    }

    pub fn to_vec(&self) -> Vec<(f64, f64)> {
        let mut vec = self.buckets.values().cloned().collect::<Vec<_>>();
        vec.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        vec
    }

    /// Merge two isotopic patterns into one.
    ///
    /// # Examples
    /// ```
    /// use molecules::molecular_formula::IsotopicPattern;
    /// let mut pattern1 = IsotopicPattern::new(vec![0.5, 0.5], vec![1.0, 2.0]);
    /// let pattern2 = IsotopicPattern::new(vec![0.5, 0.5], vec![3.0, 4.0]);
    /// let merged = pattern1.add(&pattern2);
    /// ```
    pub fn add(&mut self, other: &IsotopicPattern) {
        for (mass_index, &(mass, abundance)) in other.buckets.iter() {
            let mut bucket = *self.buckets.entry(*mass_index).or_insert((mass, 0.0));
            bucket.1 += abundance;
            bucket.0 = (bucket.0 * bucket.1 + mass * abundance) / bucket.1;
        }
    }

    /// Add a single entry to the isotopic pattern, this will merge two buckets if they have the
    /// same mass
    ///
    /// # Examples
    /// ```
    /// use molecules::molecular_formula::IsotopicPattern; let mut pattern = IsotopicPattern::new(vec![0.5, 0.5], vec![1.0, 2.0]);
    /// pattern.add_entry(1.0, 0.5);
    /// assert_eq!(pattern.buckets.len(), 2);
    /// ```
    ///
    ///
    pub fn add_entry(&mut self, mass: f64, abundance: f64) {
        let mass_index = (mass * PRECISION).round() as usize;
        let bucket = self.buckets.entry(mass_index).or_insert((mass, 0.0));
        let new_abundance = bucket.1 + abundance;
        bucket.0 = (bucket.0 * bucket.1 + mass * abundance) / new_abundance;
        bucket.1 = new_abundance;
    }

    pub fn normalize(&mut self) {
        let total_abundance: f64 = self.buckets.values().map(|&(_, a)| a).sum();
        for bucket in self.buckets.values_mut() {
            bucket.1 /= total_abundance;
        }
    }

    pub fn normalize_to_highest(&mut self) {
        let max_abundance: f64 = self
            .buckets
            .values()
            .map(|&(_, a)| a)
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        for bucket in self.buckets.values_mut() {
            bucket.1 /= max_abundance;
        }
    }

    pub fn len(&self) -> usize {
        self.buckets.len()
    }

    pub fn is_empty(&self) -> bool {
        self.buckets.is_empty()
    }

    fn combine(&self, other: &Self, n_peaks: usize) -> Self {
        let mut combined = IsotopicPattern::default();
        let min_mass = (self.buckets.keys().min().unwrap() + other.buckets.keys().min().unwrap())
            as f64
            / PRECISION;
        let max_mass = min_mass + n_peaks as f64;

        for (&_, &(mass_self, abundance_self)) in &self.buckets {
            for (&_, &(mass_other, abundance_other)) in &other.buckets {
                let combined_mass = mass_self + mass_other;
                if combined_mass <= max_mass {
                    let combined_prob = abundance_self * abundance_other;
                    combined.add_entry(combined_mass, combined_prob);
                }
            }
        }

        combined
    }
    pub fn shift(&mut self, shift: f64) {
        let mut new_buckets = IntMap::default();
        for bucket in self.buckets.values_mut() {
            let mass = bucket.0 + shift;
            let mass_index = (PRECISION * mass).round() as usize;
            new_buckets.entry(mass_index).or_insert((mass, 0.0)).1 += bucket.1;
        }
        self.buckets = new_buckets;
    }
}

#[derive(Debug, Clone)]
pub struct Isotope {
    pub mass: f64,
    pub abundance: f64,
}

impl Isotope {
    pub const fn new(mass: f64, abundance: f64) -> Isotope {
        Isotope { mass, abundance }
    }
}

impl From<[f64; 2]> for Isotope {
    fn from(isotope: [f64; 2]) -> Self {
        Isotope {
            mass: isotope[0],
            abundance: isotope[1],
        }
    }
}
