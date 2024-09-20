use rand::Rng;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Sub, Neg};

/// A 3D vector represented by its x, y, and z components.
///
/// # Example
///
/// ```
/// use molecules::vector::Vector;
///
/// let v1 = Vector::new(1.0, 2.0, 3.0);
/// let v2 = Vector::new(4.0, 5.0, 6.0);
/// let v3 = v1 + v2;
/// let v4 = v1 * 2.0;
/// let v5 = v1 / 2.0;
/// println!("Vector sum: {:?}", v3);
/// println!("Vector scalar multiplication: {:?}", v3);
/// println!("Vector scalar division: {:?}", v3);
/// ```
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
// Implement the Add trait for Vector
impl Add<Vector> for Vector {
    type Output = Vector;

    fn add(self, other: Vector) -> Vector {
        Vector {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

// Implement the Mul trait for Vector
impl Mul<f64> for Vector {
    type Output = Vector;

    fn mul(self, scalar: f64) -> Vector {
        Vector {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }
}

// Implement the Div trait for Vector
impl Div<f64> for Vector {
    type Output = Vector;

    fn div(self, scalar: f64) -> Self::Output {
        Vector {
            x: self.x / scalar,
            y: self.y / scalar,
            z: self.z / scalar,
        }
    }
}
// Implement the Sub trait for Vector
impl Sub<Vector> for Vector {
    type Output = Vector;
    fn sub(self, rhs: Vector) -> Self::Output {
        Vector::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

// Implement the AddAssign trait for Vector
impl AddAssign for Vector {
    fn add_assign(&mut self, other: Self) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

// Implement the AddAssign trait for Vector
impl MulAssign<f64> for Vector {
    fn mul_assign(&mut self, scalar: f64) {
        self.x = self.x * scalar;
        self.y = self.y * scalar;
        self.z = self.z * scalar;
    }
}

impl Default for Vector {
    /// Creates a default `Vector` being the zero vector
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v = Vector::default();
    /// println!("Vector: {:?}", v);
    /// assert_eq!(v.length(),0.0)
    /// ```
    fn default() -> Vector {
        Vector::new(0.0, 0.0, 0.0)
    }
}

impl Vector {
    /// Creates a `Vector` from an array of 3 elements, representing the x, y, and z components.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let array = [1.0, 2.0, 3.0];
    /// let v = Vector::from_vec(&array);
    /// println!("Vector: {:?}", v);
    /// ```
    pub fn from_vec(array: &[f64; 3]) -> Vector {
        Vector {
            x: array[0],
            y: array[1],
            z: array[2],
        }
    }
    /// Creates a `Vector` from x, y, and z components.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v = Vector::new(1.0, 2.0, 3.0);
    /// println!("Vector: {:?}", v);
    /// ```
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn x() -> Self {
        Self::new(1.0, 0.0, 0.0)
    }
    pub fn y() -> Self {
        Self::new(0.0, 1.0, 0.0)
    }
    pub fn z() -> Self {
        Self::new(0.0, 0.0, 1.0)
    }
    /// Calculates the difference between two `Vector` objects, this is equivalent to the vector
    /// pointing from `other` to `self`.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v1 = Vector::new(1.0, 1.0, 1.0);
    /// let v2 = Vector::new(1.0, 1.0, -1.0);
    ///
    /// let v_diff = v1.difference(&v2);
    /// assert_eq!(v_diff, Vector::new(0.0, 0.0, 2.0));
    /// assert_eq!(v_diff.length(), 2.0);
    /// ```
    pub fn difference(&self, other: &Self) -> Self {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        Self::new(dx, dy, dz)
    }
    /// Calculates the length of a vector.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    ///     let v = Vector::new(1.0, 1.0, 1.0);
    ///     println!("Vector: {:?}", v.length());
    ///     assert_eq!(v.length(), 3.0_f64.sqrt())
    ///
    /// ```
    pub fn length(&self) -> f64 {
        (self.squared_length()).sqrt()
    }

    /// Alias for length. Returns the magnitude of the vector.
    ///
    /// # Returns
    ///
    /// * `f64` - The magnitude (length) of the vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v = Vector::new(3.0, 4.0, 0.0);
    /// assert_eq!(v.magnitude(), 5.0);
    /// ```
    pub fn magnitude(&self) -> f64 {
        self.length()
    }
    
    /// Calculates the squared length of the vector.
    ///
    /// This is more efficient than `length()` when you only need to compare lengths,
    /// as it avoids the square root calculation.
    ///
    /// # Returns
    ///
    /// * `f64` - The squared length of the vector.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v = Vector::new(1.0, 2.0, 2.0);
    /// assert_eq!(v.squared_length(), 9.0);
    /// ```
    pub fn squared_length(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    /// Determines if two vectors are parallel within a given tolerance.
    ///
    /// # Arguments
    ///
    /// * `other` - The other vector to compare with.
    /// * `tolerance` - The tolerance for parallelism. Smaller values require more precise parallelism.
    ///
    /// # Returns
    ///
    /// * `bool` - True if the vectors are parallel within the given tolerance, false otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v1 = Vector::new(1.0, 0.0, 0.0);
    /// let v2 = Vector::new(2.0, 0.0, 0.0);
    /// assert!(v1.are_vectors_parallel(&v2, 1e-10));
    ///
    /// let v3 = Vector::new(1.0, 0.1, 0.0);
    /// assert!(!v1.are_vectors_parallel(&v3, 1e-10));
    /// ```
    pub fn are_vectors_parallel(&self, other: &Vector, tolerance: f64) -> bool {
        let cross_product = self.cross(other);
        cross_product.length() < tolerance * self.length() * other.length()
    }


    /// Calculates the dot product of two vectors.
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v1 = Vector::new(1.0, 2.0, 3.0);
    /// let v2 = Vector::new(4.0, 5.0, 6.0);
    /// println!("Vector: {:?}", v1.dot(&v2));
    ///
    /// ```
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    /// Calculates the distance between to vectors
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v1 = Vector::new(1.0, 2.0, 3.0);
    /// let v2 = Vector::new(4.0, 5.0, 6.0);
    /// println!("Vector: {:?}", v1.distance(&v2));
    ///
    /// ```
    pub fn distance(&self, other: &Self) -> f64 {
        self.difference(other).length()
    }
    /// Returns the squared distance between two vectors
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::vector::Vector;
    /// let v1 = Vector::new(1.0, 2.0, 3.0);
    /// let v2 = Vector::new(4.0, 5.0, 6.0);
    /// assert_eq!(v1.distance_squared(&v2),27.0);
    /// ```
    pub fn distance_squared(&self, other: &Self) -> f64 {
        let difference = *other - *self;
        difference.x.powi(2) + difference.y.powi(2) + difference.z.powi(2)
    }
    /// Generates a random unit vector
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v = Vector::random_unit_vector();
    /// println!("Vector: {:?}", v);
    /// assert_eq!(v.length().round(),1.0)
    /// ```
    ///
    pub fn random_unit_vector() -> Self {
        let mut rng = rand::thread_rng();
        let theta = rng.gen_range(0.0..std::f64::consts::PI);
        let phi = rng.gen_range(0.0..2.0 * std::f64::consts::PI);

        let x = theta.sin() * phi.cos();
        let y = theta.sin() * phi.sin();
        let z = theta.cos();

        Vector::new(x, y, z)
    }
    /// Calculates the cross product of two vectors and returns a new vector
    ///
    /// # Panics
    ///
    /// This function does not panic.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let vec1 = Vector::x();
    /// let vec2 = Vector::y();
    /// assert_eq!(vec1.cross(&vec2),Vector::z())
    /// ```
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
    /// Returns the angle between two vectors in radians.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let vec_a = Vector::new(1.0, 0.0, 0.0);
    /// let vec_b = Vector::new(0.0, 1.0, 0.0);
    ///
    /// let angle = vec_a.angle_between(&vec_b);
    /// assert_eq!(angle, std::f64::consts::FRAC_PI_2);
    /// ```
    pub fn angle_between(&self, other: &Vector) -> f64 {
        let lengths_product = self.length() * other.length();
        if lengths_product == 0.0 {
            0.0
        } else {
            let angle_cosine = (self.dot(other) / (lengths_product)).clamp(-1.0, 1.0);
            angle_cosine.acos()
        }
    }
    /// Calculates a normal vector for a given input vector
    ///
    /// # Panics
    /// This function does not panic as the case of a zero length vector is handeled.
    ///
    /// # Example
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let vec1 = Vector::new(1.0,2.0,3.0);
    /// assert_eq!(vec1.normalize().length().ceil(),1.0)
    /// ```
    pub fn normalize(&self) -> Self {
        if self.length() == 0.0 {
            *self
        } else {
            *self / self.length()
        }
    }

    /// Calculates the angle given three points in the order, start, middle, end
    ///
    /// # Panics
    /// Function does not panic
    ///
    /// # Example
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let start_point = Vector::new(1.0,0.0,0.0);
    /// let middle_point = Vector::new(0.0,0.0,0.0);
    /// let end_point = Vector::new(0.0,1.0,0.0);
    /// assert_eq!(start_point.angle_between_points(&middle_point,&end_point),1.5707963267948966);
    /// ```
    pub fn angle_between_points(&self, middle_point: &Vector, end_point: &Vector) -> f64 {
        let v1 = *middle_point - *self;
        let v2 = *end_point - *middle_point;
        v1.angle_between(&v2)
    }
    /// Calculates the triple product of three vectors
    ///
    /// # Example
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let start_point = Vector::new(1.0,0.0,0.0);
    /// let middle_point = Vector::new(0.0,1.0,0.0);
    /// let end_point = Vector::new(0.0,0.0,1.0);
    /// // Should be ~1 (right handed)
    /// let triple_product = start_point.triple_product(&middle_point, &end_point);
    /// assert!(triple_product - 1.0 < 0.01);
    ///
    /// let start_point = Vector::new(1.0,0.0,0.0);
    /// let middle_point = Vector::new(0.0,0.0,1.0);
    /// let end_point = Vector::new(0.0,1.0,0.0);
    ///
    /// // Should be ~ -1 (right handed)
    /// let triple_product = start_point.triple_product(&middle_point, &end_point);
    /// assert!(triple_product + 1.0 < 0.01);
    /// ```
    pub fn triple_product(&self, other1: &Vector, other2: &Vector) -> f64 {
        self.dot(&other1.cross(other2))
    }
    pub fn as_tuple(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }
    pub fn as_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
    pub fn as_vec(&self) -> Vec<f64> {
        vec![self.x, self.y, self.z]
    }

    /// Determines if two vectors are pointing in opposite directions within a given tolerance.
    ///
    /// # Arguments
    ///
    /// * `other` - The other vector to compare with.
    /// * `tolerance` - The tolerance for the comparison. Smaller values require more precise opposition.
    ///
    /// # Returns
    //
    /// * `bool` - True if the vectors are opposite within the given tolerance, false otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v1 = Vector::new(1.0, 0.0, 0.0);
    /// let v2 = Vector::new(-1.0, 0.0, 0.0);
    /// assert!(v1.is_opposite(&v2, 1e-10));
    ///
    /// let v3 = Vector::new(-0.99, 0.01, 0.0);
    /// assert!(!v1.is_opposite(&v3, 1e-10));
    /// assert!(v1.is_opposite(&v3, 0.02));
    /// ```
    pub fn is_opposite(&self, other: &Vector, tolerance: f64) -> bool {
        let dot_product = self.dot(other);
        let magnitudes_product = self.length() * other.length();
        
        if magnitudes_product == 0.0 {
            false // At least one vector is zero, so they can't be opposite
        } else {
            let cos_angle = dot_product / magnitudes_product;
            (cos_angle + 1.0).abs() < tolerance
        }
    }
}

// Implement the Neg trait for Vector
impl Neg for Vector {
    type Output = Vector;

    /// Negates the vector, reversing its direction.
    ///
    /// # Returns
    ///
    /// * `Vector` - A new Vector with all components negated.
    ///
    /// # Examples
    ///
    /// ```
    /// use molecules::vector::Vector;
    ///
    /// let v = Vector::new(1.0, -2.0, 3.0);
    /// let negated = -v;
    /// assert_eq!(negated, Vector::new(-1.0, 2.0, -3.0));
    /// ```
    fn neg(self) -> Self::Output {
        Vector {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}
