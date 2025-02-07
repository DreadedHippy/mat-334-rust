use crate::{determinant::Determinant, matrix_mult, scalar_mult, Matrix2X2, Matrix3X3, Matrix4X4};

pub struct MatrixInversion {
	a: Vec<Vec<f64>>,
	b: Vec<f64>,
	len: usize
}

impl Determinant for MatrixInversion {}

impl MatrixInversion {
	pub fn new(a: Vec<Vec<f64>>, b: Vec<f64>) -> Self {
		if a.len() <= 1 || a.len() > 4 || a.len() != b.len(){
			panic!("Matrix inversion can only take an n by n matrix; n ∈ [2, 4] ")
		}

		let r = a.len();

		for i in 0..r {
				if a[i].len() != r {
						panic!("Matrix inversion can only take an n by n matrix; n ∈ [2, 4]")
				}
		}

		let len = a.len();

		Self{a, b, len}
	}

	pub fn solve(&mut self) {
		let d = self.determinant();

		self.print_a();
		println!("|A| = {}\n", d);

		println!("Deriving A = minor(A)...\n");
		self.minor();
		self.print_a();

		println!("Deriving A = cofactor(A)...\n");
		self.cofactor();
		self.print_a();

		println!("Deriving A = Adj(A) = [cofactor(A)]ᵀ = ...\n");
		self.transpose();
		self.print_a();

		// converting b into a row matrix;
		let b = self.b.iter().map(|&x| vec![x]).collect::<Vec<_>>();

		println!("X = A⁻¹ • B");
		println!("= (1/|A|) • Adj(A) • B = ...\n");
		
		
		let res = match matrix_mult(&self.a, b) {
			Ok(v) => v,
			Err(e) => {
				eprintln!("{}", e);
				panic!()
			}
		};


		println!("Adj(A) • B = {}", Self::rfp(&res[0], 5));
		
		for i in 1..self.len {
			println!("             {}", Self::rfp(&res[i], 5));
		}



		let answer = scalar_mult(&res, 1.0/d);

		println!("");
		println!("(1/|A|) • Adj(A) • B");
		println!("= (1/{}) • Adj(A) • B = {}", d, Self::rfp(&answer[0], 5));

		for i in 1..self.len {
			println!("                         {}", Self::rfp(&answer[i], 5));
		}

		println!("");
		println!("Therefore...");

		for i in 0..self.len {
			println!("X{} = {:.2}", i + 1, answer[i][0])
		}

		// println!("Answer: ");
		// println!("{:?}", answer)

	}

	pub fn print_a(&self) {
		println!("A = {}",Self::rfp(&self.a[0], 5));

		for i in 1..self.len {
			println!("    {}", Self::rfp(&self.a[i], 5));
		}

		println!("")
	}

	pub fn determinant(&self) -> f64 {
		match self.len {
			2 => {Self::determinant_2x2(&Matrix2X2::from_vec(self.a.clone()).0)},
			3 => {Self::determinant_3x3(&Matrix3X3::from_vec(self.a.clone()).0)},
			4 => {Self::determinant_4x4(&Matrix4X4::from_vec(self.a.clone()).0)},
			_ => {-1.0}
		}
	}

	pub fn minor(&mut self) {
		match self.len {
			2 => {self.minor_2x2()},
			3 => {self.minor_3x3()},
			4 => {self.minor_4x4()},
			_ => {}
		}

	}

	pub fn minor_2x2(&mut self) {
		let [a, b] = [self.a[0][0], self.a[0][1]];
		let [c, d] = [self.a[1][0], self.a[1][1]];

		self.a = vec![
			vec![d, c],
			vec![b, a]
		]
	}

	pub fn minor_3x3(&mut self) {
		let [a, b, c] = [self.a[0][0], self.a[0][1], self.a[0][2]];
		let [d, e, f] = [self.a[1][0], self.a[1][1], self.a[1][2]];
		let [g, h, i] = [self.a[2][0], self.a[2][1], self.a[2][2]];

		let na = Self::determinant_2x2(&[[e, f], [h, i]]);
		let nb = Self::determinant_2x2(&[[d, f], [g, i]]);
		let nc = Self::determinant_2x2(&[[d, e], [g, h]]);
		let nd = Self::determinant_2x2(&[[b, c], [h, i]]);
		let ne = Self::determinant_2x2(&[[a, c], [g, i]]);
		let nf = Self::determinant_2x2(&[[a, b], [g, h]]);
		let ng = Self::determinant_2x2(&[[b, c], [e, f]]);
		let nh = Self::determinant_2x2(&[[a, c], [d, f]]);
		let ni = Self::determinant_2x2(&[[a, b], [d, e]]);

		self.a = vec![
			vec![na, nb, nc],
			vec![nd, ne, nf],
			vec![ng, nh, ni],
		]
	}

	pub fn minor_4x4(&mut self) {
		let [a, b, c, d] = [self.a[0][0], self.a[0][1], self.a[0][2], self.a[0][3]];
		let [e, f, g, h] = [self.a[1][0], self.a[1][1], self.a[1][2], self.a[1][3]];
		let [i, j, k, l] = [self.a[2][0], self.a[2][1], self.a[2][2], self.a[2][3]];
		let [m, n, o, p] = [self.a[3][0], self.a[3][1], self.a[3][2], self.a[3][3]];

		let na = Self::determinant_3x3(&[[f, g, h], [j, k, l], [n, o, p]]);
		let nb = Self::determinant_3x3(&[[e, g, h], [i, k, l], [m, o, p]]);
		let nc = Self::determinant_3x3(&[[e, f, h], [i, j, l], [m, n, p]]);
		let nd = Self::determinant_3x3(&[[e, f, g], [i, j, k], [m, n, o]]);
		let ne = Self::determinant_3x3(&[[b, c, d], [j, k, l], [n, o, p]]);
		let nf = Self::determinant_3x3(&[[a, c, d], [i, k, l], [m, o, p]]);
		let ng = Self::determinant_3x3(&[[a, b, d], [i, j, l], [m, n, p]]);
		let nh = Self::determinant_3x3(&[[a, b, c], [i, j, k], [m, n, o]]);
		let ni = Self::determinant_3x3(&[[b, c, d], [f, g, h], [n, o, p]]);
		let nj = Self::determinant_3x3(&[[a, c, d], [e, g, h], [m, o, p]]);
		let nk = Self::determinant_3x3(&[[a, b, d], [e, f, h], [m, n, p]]);
		let nl = Self::determinant_3x3(&[[a, b, c], [e, f, g], [m, n, o]]);
		let nm = Self::determinant_3x3(&[[b, c, d], [f, g, h], [j, k, l]]);
		let nn = Self::determinant_3x3(&[[a, c, d], [e, g, h], [i, k, j]]);
		let no = Self::determinant_3x3(&[[a, b, d], [e, f, h], [i, j, l]]);
		let np = Self::determinant_3x3(&[[a, b, c], [e, f, g], [i, j, k]]);

		self.a = vec![
			vec![na, nb, nc, nd],
			vec![ne, nf, ng, nh],
			vec![ni, nj, nk, nl],
			vec![nm, nn, no, np],
		]

	}

	pub fn cofactor(&mut self) {
		let n = self.len;

		for i in 0..n {
			for j in 0..n {
				if i+j & 1 == 1 {
					self.a[i][j] = -self.a[i][j];
				}
			}
		}
	}

	pub fn transpose(&mut self) {
		let n = self.len;

		for i in 0..n {
			for j in i..n {
				(self.a[i][j], self.a[j][i]) = (self.a[j][i], self.a[i][j])
			}
		}
	}
}