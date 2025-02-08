use crate::determinant::Determinant;

pub struct LUDecomposition {
	a: Vec<Vec<f64>>,
	b: Vec<f64>,
	len: usize
}

impl Determinant for LUDecomposition {}

impl LUDecomposition {
	pub fn new(a: Vec<Vec<f64>>, b: Vec<f64>) -> Self {
		if a.len() <= 1 || a.len() > 4 || a.len() != b.len(){
			panic!("L-U Decomposition can only take an n by n matrix; n ∈ [2, 4] ")
		}

		let r = a.len();

		for i in 0..r {
			if a[i].len() != r {
				panic!("L-U Decomposition can only take an n by n matrix; n ∈ [2, 4]")
			}
		}

		// Check it is L-U factorizable
		if !Self::is_factorizable(&a) {
			eprintln!("The provided matrix is not L-U Factorizable");
			panic!()
		}

		let len = a.len();
		Self{a, b, len}
	}

	pub fn is_factorizable(a: &Vec<Vec<f64>>) -> bool{
		let n = a.len();

		for i in 0..n {
			match i+1 {
				1 => {
					if a[i][i] == 0.0 { return false}
				},
				2 => {
					let v = [
						[a[0][0], a[0][1]],
						[a[1][0], a[1][1]],
					];

					let d = Self::determinant_2x2(&v);

					if d == 0.0 { return false}
					
				},
				3 => {
					let v = [
						[a[0][0], a[0][1], a[0][2]],
						[a[1][0], a[1][1], a[1][2]],
						[a[2][0], a[2][1], a[2][2]],
					];

					let d = Self::determinant_3x3(&v);

					if d == 0.0 { return false}

				},
				4 => {
					let v = [
						[a[0][0], a[0][1], a[0][2], a[0][3]],
						[a[1][0], a[1][1], a[1][2], a[1][3]],
						[a[2][0], a[2][1], a[2][2], a[2][3]],
						[a[3][0], a[3][1], a[2][2], a[3][3]],
					];

					let d = Self::determinant_4x4(&v);

					if d == 0.0 { return false}
				},
				_ => { return false }
			}
		}

		return true
	}

	pub fn solve(&mut self) {
		// Check it is L-U factorizable
		if !Self::is_factorizable(&self.a) {
			eprintln!("The provided matrix is not L-U Factorizable");
			return
		}

		println!("The provided matrix is L-U factorizable, proceeding to solve...\n");



		match self.len {
			3 => {self.solve_3x3()},
			k => {eprintln!("Cannot solve a {} by {} matrix", k, k)}
		}
		
	}

	pub fn solve_3x3(&mut self) {
		let padding = 7;
		let mut l = vec![vec![0.0; self.len]; self.len];
		let mut u = vec![vec![0.0; self.len]; self.len];

		println!("The given matrix A can be expressed as LU where...\n");
		
		println!("A = |A11 A12 A12 |, L = | L22 0   0   |,  U = | 1   U12 U13 |");
		println!("A = |A21 A22 A23 |  L = | L21 L22 0   |   U = | 0   1   U23 |");
		println!("A = |A31 A32 A33 |  L = | L31 L32 L33 |   U = | 0   0   1   |");

		println!("");
		println!("Multiplying L by U and changing the subject of formulae...\n");

		l[0][0] = self.a[0][0];
		println!("L11 = A11 = {}", l[0][0]);
		l[1][0] = self.a[1][0];
		println!("L21 = A21 = {}", l[1][0]);
		l[2][0] = self.a[2][0];
		println!("L31 = A31 = {}", l[2][0]);
		u[0][1] = self.a[0][1]/l[0][0];
		println!("U12 = A12/L11 = {}", u[0][1]);
		u[0][2] = self.a[0][2]/l[0][0];
		println!("U13 = A13/L11 = {}", u[0][2]);
		l[1][1] = self.a[1][1] - (l[1][0] * u[0][1]);
		println!("L22 = A22 - (L21 • U12) = {}", l[1][1]);
		u[1][2] = (self.a[1][2] - (l[1][0] * u[0][2]))/l[1][1];
		println!("U23 = (A23 - (L21 • U13))/L22 = {}", u[1][2]);
		l[2][1] = self.a[2][1] - (l[2][0] * u[0][1]);
		println!("L32 = A32 - (L31 • U12) = {}", l[2][1]);
		l[2][2] = self.a[2][2] - ((l[2][0] * u[0][2]) + (l[2][1] * u[1][2]));
		println!("L33 = A33 - ((L31 • U13) + (L32 • U23) = {}", l[2][2]);

		println!("");

		for i in 0..3 {
			u[i][i] = 1.0;
		}

		println!("Now, L = ...");

		for row in &l {
			println!("{}", Self::rfp(row, padding));
		}

		println!("");
		println!("And U = ...");
		for row in &u {
			println!("{}", Self::rfp(row, padding));
		}

		// Solving for P
		println!("");
		println!("AX = B");
		println!("∵ A = LU");
		println!("LUX = B");
		println!("Let UX = P");
		println!("Then LP = B");
		println!("");
		println!("Solving for P...");
		println!("");

		println!("Let P = |p|");
		println!("        |q|");
		println!("        |r|");

		println!("");

		println!("Recall B = |{:>5.2}|", self.b[0]);
		println!("           |{:>5.2}|", self.b[1]);
		println!("           |{:>5.2}|", self.b[2]);

		println!("");

		println!("Evaluating from top to bottom...");
		let p = self.b[0]/l[0][0];
		println!("p = B11/L11 = {:>5.2}", p);
		let q = (self.b[1] - (l[1][0] * p))/l[1][1];
		println!("q = (B21 - (L21 • p))/L22 = {:>5.2}", q);
		let r = (self.b[2] - ((l[2][0]* p) + (l[2][1] * q)))/l[2][2];
		println!("r = (B31 - ((L31 • p) + (L32 • q)))/L33 = {:>5.2}", r);

		println!("");

		println!("P = |{:>5.2}|", p);
		println!("    |{:>5.2}|", q);
		println!("    |{:>5.2}|", r);

		println!("");
		println!("Recall UX = P; Solving for X...");
		println!("");

		let p = [p, q, r];

		println!("Let X = |x|");
		println!("        |y|");
		println!("        |z|");

		println!("");

		println!("Evaluating from bottom to top...");
		let z = p[2]/u[2][2];
		println!("z = P31/U33 = {:>5.2}", z);
		let y = (p[1] - (u[1][2] * z))/u[1][1];
		println!("y = (P21 - (U23 • z))/U22 = {:>5.2}", y);
		let x = (p[0] - ((u[0][1] * y) + (u[0][2]* z)))/u[0][0];
		println!("x = (P11 - ((U12 • y) + (U13 • z)))/U11 = {:>5.2}", x);
		println!("");

		println!("X = |{:>5}|", format!("{:.2}", x));
		println!("    |{:>5}|", format!("{:.2}", y));
		println!("    |{:>5}|", format!("{:.2}", z));

	}
}
