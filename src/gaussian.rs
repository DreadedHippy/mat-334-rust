use crate::{Matrix2X2, Matrix3X3};
#[allow(unused_imports)]
use crate::{format_as_matrix_row, rf, solve_4x4, Matrix4X4};

pub struct Gaussian {
	a: Vec<Vec<f64>>,
	b: Vec<f64>
}

impl Gaussian {
	pub fn new(a: Vec<Vec<f64>>, b: Vec<f64>) -> Self {
        
		if a.len() <= 1 || a.len() > 4 || a.len() != b.len(){
				panic!("Gaussian Elimination can only take an n by n matrix; n ∈ [2, 4] ")
		}

		let r = a.len();

		for i in 0..r {
				if a[i].len() != r {
						panic!("Gaussian Elimination can only take an n by n matrix; n ∈ [2, 4]")
				}
		}

		
		Self{a, b}
	}

	pub fn solve(&self) {
		let n = self.a.len();

		match n {
			4 => Self::solve_4x4(Matrix4X4::from_vec(self.a.clone()).0, self.b.clone().try_into().unwrap()),
			3 => Self::solve_3x3(Matrix3X3::from_vec(self.a.clone()).0, self.b.clone().try_into().unwrap()),
			2 => Self::solve_2x2(Matrix2X2::from_vec(self.a.clone()).0, self.b.clone().try_into().unwrap()),
			_ => {}
		}
	}

	pub fn c(i: f64) -> String {
		if i == 1.0 {String::new()} else if i == -1.0 {"-".to_string()} else {i.to_string()}
	}

	pub fn solve_4x4(mut a: [[f64; 4]; 4], mut b: [f64; 4]) {
		let mut y = 0;
		// step 1
		// let (r2, r3, r4) = (lcm(a[1][y], a[y][y]), lcm(a[2][y], a[y][y]), lcm(a[3][y], a[y][y]));
		let r2c = (a[1][y], a[y][y]);
		let r3c = (a[2][y], a[y][y]);
		let r4c = (a[3][y], a[y][y]);

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		for i in (y+1)..4 {
			println!("{} ¦ {} <-- R{} = {}R{} - {}R{}", rf(&a[i]), rf(&[b[i]]), i+1, Self::c(a[i][y]), 1, Self::c(a[y][y]), i + 1)
		}

		for i in y..4 {
			a[1][i] = (r2c.0 * a[y][i]) - (r2c.1 * a[1][i]);
			a[2][i] = (r3c.0 * a[y][i]) - (r3c.1 * a[2][i]);
			a[3][i] = (r4c.0 * a[y][i]) - (r4c.1 * a[3][i]);
		}
		b[1] = (r2c.0 * b[y]) - (r2c.1 * b[1]);
		b[2] = (r3c.0 * b[y]) - (r3c.1 * b[2]);
		b[3] = (r4c.0 * b[y]) - (r4c.1 * b[3]);
		
		
		// step 2
		y = 1;
		let r3c = (a[2][y], a[y][y]);
		let r4c = (a[3][y], a[y][y]);

		println!("");

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		for i in (y+1)..4 {
			println!("{} ¦ {} <-- R{} = {}R{} - {}R{}", rf(&a[i]), rf(&[b[i]]), i+1, Self::c(a[i][y]), y+1, Self::c(a[y][y]), i + 1);
		}

		for i in y..4 {
			a[2][i] = (r3c.0 * a[y][i]) - (r3c.1 * a[2][i]);
			a[3][i] = (r4c.0 * a[y][i]) - (r4c.1 * a[3][i]);
		}
		b[2] = (r3c.0 * b[y]) - (r3c.1 * b[2]);
		b[3] = (r4c.0 * b[y]) - (r4c.1 * b[3]);

				
		// step 3
		y = 2;
		let r4c = (a[3][y], a[y][y]);

		

		println!("");

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		for i in (y+1)..4 {
			println!("{} ¦ {} <-- R{} = {}R{} - {}R{}", rf(&a[i]), rf(&[b[i]]), i+1, Self::c(a[i][y]), y+1, Self::c(a[y][y]), i + 1)
		}

		for i in y..4 {
			a[3][i] = (r4c.0 * a[y][i]) - (r4c.1 * a[3][i]);
		}
		b[3] = (r4c.0 * b[y]) - (r4c.1 * b[3]);

		println!("");

		y += 1;

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		println!("");

		
		
		let x4 = b[3]/a[3][3];
		let x3 = (b[2] - (a[2][3] * x4))/a[2][2];
		let x2 = (b[1] - ((a[1][3] * x4) + (a[1][2] * x3)))/a[1][1];
		let x1 = (b[0] - ((a[0][3] * x4) + (a[0][2] * x3) + (a[0][1] * x2)))/a[0][0];

		println!("Evaluating from bottom to top...");

		let solves = [
			format!("X1 = ({} - (({})X2 + ({})X3 + ({})X4))/{} = {}", b[0], a[0][3], a[0][2], a[0][1], a[0][0], x1),
			format!("X2 = ({} - (({})X3 + ({})X4))/{} = {}", b[1], a[1][2], a[1][2], a[1][1], x2),
			format!("X3 = ({} - ({})X4)/{} = {}", b[2], a[2][3], a[2][2], x3),
			format!("X4 = {}/{} = {}", b[3], a[3][3], x4),
		];

		for i in (0..=y).rev() {
			println!("{} | {} ===> {}", rf(&a[i]), rf(&[b[i]]), solves[i]);
		}

		println!("");

		println!("X1 = {}", x1);
		println!("X2 = {}", x2);
		println!("X3 = {}", x3);
		println!("X4 = {}", x4);
	}

	

	pub fn solve_3x3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) {
		let n = 3;
		let mut y = 0;
		// step 1
		let r2c = (a[1][y], a[y][y]);
		let r3c = (a[2][y], a[y][y]);

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		for i in (y+1)..n {
			println!("{} ¦ {} <-- R{} = {}R{} - {}R{}", rf(&a[i]), rf(&[b[i]]), i+1, Self::c(a[i][y]), 1, Self::c(a[y][y]), i + 1)
		}

		for i in y..n {
			a[1][i] = (r2c.0 * a[y][i]) - (r2c.1 * a[1][i]);
			a[2][i] = (r3c.0 * a[y][i]) - (r3c.1 * a[2][i]);
		}
		b[1] = (r2c.0 * b[y]) - (r2c.1 * b[1]);
		b[2] = (r3c.0 * b[y]) - (r3c.1 * b[2]);
		
		
		// step 2
		y = 1;
		let r3c = (a[2][y], a[y][y]);

		println!("");

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		for i in (y+1)..n {
			println!("{} ¦ {} <-- R{} = {}R{} - {}R{}", rf(&a[i]), rf(&[b[i]]), i+1, Self::c(a[i][y]), y+1, Self::c(a[y][y]), i + 1);
		}

		for i in y..n {
			a[2][i] = (r3c.0 * a[y][i]) - (r3c.1 * a[2][i]);
		}
		b[2] = (r3c.0 * b[y]) - (r3c.1 * b[2]);
		println!("");

		y += 1;

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		println!("");


		
		let x3 = b[2]/a[2][2];
		let x2 = (b[1] - (a[1][2] * x3))/a[1][1];
		let x1 = (b[0] - ((a[0][2] * x3) + (a[0][1] * x2)))/a[0][0];

		println!("Evaluating from bottom to top...");

		let solves = [
			format!("X1 = ({} - (({})X2 + ({})X3)/{} = {}", b[0], a[0][1], a[0][2], a[0][0], x1),
			format!("X2 = ({} - ({})X3)/{} = {}", b[1], a[1][2], a[1][1], x2),
			format!("X3 = {}/{} = {}", b[2], a[2][2], x3),
		];

		for i in (0..=y).rev() {
			println!("{} ¦ {} ===> {}", rf(&a[i]), rf(&[b[i]]), solves[i]);
		}

		println!("");

		println!("X1 = {}", x1);
		println!("X2 = {}", x2);
		println!("X3 = {}", x3);
		// println!("X4 = {}", x4);
	}

	pub fn solve_2x2(mut a: [[f64; 2]; 2], mut b: [f64; 2]) {
		let n = 2;
		let mut y = 0;
		// step 1
		let r2c = (a[1][y], a[y][y]);

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		for i in (y+1)..n {
			println!("{} ¦ {} <-- R{} = {}R{} - {}R{}", rf(&a[i]), rf(&[b[i]]), i+1, Self::c(a[i][y]), 1, Self::c(a[y][y]), i + 1)
		}

		for i in y..n {
			a[1][i] = (r2c.0 * a[y][i]) - (r2c.1 * a[1][i]);
		}
		b[1] = (r2c.0 * b[y]) - (r2c.1 * b[1]);
		
		
		// step 2
		y = 1;

		println!("");

		for i in 0..=y {
			println!("{} ¦ {}", rf(&a[i]), rf(&[b[i]]));
		}

		println!("");


		
		let x2 = b[1]/a[1][1];
		let x1 = (b[0] - (a[0][1] * x2))/a[0][0];
		// let x1 = (b[0] - ((a[0][2] * x3) + (a[0][1] * x2)))/a[0][0];

		println!("Evaluating from bottom to top...");

		let solves = [
			format!("X1 = ({} - ({})X3)/{} = {}", b[0], a[0][1], a[0][0], x1),
			format!("X2 = {}/{} = {}", b[1], a[1][1], x2),
		];

		for i in (0..=y).rev() {
			println!("{} ¦ {} ===> {}", rf(&a[i]), rf(&[b[i]]), solves[i]);
		}

		println!("");

		println!("X1 = {}", x1);
		println!("X2 = {}", x2);
		// println!("X4 = {}", x4);
	}
}