pub trait Determinant {
	
	fn determinant_4x4(matrix: &[[f64; 4]; 4]) -> f64{
		let (a, b, c, d) = (matrix[0][0], matrix[0][1], matrix[0][2], matrix[0][3]);
		let (e, f, g, h) = (matrix[1][0], matrix[1][1], matrix[1][2], matrix[1][3]);
		let (i, j, k, l) = (matrix[2][0], matrix[2][1], matrix[2][2], matrix[2][3]);
		let (m, n, o, p) = (matrix[3][0], matrix[3][1], matrix[3][2], matrix[3][3]);

		let ax = [[f, g, h], [j, k, l], [n, o, p]];
		let bx = [[e, g, h], [i, k, l], [m, o, p]];
		let cx = [[e, f, h], [i, j, l], [m, n, p]];
		let dx = [[e, f, g], [i, j, k], [m, n, o]];		

		let (ar, br, cr, dr) = (Self::determinant_3x3(&ax), Self::determinant_3x3(&bx), Self::determinant_3x3(&cx), Self::determinant_3x3(&dx));

		let res = (a * ar) - (b * br) + ( c * cr) - (d * dr);

		println!("");
		return res;
		
	}

	#[allow(unused)]
	fn rf(row: &[f64]) -> String {
		Self::format_as_matrix_row(row)
	}
	
	#[allow(unused)]
	fn format_as_matrix_row(row: &[f64]) -> String {
		return format!("|{}|", row.iter().map(|x| format!("{x:>3} ")).collect::<String>().trim_end())
	}

	
	fn rfp(row: &[f64], p: usize) -> String {
		Self::format_as_matrix_row_p(row, p)
	}

	fn format_as_matrix_row_p(row: &[f64], p: usize) -> String {
		let format = row
			.iter()
			.map(|x|
				if x.fract() == 0.0 {
					format!("{}", x)
				} else {
					format!("{:.2}", x)
				}
			)
			.map(|x| format!("{x:>p$}", p = p)).collect::<String>().trim_end().to_string();

		return format!("|{}|", format)
	}

	fn determinant_3x3(matrix: &[[f64; 3]; 3]) -> f64{
		let (a, b, c) = (matrix[0][0], matrix[0][1], matrix[0][2]);
		let (d, e, f) = (matrix[1][0], matrix[1][1], matrix[1][2]);
		let (g, h, i) = (matrix[2][0], matrix[2][1], matrix[2][2]);

		let ax = [[e, f], [h, i]];
		let bx = [[d, f], [g, i]];
		let cx = [[d, e], [g, h]];

		let (ar, br, cr) = (Self::determinant_2x2(&ax), Self::determinant_2x2(&bx), Self::determinant_2x2(&cx));

		let res = (a * ar) - (b * br) + ( c * cr);

		return res;
	}

	fn determinant_2x2(matrix: &[[f64; 2]; 2]) -> f64{
		let (a, b) = (matrix[0][0], matrix[0][1]);
		let (c, d) = (matrix[1][0], matrix[1][1]);

		let res = (a * d) - (b * c);

		return res;
	}
}