use std::{io::{Read, Write}, panic};

use gaussian::Gaussian;
use inversion::MatrixInversion;
use lu_decomposition::LUDecomposition;
#[allow(dead_code)]

mod gaussian;
mod inversion;
mod determinant;
mod lu_decomposition;

fn main() {
    // println!("Hello, world!");
    let args = std::env::args().collect::<Vec<_>>();

    if args.len() >= 2 {
        run(args[1].to_owned(), args[2].to_owned());
    } else {
        let mut command = String::new();
        print!("Input a command: ");
        std::io::stdout().flush().expect("Unable to flush stdout");

        std::io::stdin().read_line(&mut command).expect("Unable to read input");

        let mut n = String::new();
        print!("Input the value of n for your n x n matrix: ");
        std::io::stdout().flush().expect("Unable to flush stdout");

        std::io::stdin().read_line(&mut n).expect("Unable to read input");
        run(command.trim().to_string(), n.trim().to_string());
    }
    

        
    println!("----------------------");
    print!("Press any key to exit!");
    std::io::stdout().flush().expect("Unable to flush stdout");
    let _ = std::io::stdin().read(&mut []).expect("Unable to read input");
}

// fn run (args[1])
pub fn check_invalid_matrix_size(size_arg: String) -> Result<usize, ()> {
    
    let size = size_arg.parse::<usize>();
    
    if size.is_err(){
        eprintln!("Pass in a valid matrix size [1, 4]");
        return Err(())
    }

    let n = size.unwrap();

    if n < 2 || n > 4 {
        eprintln!("Pass in a valid matrix size [2, 4]");
        return Err(())
    }

    Ok(n)
}

pub fn get_b_input(size_arg: String) -> Vec<f64> {
    let n = if let Ok(size) = check_invalid_matrix_size(size_arg) {
        size
    } else { panic!() };

    print!("Input 'b' as a row of numbers: ");
    let _ = std::io::stdout().flush();

    

    let mut input = String::new();
    let _ = std::io::stdin().read_line(&mut input);

    let v = input.split_whitespace().map(|x| x.parse::<f64>().expect("Invalid integer passed")).collect::<Vec<f64>>();
    if v.len() != n {
        eprintln!("Your input must have a length of {}", n);
        panic!()
    }

    return v


}

pub fn get_input(size_arg: String) -> Vec<Vec<f64>>{
    let n = if let Ok(size) = check_invalid_matrix_size(size_arg) {
        size
    } else { panic!() };

    let mut matrix = Vec::new();
    for i in 0..n {
        print!("Input row {}: ", i + 1);
        let _ = std::io::stdout().flush();

        let mut input = String::new();
        let _ = std::io::stdin().read_line(&mut input);

        let v = input.split_whitespace().map(|x| x.parse::<f64>().expect("Invalid integer passed")).collect::<Vec<f64>>();
        if v.len() != n {
            eprintln!("Your input must have a length of {}", n);
            panic!()
        }

        matrix.push(v);
    }

    return matrix;
}
pub struct Cramer {
    a: Vec<Vec<f64>>,
    b: Vec<f64>
}

impl Cramer {
    pub fn new(a: Vec<Vec<f64>>, b: Vec<f64>) -> Self {
        
        if a.len() <= 1 || a.len() > 4 || a.len() != b.len(){
            panic!("Cramer can only take an n by n matrix; n ∈ [2, 4] ")
        }

        let r = a.len();

        for i in 0..r {
            if a[i].len() != r {
                panic!("Cramer can only take an n by n matrix; n ∈ [2, 4]")
            }
        }

        
        Self{a, b}
    }

    pub fn solve(&self) {
        let n = self.a.len();
        let det = Self::get_determinant(self.a.clone());

        let mut x = vec![-1.0; n];

        for i in 0..n {
            let current = self.get_a_b_substitution(i);
            
            let detx = Self::get_determinant(current);
            
            let res = detx/det;
            println!("X{} = det(A{})/det(A) = {}/{} = {}", i + 1, i +1, detx, det, res);
            x[i] = res;

        }

        println!("Result:");
        println!("{:#?}", x);
    }

    fn get_determinant(matrix: Vec<Vec<f64>>) -> f64{
        let n = matrix.len();

        let e = match n {
            2 => {solve_2x2(&Matrix2X2::from_vec(matrix).0)},
            3 => {solve_3x3(&Matrix3X3::from_vec(matrix).0)},
            4 => {solve_4x4(&Matrix4X4::from_vec(matrix).0)},
            _ => { -1.0}
        };

        return e
    }

    fn get_a_b_substitution(&self, index: usize) -> Vec<Vec<f64>> {
        let mut new = self.a.clone();
        let n = new.len();

        for i in 0..n {
            new[i][index] = self.b[i]
        }

        return new
    }
}


pub fn run(command: String, n: String) {
    
        
    match command.as_str() {
        "c" | "cramer"  => {
            let a = get_input(n.to_owned());
            let b = get_b_input(n.to_owned());
            
            println!("");
            
            // println!("{:?} {:?}", a, b);
            let c = Cramer::new(a, b);
            
            c.solve();
        },
        "g" | "gauss" | "gaussian" => {
            let a = get_input(n.to_owned());
            let b = get_b_input(n.to_owned());
            
            println!("");
            
            // println!("{:?} {:?}", a, b);
            let r = Gaussian::new(a, b);
            
            r.solve();
            
        },
        "i" | "inversion" | "inverse" => {
            let a = get_input(n.to_owned());
            let b = get_b_input(n.to_owned());
            
            println!("");
            
            let mut r = MatrixInversion::new(a, b);
            
            r.solve();
            
        },
        "f" | "lu" | "factorization" | "lud" => {
            let a = get_input(n.to_owned());
            let b = get_b_input(n.to_owned());
            
            println!("");
            
            let mut r = LUDecomposition::new(a, b);
            
            r.solve();
            
        },
        k => {
            println!("{k:?}");
            eprintln!("Invalid command")
        }
    }


}

pub struct Matrix4X4(pub [[f64; 4]; 4]);
pub struct Matrix3X3(pub [[f64; 3]; 3]);
pub struct Matrix2X2(pub [[f64; 2]; 2]);

impl Matrix2X2 {
    pub fn from_vec(v: Vec<Vec<f64>>) -> Self {
        Self([
            [v[0][0], v[0][1]],
            [v[1][0], v[1][1]]
        ])
    }
}


impl Matrix3X3 {
    pub fn from_vec(v: Vec<Vec<f64>>) -> Self {
        Self([
            [v[0][0], v[0][1], v[0][2]],
            [v[1][0], v[1][1], v[1][2]],
            [v[2][0], v[2][1], v[2][2]],
        ])
    }
}

impl Matrix4X4 {
    pub fn from_vec(v: Vec<Vec<f64>>) -> Self {
        Self([
            [v[0][0], v[0][1], v[0][2], v[0][3]],
            [v[1][0], v[1][1], v[1][2], v[1][3]],
            [v[2][0], v[2][1], v[2][2], v[2][3]],
            [v[3][0], v[3][1], v[3][2], v[3][3]],
        ])
    }
}

pub fn solve_4x4(matrix: &[[f64; 4]; 4]) -> f64{
    let (a, b, c, d) = (matrix[0][0], matrix[0][1], matrix[0][2], matrix[0][3]);
    let (e, f, g, h) = (matrix[1][0], matrix[1][1], matrix[1][2], matrix[1][3]);
    let (i, j, k, l) = (matrix[2][0], matrix[2][1], matrix[2][2], matrix[2][3]);
    let (m, n, o, p) = (matrix[3][0], matrix[3][1], matrix[3][2], matrix[3][3]);

    let ax = [[f, g, h], [j, k, l], [n, o, p]];
    let bx = [[e, g, h], [i, k, l], [m, o, p]];
    let cx = [[e, f, h], [i, j, l], [m, n, p]];
    let dx = [[e, f, g], [i, j, k], [m, n, o]];

    println!("{} = {:>3} × {} - {:>3} × {} + {:>3} × {} - {:>3} × {}", rf(&matrix[0]), a, rf(&ax[0]), b, rf(&bx[0]), c, rf(&cx[0]), d, rf(&dx[0]));
    
    for i in 1..3 {
        println!("{}   {:<3}   {}   {:<3}   {}   {:<3}   {}   {:<3}   {}", rf(&matrix[i]), "", rf(&ax[i]), "", rf(&bx[i]), "", rf(&cx[i]), "", rf(&dx[i]));
    }

    println!("{}", rf(&matrix[3]));

    println!("");

    

    let (ar, br, cr, dr) = (solve_3x3(&ax), solve_3x3(&bx), solve_3x3(&cx), solve_3x3(&dx));

    let res = (a * ar) - (b * br) + ( c * cr) - (d * dr);

    println!("{} = ({} × {}) - ({} × {}) + ({} × {}) - ({} × {}) = {}", rf(&matrix[0]), a, ar, b, br, c, cr, d, dr, res);

    for i in 1..=3 {
        println!("{}", rf(&matrix[i]))
    }

    println!("");
    return res;
    
}

fn rf(row: &[f64]) -> String {
    format_as_matrix_row(row)
}

fn format_as_matrix_row(row: &[f64]) -> String {
    return format!("|{}|", row.iter().map(|x| format!("{x:>3} ")).collect::<String>().trim_end())
}

pub fn solve_3x3(matrix: &[[f64; 3]; 3]) -> f64{
    let (a, b, c) = (matrix[0][0], matrix[0][1], matrix[0][2]);
    let (d, e, f) = (matrix[1][0], matrix[1][1], matrix[1][2]);
    let (g, h, i) = (matrix[2][0], matrix[2][1], matrix[2][2]);

    let ax = [[e, f], [h, i]];
    let bx = [[d, f], [g, i]];
    let cx = [[d, e], [g, h]];

    println!("{} = {:>3} × {} - {:>3} × {} + {:>3} × {}", rf(&matrix[0]), a, rf(&ax[0]), b, rf(&bx[0]), c, rf(&cx[0]));
    
    for i in 1..2 {
        println!("{}   {:>3}   {}   {:>3}   {}   {:>3}   {}", rf(&matrix[i]), "", rf(&ax[i]), "", rf(&bx[i]), "", rf(&cx[i]));
    }

    println!("{}", rf(&matrix[2]));

    println!("");

    let (ar, br, cr) = (solve_2x2(&ax), solve_2x2(&bx), solve_2x2(&cx));

    let res = (a * ar) - (b * br) + ( c * cr);

    println!("{} = ({} × {}) - ({} × {}) + ({} × {}) = {}", rf(&matrix[0]), a, ar, b, br, c, cr, res);

    for i in 1..=2 {
        println!("{}", rf(&matrix[i]))
    }

    println!("");
    return res;
}

pub fn solve_2x2(matrix: &[[f64; 2]; 2]) -> f64{
    let (a, b) = (matrix[0][0], matrix[0][1]);
    let (c, d) = (matrix[1][0], matrix[1][1]);

    let res = (a * d) - (b * c);
    println!("{} = ({} × {}) - ({} × {}) = {}", rf(&matrix[0]), a, d, b, c, res);
    println!("{}", rf(&matrix[1]));

    println!("");

    return res;

}

pub fn matrix_mult(a: &Vec<Vec<f64>>, b: Vec<Vec<f64>>) -> Result<Vec<Vec<f64>>, String>{
    if a[0].len() != b.len() {
        return Err("Number of columns in 'a' must be equal to number of rows in 'b'".to_string())
    }

    let ar = a.len();
    let bc = b[0].len();

    let mut result = Vec::new();

    let p = b.len();

    for ai in 0..ar {
        let mut result_row = Vec::new();

        for bj in 0..bc {
            let mut value = 0.0;

            for p in 0..p {
                value += a[ai][p] * b[p][bj];
            }

            result_row.push(value);
        }

        result.push(result_row);
    }

    Ok(result)
}

pub fn scalar_mult(matrix: &Vec<Vec<f64>>, s: f64) -> Vec<Vec<f64>> {
    let r = matrix.len();
    let c = matrix[0].len();

    let mut result = Vec::new();
    for i in 0..r {
        let mut row = Vec::new();
        for j in 0..c {
            row.push(matrix[i][j] * s);
        }
        result.push(row);
    }

    return result;
}