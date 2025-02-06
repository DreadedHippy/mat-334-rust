# Matrix Calculator

This is a Rust project with 0 external dependencies that solves some math problems, including:

- Finding the solutions to linear equations using Gaussian Elimination and Cramer's rule of matrices.

## Getting Started

### Prerequisites

- Rust programming language installed on your machine. You can download it from [here](https://www.rust-lang.org/tools/install).

### Running the Project

To run this project, enter the command:

```sh
cargo run -- <solution method> <size of nxn matrix>
```

Replace `<solution method>` with `gaussian`, `cramer`, or one of the other [commands](#command-aliases) depending on the method you want to use, and `<size of nxn matrix>` with the size of the matrix you are working with.

### Example

```sh
cargo run -- gaussian 3
```

This command will run the Gaussian Elimination method on a 3x3 matrix.

## Note
For equations with no solution, the program crashes

### Command Aliases

| Operation  | Command Aliases |
| ------------- | ------------- |
| Cramer's Rule | `cramer`, `c`  |
| Gaussian Elimination  | `gaussian`, `g`, `gauss` |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Rust documentation and community for their support and resources.
