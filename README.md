# PolynomialRootFinding

This project is aims to build program to isolate real roots. We applied Budan's
theorem and continued fraction method. We also applied Yun's algorithm to
perform square-free decomposition

## Usage
Pls make sure you installed [CMake](https://cmake.org/) on your machine. 

Using following command to run the project:

```bash
git clone https://github.com/willyii/PolynomialRootFinding
cd PolynomialRootFinding
make clean
make build
```

There are several ways to run this program:

1. 

```bash
./build/Poly
```

This will generate 1000 random polynomials in file **data/random_poly.txt**. And
this program will try to solve them and return the running time.

2.

```bash
./build/Poly valid
```

This will generate 1000 random polynomials and their solutions. This program
will print out the exact answer and the answer of this program. User can check
if it can return correct answer.

3.

```bash
./build/Poly -path_to_test_file valid
```

The test polynomials in test_file should be in form like a\*x\^3+b\*x\^2+c. It
will return the answer of these polynomials


