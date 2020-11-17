# PolynomialRootFinding

This projct is to build a application that can solve polynomial root finding 
problem efficiently. 

The main method applied here is root isolation. Which is continually minize thes range of roots. When the range in small enough, apply Newton Method in that range to find the root with in that range.

This project will apply three methods and compare them together: Budan's threom, Vincent Theorem and Bisection method.

## Usage
Pls make sure you installed [CMake](https://cmake.org/) on your machine. 

Using following command to run the project:

```bash
git clone https://github.com/willyii/PolynomialRootFinding
cd PolynomialRootFinding
make clean
make build
./valid.sh
```

Last command will run the **valid.py** first. Which will generate 10000 random polynomial coefficients and coresponding roots save in **test/validation.test**
(Polynomial with no root will use # to mark). Then it will call this project to compute the roots and valid if it is correct or not.


