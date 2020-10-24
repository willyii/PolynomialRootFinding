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


## Update log:

### 10/03/2020
- Finished GCD for two polynomial
- Pass the test of GCD
- Move some function into util namespace
- Finished square free decomposition
- Finished simple test case of decomposition
- Finished some complicated case for square free decomposition

### 10/05/2020
- Finished the basic operator for ceof and Poly
- Reformated the code, abstract the isolate class and the other class should inherated from this one

### 10/06/2020
- Add Monic function to poly, make it robust to the really close roots
- Finished budan's theorem based root isolation
- Pass the simple tests for the polynomial with no duplicated roots.

### 10/07/2020
- Tried some random case.
- Redfined the newtwon method stop condition
- Redfined the bisection stop condition
- Make the Newtown Method robust to the no root condition
- Combine the decomposition with the root finding.
- Refine the freeSqureDecompo process

### 10/08/2020
- Move the develop process into project field.
- Fix the bug in the roots boundry
- Meet precision problem, try to fix it,

### 10/21/2020
- Finished the Budan method. 
- Finished Test Cases

