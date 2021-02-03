// ----------------------------------------------------------------------------
//
// FILENAME: param.h
//
// DESCRIPTION:
//    This file defines some parameters that will be used for the program
//
// AUTHOR: Xinlong Yi
//
// ----------------------------------------------------------------------------

#ifndef POLY_PARAM_H
#define POLY_PARAM_H

// Regard number as 0 when it's smaller than kEPSILON
double kEPSILON = 1e-4;

// Maxium possible degree of polynomials
static const int kMAXDEGREE = 7;

// Minimum range
static const double kMINRANGE = 1e-13;

// Store file path for generated random polynomials
static const char *kRANDOM_FILE = "data/random_poly.txt";
static const char *kRANDOM_FILE_SOL = "data/random_poly_sol.txt";

// Number of polynomials generated
static const int kNUMPOLY = 1000;

#endif // POLY_PARAM_H
