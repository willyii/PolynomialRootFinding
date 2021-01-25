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
static const double kEPSILON = 1e-9;

// Maxium possible degree of polynomials
static const int kMAXDEGREE = 7;

// Minimum range
static const double kMINRANGE = 1e-6;

// Store file path for generated random polynomials
static const char *kRANDOM_FILE = "src/random_poly.txt";

#endif // POLY_PARAM_H
