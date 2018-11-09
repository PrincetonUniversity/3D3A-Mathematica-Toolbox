//
//  makeMinimumPhase.c
//  Minimum-phase conversion routine for real valued, 1-D lists in Mathematica.
//
//  Version 1.0: adapted from BacchDSP.c in LibBacchAD_v1.65.
//
//  Created by Joe Tylka on 19 March 2016.
//  Copyright (c) 2016 Joe Tylka. All rights reserved.
//

#define EPS 0.000000000000000152655665885959
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>
#include "mathlink.h"

extern void makeMinimumPhase(int p);
void rCeps(double* x, double* yHat, int IRLen);
void rotate_IR(double* vec, int len, int delay);

int main(int argc, char *argv[])
{
    MLMain(argc, argv);
}

void makeMinimumPhase(int p)
{
    /* Inputs */
    int IRLen = 0;
    double *x;
    
    /* Local Variables */
    int ii = 0;
    int TFLen = 0;
    int bigSize = 0;
    int size1 = 0;
    double *x_temp, *y_temp, *y_big, *y_hat;
    fftw_complex *X_temp, *X_zero;
    fftw_plan myfftPLAN1, myifftPLAN1;
    
    /* Get Inputs from Mathematica */
    MLGetReal64List(stdlink, &x, &IRLen);
    
    /* Outputs */
    double *outvec = malloc(sizeof(double) * IRLen);
    
    /* Begin Real Code */
    TFLen = IRLen/2 + 1;
    bigSize = (int) pow(2.0, ceil(log2((double) IRLen)) + (double) p);
    size1 = (bigSize - IRLen)/2;
    
    // Allocate memory for FFTW computations
    x_temp = (double*) fftw_malloc(sizeof(double) * IRLen);
    y_temp = (double*) fftw_malloc(sizeof(double) * IRLen);
    X_temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TFLen);
    X_zero = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TFLen);
    
    // Create FFT Plans
    myfftPLAN1 = fftw_plan_dft_r2c_1d(IRLen, x_temp, X_temp, FFTW_ESTIMATE);
    myifftPLAN1 = fftw_plan_dft_c2r_1d(IRLen, X_zero, y_temp, FFTW_ESTIMATE);
    
    // Copy over input IR
    for(ii = 0; ii < IRLen; ii++)
    {
        x_temp[ii] = x[ii];
    }
    
    fftw_execute(myfftPLAN1); // X_temp = FFT(x_temp);
    
    // Compute zero-phase frequency response
    for (ii = 0; ii < TFLen; ii++)
    {
        X_zero[ii] = cabs(X_temp[ii]); // Warning: implicit cast to complex
    }
    
    fftw_execute(myifftPLAN1); // y_temp = IFFT(X_zero);
    
    // Center-align IR in array
    rotate_IR(y_temp, IRLen, ((int) floor(IRLen/2)));
    
    // Allocate bigSize IRs here
    // DEBUG: try with just one array
    y_big = malloc(sizeof(double) * bigSize);
    y_hat = malloc(sizeof(double) * bigSize);
    
    // Zero-pad to prevent Hilbert artifacts
    for (ii = 0; ii < size1; ii++)
    {
        y_big[ii] = 0.0;
    }
    for (ii = size1; ii < size1 + IRLen; ii++)
    {
        y_big[ii] = y_temp[ii - size1] / ((double) IRLen);
    }
    for (ii = size1 + IRLen; ii < bigSize; ii++)
    {
        y_big[ii] = 0.0;
    }
    
    // Compute minimum-phase signal with same real cepstrum
    rCeps(y_big, y_hat, bigSize);
    
    // Pass back outside
    for (ii = 0; ii < IRLen; ii++)
    {
        outvec[ii] = y_hat[ii];
    }
    
    /* Send Output to Mathematica */
    MLPutReal64List(stdlink, outvec, IRLen);
    
    /* Free Memory */
    free(outvec);
    free(y_big);
    free(y_hat);
    fftw_destroy_plan(myfftPLAN1);
    fftw_destroy_plan(myifftPLAN1);
    fftw_free(X_temp);
    fftw_free(X_zero);
    fftw_free(x_temp);
    fftw_free(y_temp);
    void fftw_cleanup(void);
    MLReleaseReal64List(stdlink, x, IRLen);
    return;
}

// Computes and returns yHat, the unique minimum-phase sequence with the same real cepstrum as x
void rCeps(double* x, double* yHat, int IRLen)
{
    // Declare local variables
    int ii = 0;
    int TFLen = IRLen/2 + 1;
    double temp;
    double *x_temp, *x_ceps, *z_ceps, *y_ceps;
    fftw_complex *X_temp, *X_ceps, *Z_ceps, *Y_ceps;
    fftw_plan cepsFFTPlan1, cepsFFTPlan2, cepsFFTPlan3, cepsFFTPlan4;
    
    // Allocate memory for first FFT
    x_temp = (double*) fftw_malloc(sizeof(double) * IRLen);
    X_temp = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TFLen);
    cepsFFTPlan1 = fftw_plan_dft_r2c_1d(IRLen, x_temp, X_temp, FFTW_ESTIMATE);
    
    // Copy over input IR
    for(ii = 0; ii < IRLen; ii++)
    {
        x_temp[ii] = x[ii];
    }
    
    fftw_execute(cepsFFTPlan1); // X_temp = FFT(x_temp);
    
    // Allocate memory for second FFT
    x_ceps = (double*) fftw_malloc(sizeof(double) * IRLen);
    X_ceps = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TFLen);
    cepsFFTPlan2 = fftw_plan_dft_c2r_1d(IRLen, X_ceps, x_ceps, FFTW_ESTIMATE);
    
    // Compute logarithm of magnitude response
    for (ii = 0; ii < TFLen; ii++)
    {
        temp = cabs(X_temp[ii]);
        if(temp == 0.0)
        {
            fprintf(stderr, "WARNING: rCeps: zero found in FFT of signal at index %i \n", ii);
            temp = EPS; // IS THIS LEGAL??? no.
        }
        X_ceps[ii] = log(temp); // Warning: implicit cast to complex
    }
    
    // apparently this IFFT is outputting infinite values! that's the problem
    // TODO: fix it?
    
    // Clean up first FFT
    fftw_destroy_plan(cepsFFTPlan1);
    fftw_free(x_temp);
    fftw_free(X_temp);
    
    fftw_execute(cepsFFTPlan2); // x_ceps = IFFT(X_ceps) = IFFT(log|FFT(x)|);
    
    // Allocate memory for third FFT
    z_ceps = (double*) fftw_malloc(sizeof(double) * IRLen);
    Z_ceps = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TFLen);
    cepsFFTPlan3 = fftw_plan_dft_r2c_1d(IRLen, z_ceps, Z_ceps, FFTW_ESTIMATE);
    
    // Take Hilbert transform
    z_ceps[0] = 1.0 * creal(x_ceps[0] / ((double) IRLen));
    for (ii = 1; ii < IRLen/2; ii++)
    {
        z_ceps[ii] = 2.0 * creal(x_ceps[ii] / ((double) IRLen));
    }
    z_ceps[IRLen/2] = 1.0 * creal(x_ceps[IRLen/2] / ((double) IRLen));
    for (ii = IRLen/2 + 1; ii < IRLen; ii++)
    {
        z_ceps[ii] = 0.0 * creal(x_ceps[ii] / ((double) IRLen));
    }
    
    // Clean up second FFT
    fftw_destroy_plan(cepsFFTPlan2);
    fftw_free(x_ceps);
    fftw_free(X_ceps);
    
    fftw_execute(cepsFFTPlan3); // Z_ceps = FFT(z_ceps);
    
    // Allocate memory for fourth FFT
    y_ceps = (double*) fftw_malloc(sizeof(double) * IRLen);
    Y_ceps = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TFLen);
    cepsFFTPlan4 = fftw_plan_dft_c2r_1d(IRLen, Y_ceps, y_ceps, FFTW_ESTIMATE);
    
    // Compute complex exponential
    for (ii = 0; ii < TFLen; ii++)
    {
        Y_ceps[ii] = cexp(Z_ceps[ii]);
    }
    
    // Clean up third FFT
    fftw_destroy_plan(cepsFFTPlan3);
    fftw_free(z_ceps);
    fftw_free(Z_ceps);
    
    fftw_execute(cepsFFTPlan4); // y_ceps = IFFT(Y_ceps);
    
    // Normalize and pass outside
    for (ii = 0; ii < IRLen; ii++)
    {
        yHat[ii] = creal(y_ceps[ii] / ((double) IRLen));
    }
    
    // Clean up fourth FFT
    fftw_destroy_plan(cepsFFTPlan4);
    fftw_free(y_ceps);
    fftw_free(Y_ceps);
    
    // Finish clean up
    void fftw_cleanup(void);
}

void rotate_IR(double* vec, int len, int delay)
{
    int jj = 0;
    int N = delay%len; // N will have the same sign as delay
    double tempVec[len];
    
    if(N == 0)
    {
        // No rotation needed
        return;
    }
    else
    {
        // Add len to N to make N positive
        if(N < 0)
        {
            N = N + len;
        }
        
        // Perform rotation
        for(jj = 0; jj < (len - N); jj++)
        {
            tempVec[jj + N] = vec[jj];
        }
        for(jj = len - N; jj < len; jj++)
        {
            tempVec[jj - (len - N)] = vec[jj];
        }
        
        // Pass rotated vector outside
        for(jj = 0; jj < len; jj++)
        {
            vec[jj] = tempVec[jj];
        }
    }
} // Verified interpolated IRs and resulting .bac files 3/18/14 - passed