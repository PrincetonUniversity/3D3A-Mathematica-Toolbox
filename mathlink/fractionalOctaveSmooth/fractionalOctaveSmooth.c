//
//  fractionalOctaveSmooth.c
//  Fractional octave smoothing routine for real valued, 1-D lists in Mathematica.
//

// ==============================================================================
// This file is part of the 3D3A Mathematica Toolbox.
//
// Joseph G. Tylka <josephgt@princeton.edu>
// 3D Audio and Applied Acoustics (3D3A) Laboratory
// Princeton University, Princeton, New Jersey 08544, USA
//
// MIT License
//
// Copyright (c) 2018 Princeton University
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
// ==============================================================================

//  Version History:
//
//  Version 1.2: rearrange input arguments; dynamically allocate char arrays.
//
//  Version 1.1: add 'scale' input argument for smoothing magnitude, power or dB
//      spectra.
//
//  Version 1.0: created from version 2.2 of "realOctaveSmooth" program and
//      version 1.1 of "logCompensatedOctaveSmooth" program.
//
//  Created by Joe Tylka on 9 November 2018.


#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

extern void fractionalOctaveSmooth();
void computeSmoothingMatrix(int FFTLen, double frac, char* method, char* winType, double** M);
void smoothingMatrix_eye(int FFTLen, double** M);
void smoothingMatrix_tylka(int FFTLen, double frac, char* winType, double** M);
void smoothingMatrix_hatz(int FFTLen, double frac, char* winType, double** M);

int main(int argc, char *argv[])
{
    MLMain(argc, argv);
}

void fractionalOctaveSmooth()
{
    /* Inputs */
    double *H;
    double frac;
    const char *methodIn;
    const char *winTypeIn;
    const char *scaleIn;
    int FFTLen = 0;
    char *method;
    char *winType;
    char *scale;
    
    /* Get inputs from Mathematica */
    MLGetReal64List(stdlink, &H, &FFTLen);
    MLGetReal64(stdlink, &frac);
    MLGetString(stdlink, &methodIn);
    MLGetString(stdlink, &winTypeIn);
    MLGetString(stdlink, &scaleIn);
    
    if (frac == 0.0)
    {   // No smoothing case
        MLPutReal64List(stdlink, H, FFTLen);
        MLReleaseReal64List(stdlink, H, FFTLen);
        return;
    }
    
    /* Make strings lower case */
    method = (char*) malloc (sizeof(char) * (strlen(methodIn) + 1));
    for (int ii = 0; ii < strlen(methodIn); ii++)
    {
        method[ii] = tolower(methodIn[ii]);
    }
    method[strlen(methodIn)] = '\0';
    
    winType = (char*) malloc (sizeof(char) * (strlen(winTypeIn) + 1));
    for (int ii = 0; ii < strlen(winTypeIn); ii++)
    {
        winType[ii] = tolower(winTypeIn[ii]);
    }
    winType[strlen(winTypeIn)] = '\0';
    
    scale = (char*) malloc (sizeof(char) * (strlen(scaleIn) + 1));
    for (int ii = 0; ii < strlen(scaleIn); ii++)
    {
        scale[ii] = tolower(scaleIn[ii]);
    }
    scale[strlen(scaleIn)] = '\0';
    
    /* Release input strings */
    MLReleaseString(stdlink, methodIn);
    MLReleaseString(stdlink, winTypeIn);
    MLReleaseString(stdlink, scaleIn);
    
    /* Local variables */
    int specLen = 1 + FFTLen/2;
    double *Hin, *Hout;
    double **M;
    
    /* Outputs */
    double *Hsm;
    
    /* Allocate memory */
    Hin = (double*) malloc (sizeof(double) * FFTLen);
    Hout = (double*) malloc (sizeof(double) * specLen);
    Hsm = (double*) malloc (sizeof(double) * FFTLen);
    M = (double**) malloc (sizeof(double*) * specLen);
    for (int ii = 0; ii < specLen; ii++)
    {
        M[ii] = (double*) malloc (sizeof(double) * FFTLen);
        for (int jj = 0; jj < FFTLen; jj++)
        {
            M[ii][jj] = 0.;
        }
    }
    
    /* Begin real code */
    for (int jj = 0; jj < FFTLen; jj++)
    {   // Re-scale input
        if (strcmp(scale,"raw") == 0)
        {
            Hin[jj] = H[jj];
        }
        else if (strcmp(scale,"power") == 0)
        {
            Hin[jj] = fabs(H[jj]) * fabs(H[jj]);
        }
        else if (strcmp(scale,"db") == 0)
        {
            Hin[jj] = log(fabs(H[jj]));
        }
    }
    
    /* Release input spectrum */
    MLReleaseReal64List(stdlink, H, FFTLen);
    
    /* Compute smoothing matrix */
    computeSmoothingMatrix(FFTLen, frac, method, winType, M);
    
    for (int ii = 0; ii < specLen; ii++)
    {
        Hout[ii] = 0.;
        for (int jj = 0; jj < FFTLen; jj++)
        {   // Perform matrix multiplication
            Hout[ii] = Hout[ii] + (M[ii][jj] * Hin[jj]);
        }
        
        // Re-scale output
        if (strcmp(scale,"raw") == 0)
        {
            Hsm[ii] = Hout[ii];
        }
        else if (strcmp(scale,"power") == 0)
        {
            Hsm[ii] = sqrt(Hout[ii]);
        }
        else if (strcmp(scale,"db") == 0)
        {
            Hsm[ii] = exp(Hout[ii]);
        }
    }
    for (int ii = specLen; ii < FFTLen; ii++)
    {   // Flip for symmetry above Nyquist
        Hsm[ii] = Hsm[FFTLen - ii];
    }
    
    /* Send output to Mathematica */
    MLPutReal64List(stdlink, Hsm, FFTLen);
    
    /* Free memory */
    for (int ii = 0; ii < specLen; ii++)
    {
        free(M[ii]);
    }
    free(M);
    free(Hsm);
    free(Hin);
    free(Hout);
    free(method);
    free(winType);
    free(scale);
    return;
}

void computeSmoothingMatrix(int FFTLen, double frac, char* method, char* winType, double** M)
{
    int specLen = 1 + FFTLen/2;
    
    if (frac == 0.)
    {   // No smoothing applied.
        smoothingMatrix_eye(FFTLen, M);
        return;
    }
    
    if (strcmp(method,"tylka") == 0)
    {
        smoothingMatrix_tylka(FFTLen, frac, winType, M);
    }
    else if ((strcmp(method,"hatz") == 0) || (strcmp(method,"hatziantoniou") == 0))
    {
        smoothingMatrix_hatz(FFTLen, frac, winType, M);
    }
    else
    {
        smoothingMatrix_eye(FFTLen, M);
    }
    return;
}

void smoothingMatrix_eye(int FFTLen, double** M)
{
    int specLen = 1 + FFTLen/2;
    for (int ii = 0; ii < specLen; ii++)
    {
        for (int jj = 0; jj < FFTLen; jj++)
        {
            if (ii == jj)
            {
                M[ii][jj] = 1.;
            }
            else
            {
                M[ii][jj] = 0.;
            }
        }
    }
    return;
}

void smoothingMatrix_tylka(int FFTLen, double frac, char* winType, double** M)
{
    int specLen = 1 + FFTLen/2;
    for (int ii = 0; ii < specLen; ii++)
    {
        for (int jj = 0; jj < FFTLen; jj++)
        {   // Initialize M to zero
            M[ii][jj] = 0.;
        }
    }
    
    /* Local variables */
    int nL = 0;
    int nH = 0;
    double fL = 0.;
    double fH = 0.;
    int winLen = 0;
    double winTemp = 0.;
    double mL = 0.;
    double mH = 0.;
    double winSum = 0.;
    double term1 = 0.;
    double term2 = 0.;
    double term2a = 0.;
    double term2b = 0.;
    
    M[0][0] = 1.; // DC value is never smoothed
    for (int n = 1; n < specLen-1; n++)
    {
        fL = ((double) n) * exp2(-0.5/frac);
        fH = ((double) n) * exp2(0.5/frac);
        nL = (int) floor(fL);
        nH = (int) ceil(fH);
        
        // Upper limit correction
        if (nH > (specLen - 1))
        {
            nH = specLen - 1;
            nL = (int) floor(pow((double) n, 2.) / nH);
        }
        if (fH > (specLen - 1.))
        {
            fH = specLen - 1.;
            fL = pow((double) n, 2.) / fH;
        }
        
        // Window function computation
        winLen = nH - nL + 1;
        if (winLen == 1)
        {
            M[n][n] = 1.;
            continue;
        }
        
        winSum = 0.;
        for (int m = nL; m <= nH; m++)
        {
            mL = fmin(fmax(((double) m) - 0.5, fL), fH);
            mH = fmin(fmax(((double) m) + 0.5, fL), fH);
            if (mL == mH) { continue; }
            
            term1 = log2(mH / ((double) n)) - log2(mL / ((double) n));
            if ((strcmp(winType,"rectangular") == 0) || (strcmp(winType,"rect") == 0))
            {
                M[n][m] = term1;
            }
            else if ((strcmp(winType,"hanning") == 0) || (strcmp(winType,"hann") == 0))
            {
                term2a = sin(2. * M_PI * frac * log2(mH / ((double) n)));
                term2b = sin(2. * M_PI * frac * log2(mL / ((double) n)));
                term2 = term2a - term2b;
                M[n][m] = (term1 / 2.) + (term2 / (4. * M_PI * frac));
            }
            winSum = winSum + M[n][m];
        }
        
        // Normalize window
        for (int m = nL; m <= nH; m++)
        {
            M[n][m] = M[n][m] / winSum;
        }
    }
    M[specLen-1][specLen-1] = 1.; // Nyquist value is never smoothed
    
    return;
}

void smoothingMatrix_hatz(int FFTLen, double frac, char* winType, double** M)
{
    int specLen = 1 + FFTLen/2;
    for (int ii = 0; ii < specLen; ii++)
    {
        for (int jj = 0; jj < FFTLen; jj++)
        {   // Initialize M to zero
            M[ii][jj] = 0.;
        }
    }
    
    /* Local variables */
    int m = 0;
    int mMax = (FFTLen / 2) - 1;
    int winLen = 0;
    double Q = 1. / (exp2(0.5 / frac) - exp2(-0.5 / frac));
    double winSum = 0.;
    
    M[0][0] = 1.; // DC value is never smoothed
    for (int n = 1; n < specLen-1; n++)
    {
        // Window half-width calculation
        m = (int) fmin(floor((0.5 * ((double) n)) / Q), mMax);
        
        // Window function computation
        winLen = (2 * m) + 1;
        if (winLen == 1)
        {
            M[n][n] = 1.;
            continue;
        }
        
        winSum = 0.;
        for (int k = (n - m); k <= (n + m); k++)
        {
            if ((strcmp(winType,"rectangular") == 0) || (strcmp(winType,"rect") == 0))
            {
                M[n][k] = 1.;
            }
            else if ((strcmp(winType,"hanning") == 0) || (strcmp(winType,"hann") == 0))
            {
                M[n][k] = 0.5 * (1. - cos(2. * M_PI * (k - (n - m)) / (winLen - 1.)));
            }
            winSum = winSum + M[n][k];
        }
        
        // Normalize window
        for (int k = (n - m); k <= (n + m); k++)
        {
            M[n][k] = M[n][k] / winSum;
        }
    }
    
    return;
}