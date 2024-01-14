//
//  fractionalOctaveSmooth.c
//  Fractional Octave Smoothing routine for real valued, 1-D lists in Mathematica.
//
//  Version 2.0: Changed lowlposv and uplposv to lowlpos and uplpos - no longer length Len vectors, but single integers. Freed allocated memory for outvec and win.
//
//  Version 1.0: fractionDenominator must be of type double. windowType = 0 corresponds to a rectangular window, while any other value corresponds to a Hanning window. Other window types are easily added by adding a new function and defining a type number.
//
//  Created by Joe Tylka on 8/20/13.
//  Copyright (c) 2013 Joe Tylka. All rights reserved.
//

#include <mathlink.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

void realOctaveSmooth(double fractionDenominator, int windowType);
void createHanning(int N, double* win);
void createRectangular(int N, double* win);

int main(int argc, char *argv[])
{
    MLMain(argc, argv);
}

void realOctaveSmooth(double fractionDenominator, int windowType)
{
    /* Inputs */
    int Len = 0;
    double *data;
    
    /* Local Variables */
    int ii = 0;
    int jj = 0;
    int lowlpos = 0;
    int uplpos = 0;
    int tempLen = 0;
    double temp = 0.0;
    double partsum = 0.0;
    double indexRatio = exp2(1.0/fractionDenominator);
        
    /* Get Inputs from Mathematica */
    MLGetReal64List(stdlink, &data, &Len);
    
    if(fractionDenominator == 0.0)
    {
        MLPutReal64List(stdlink, data, Len);
        MLReleaseReal64List(stdlink, data, Len);
        return;
    }
    
    /* More Local Variables */
    double *win = calloc(Len, sizeof(double));
    
    /* Outputs */
    double *outvec = calloc(Len, sizeof(double));
    double *out;
    
    /* Begin Real Code */
    for(ii = 0; ii < Len; ii = ii + 1)
    {
        // Compute 1/N octave bounds
        lowlpos = floor(ii / indexRatio);
        uplpos = ceil(ii * indexRatio);
        if(uplpos >= Len)
        {
            uplpos = Len - 1;
        }
        
        tempLen = uplpos - lowlpos + 1;
        
        // Create Window
        switch(windowType)
        {
            case 0:
                createRectangular(tempLen,win);
                break;
            default:
                createHanning(tempLen,win);
                break;
        }
        
        // Compute sum of window x data
        for(jj = 0; jj < tempLen; jj = jj + 1)
        {
            partsum = partsum + data[lowlpos + jj] * win[jj];
        }
        outvec[ii] = partsum;
        partsum = 0.0;
    }
    
    /* Format Output */
    out = outvec;
    
    /* Send Output to Mathematica */
    MLPutReal64List(stdlink, outvec, Len);
    
    /* Free Memory */
    free(win);
    free(outvec);
    MLReleaseReal64List(stdlink, data, Len);
    return;
}

void createHanning(int N, double* win)
{
    int ii = 0;
    double n = 0.5;
    double Ntemp = (double) N;
    double partsum = 0.0;
    double twopi = 2.0 * M_PI;
    
    for(ii = 0; ii < N; ii = ii + 1)
    {
        win[ii] = 0.5*(1 - cos(twopi * n / Ntemp));
        partsum = partsum + win[ii];
        n = n + 1.0;
    }
    n = 0.0;
    
    for(ii = 0; ii < N; ii = ii + 1)
    {
        win[ii] = win[ii]/partsum;
    }
    
}

void createRectangular(int N, double* win)
{
    int ii = 0;
    double Ntemp = (double) N;
    
    for(ii = 0; ii < N; ii = ii + 1)
    {
        win[ii] = 1.0/Ntemp;
    }
}