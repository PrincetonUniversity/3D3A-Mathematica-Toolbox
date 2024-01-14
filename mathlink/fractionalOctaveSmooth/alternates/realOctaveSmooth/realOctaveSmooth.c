//
//  fractionalOctaveSmooth.c
//  Fractional Octave Smoothing routine for real valued, 1-D lists in Mathematica.
//
//  Version 2.2: Removed window creating functions, instead compute the value of the point in the window before each computation. Change uplpos computation to only calculate those values which are less than Len.
//
//  Version 2.1: Changed calloc to malloc. Removed "out" pointer. Removed n = 0.0 in createHanning(). Added exact window normalization instead of computing it in createHanning().
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
    int winLen = 0;
    int firstCap = 0;
    double n = 0.5;
    double iiDouble = 0.0;
    double LenDouble = 0.0;
    double partsum = 0.0;
    double win = 0.0;
    double winLenDouble = 0.0;
    double twopi = 2.0 * M_PI;
    double indexRatio = exp2(1.0/fractionDenominator);
        
    /* Get Inputs from Mathematica */
    MLGetReal64List(stdlink, &data, &Len);
    
    // No Smoothing Case
    if(fractionDenominator == 0.0)
    {
        MLPutReal64List(stdlink, data, Len);
        MLReleaseReal64List(stdlink, data, Len);
        return;
    }
    
    /* Outputs */
    double *outvec = malloc(sizeof(double) * Len);
    
    /* Begin Real Code */
    LenDouble = (double) Len;
    firstCap = floor(LenDouble / indexRatio);
    for(ii = 0; ii < Len; ii = ii + 1)
    {
        // Compute 1/N octave bounds
        lowlpos = floor(iiDouble / indexRatio);
        if(ii < firstCap)
        {
            uplpos = ceil(iiDouble * indexRatio);
        }
        else
        {
            uplpos = Len - 1;
        }
        
        winLen = uplpos - lowlpos + 1;
        winLenDouble = (double) winLen;
        
        // Compute sum of window x data
        switch(windowType)
        {
            case 0: // Rectangular Window
                win = 1.0/winLenDouble;
                for(jj = 0; jj < winLen; jj = jj + 1)
                {
                    partsum = partsum + data[lowlpos + jj] * win;
                }
                break;
            default: // Hanning Window
                for(jj = 0; jj < winLen; jj = jj + 1)
                {
                    win = (1 - cos(twopi * n / winLenDouble))/winLenDouble;
                    partsum = partsum + data[lowlpos + jj] * win;
                    n = n + 1.0;
                }
                n = 0.5;
                break;
        }
        outvec[ii] = partsum;
        partsum = 0.0;
        iiDouble = iiDouble + 1.0;
    }
    
    /* Send Output to Mathematica */
    MLPutReal64List(stdlink, outvec, Len);
    
    /* Free Memory */
    free(outvec);
    MLReleaseReal64List(stdlink, data, Len);
    return;
}