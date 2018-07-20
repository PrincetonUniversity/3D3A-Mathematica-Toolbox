//
//  logCompensatedOctaveSmooth.c
//  Fractional Octave Smoothing routine using Logarithmically-Compensated Weights for real valued, 1-D lists in Mathematica.
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
//  Version 1.1: added wTotal variable to normalize each window.
//
//  Version 1.0: created from version 2.2 of "realOctaveSmooth" program.
//
//  Created by Joe Tylka on 18 March 2016.


#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

extern void logCompensatedOctaveSmooth(double fractionDenominator, int windowType);

int main(int argc, char *argv[])
{
    MLMain(argc, argv);
}

void logCompensatedOctaveSmooth(double fractionDenominator, int windowType)
{
    /* Inputs */
    int fullLen = 0;
    double *data;
    
    /* Local Variables */
    int n = 1;
    int m = 0;
    int nL = 0;
    int nH = 0;
    int halfLen = 0;
    double w = 0.0;
    double fL = 0.0;
    double fH = 0.0;
    double mL = 0.0;
    double mH = 0.0;
    double logmL = 0.0;
    double logmH = 0.0;
    double wTotal = 0.0;
    double nDouble = 1.0;
    double mDouble = 0.0;
    double partsum = 0.0;
    double twopi = 2.0 * M_PI;
    double twopiFrac = twopi * fractionDenominator;
    double indexRatio = exp2(1.0/(2.0*fractionDenominator));
    
    /* Get Inputs from Mathematica */
    MLGetReal64List(stdlink, &data, &fullLen);
    
    // No Smoothing Case
    if(fractionDenominator == 0.0)
    {
        MLPutReal64List(stdlink, data, fullLen);
        MLReleaseReal64List(stdlink, data, fullLen);
        return;
    }
    
    /* Outputs */
    double *outvec = malloc(sizeof(double) * fullLen);
    
    /* Begin Real Code */
    halfLen = fullLen / 2;
    outvec[0] = data[0]; // DC value is never smoothed
    for(n = 1; n < halfLen; n = n + 1)
    {
        // Compute 1/N octave bounds
        fL = nDouble / indexRatio;
        nL = (int) floor(fL);
        fH = nDouble * indexRatio;
        nH = (int) ceil(fH);
        
        // Upper limit correction
        if(nH > halfLen)
        {
            nH = halfLen;
            nL = floor(pow(nDouble, 2.0) / nH);
        }
        if(fH > halfLen)
        {
            fH = halfLen;
            fL = pow(nDouble, 2.0) / fH;
        }
        
        // Compute sum of window x data
        for(m = nL; m <= nH; m = m + 1)
        {
            mDouble = (double) m;
            mL = mDouble - 0.5; // Lower limit of integral
            if(mL < fL) { mL = fL; }
            if(mL > fH) { mL = fH; }
            mH = mDouble + 0.5; // Upper limit of integral
            if(mH < fL) { mH = fL; }
            if(mH > fH) { mH = fH; }
            
            w = 0.0;
            if(mL != mH)
            {
                logmL = log2(mL / nDouble);
                logmH = log2(mH / nDouble);
                switch(windowType)
                {
                    case 0: // Rectangular Window
                        w = logmH - logmL;
                        break;
                    case 1: // Hanning Window
                        w = (logmH - logmL + (sin(twopiFrac * logmH) - sin(twopiFrac * logmL)) / twopiFrac) / 2.0;
                        break;
                }
            }
            wTotal = wTotal + w;
            partsum = partsum + data[m] * w;
        }
        outvec[n] = partsum / wTotal;
        wTotal = 0.0;
        partsum = 0.0;
        nDouble = nDouble + 1.0;
    }
    outvec[halfLen] = data[halfLen]; // Nyquist value is never smoothed
    for(n = halfLen + 1; n < fullLen; n = n + 1)
    {
        outvec[n] = outvec[fullLen - n];
    }
    
    /* Send Output to Mathematica */
    MLPutReal64List(stdlink, outvec, fullLen);
    
    /* Free Memory */
    free(outvec);
    MLReleaseReal64List(stdlink, data, fullLen);
    return;
}