//
//  HatzOctaveSmooth.c
//  Fractional Octave Smoothing routine for real valued, 1-D lists in Mathematica.
//
//  Version 1.0:
//
//  Created by Joe Tylka on 7/28/14.
//  Copyright (c) 2014 Joe Tylka. All rights reserved.
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

void HatzSmooth(double SmoothFrac, int WinType)
{
    // Inputs
    int HLen = 0;
    double *H;
    
    // Get Inputs from Mathematica
    MLGetReal64List(stdlink, &H, &HLen);
    
    // No Smoothing Case
    if(SmoothFrac == 0.0)
    {
        MLPutReal64List(stdlink, H, HLen);
        MLReleaseReal64List(stdlink, data, HLen);
        return;
    }
    
    // Outputs
    double *Hsm = malloc(sizeof(double) * HLen);
    Hsm[0] = H[0];
    Hsm[HLen - 1] = H[HLen - 1];
    
    // Local Variables
    int FullLen = 2 * (HLen - 1);
    int M = (FullLen / 2) - 1;
    double Q = 1 / (exp2(0.5 / SmoothFrac) - exp2(-0.5 / SmoothFrac));
    double b = 0.5; // Only needed for Hanning window
    double *HFull = malloc(sizeof(double) * FullLen);
    int ii = 0;
    double iidouble = 0.0;
    int jj = 0;
    double jjdouble = 0.0;
    int m = 0;
    double mdouble = 0.0;
    int indx = 0;
    double win = 0.0;
    double partsum = 0.0;
    
    HFull[0] = H[0];
    for(ii = 1; ii < HLen - 1; ii = ii + 1)
    {
        HFull[ii] = H[ii];
        HFull[FullLen - ii] = H[ii];
    }
    HFull[HLen - 1] = H[HLen - 1];
    
    // Begin Real Code
    for(ii = 1; ii < HLen - 1; ii = ii + 1)
    {
        iidouble = (double) ii;
        m = (int) floor((0.5 * iidouble) / Q);
        if(m < 1)
        {
            m = 1;
        }
        if(m > M)
        {
            m = M;
        }
        mdouble = (double) m;
        
        partsum = 0.0;
        for(jj = -m; jj <= m; jj = jj + 1)
        {
            switch(WinType)
            {
                case 0: // Rectangular Window
                    win = 1.0 / (2 * mdouble + 1);
                    break;
                default: // Hanning Window
                    jjdouble = (double) jj;
                    win = (b - (b - 1) * cos(jjdouble * M_PI / mdouble)) / (2 * b * mdouble + 1);
                    break;
            }
            if(ii + jj < 0)
            {
                indx = ((ii + jj) % FullLen + FullLen) % FullLen;
                partsum = partsum + HFull[indx] * win;
            }
            else
            {
                partsum = partsum + HFull[ii + jj] * win;
            }
        }
        Hsm[ii] = partsum;
    }
    
    // Send Output to Mathematica
    MLPutReal64List(stdlink, Hsm, HLen);
    
    // Free Memory
    free(Hsm);
    MLReleaseReal64List(stdlink, H, HLen);
    return;
}