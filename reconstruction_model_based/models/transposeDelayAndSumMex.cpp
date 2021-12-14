/*********************************************************************
* Keep in mind:
* <> Use 0-based indexing as always in C or C++
* <> Indexing is column-based as in Matlab (not row-based as in C)
* <> Use linear indexing.  [x*dimy+y] instead of [x][y]
*
********************************************************************/
#include "mex.h"
#include <omp.h>
#include "immintrin.h" // for AVX 
#include <cmath>
/*
 * Compiling instructions are available in ths script ./compile_mex_files.m
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Expected inputs:
     * sigMat: numWavelengths x numTimesteps x numDetectors double
     * travelTimes: numPixels x numDetectors double
     * circularInterpolationRadiusInSamples: double
     * 
     * Output: numWavelengths x numPixels double
     */
    
    // Associate input variables
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","sigMat must be type double.");
    }
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","travelTimes must be type double.");
    }
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","circularInterpolationRadiusInSamples must be type double.");
    }
    
    double *sigMat;
    sigMat = mxGetDoubles(prhs[0]);
    
    double *travelTimes;
    travelTimes = mxGetDoubles(prhs[1]);
    
    double circularInterpolationRadiusInSamples;
    circularInterpolationRadiusInSamples = mxGetScalar(prhs[2]);
    
    // Get required dimensions from input data (no consistency checking)
    int numTimesteps, numPixels, numWavelengths, numDetectors;
    numTimesteps = mxGetDimensions(prhs[0])[1];
    numPixels = mxGetM(prhs[1]);
    numDetectors = mxGetN(prhs[1]);
    
    numWavelengths = mxGetM(prhs[0]);
    
    /*mexPrintf(" numPixels: %i \n", numPixels);
    mexPrintf(" numDetectors: %i \n", numDetectors);
    mexPrintf(" numWavelengths: %i \n", numWavelengths);
    mexPrintf(" numTimesteps: %i \n", numTimesteps);*/
    
   // Associate output variable
    double *outputImage;
    plhs[0] = mxCreateDoubleMatrix(numWavelengths, numPixels, mxREAL);
    outputImage = mxGetDoubles(plhs[0]);
    int wavelength, detector, x;
    double travelTimeFromTransducerToPixel, factor;
    int currentRecievingSignalFromPixel, firstRecievingSignalFromPixel, lastRecievingSignalFromPixel;
    __m256d image_vec, signal_vec, factor_vec;
        
    for(detector=0; detector<numDetectors; detector++) {
         
        #pragma omp parallel for private(travelTimeFromTransducerToPixel, factor, currentRecievingSignalFromPixel, firstRecievingSignalFromPixel, lastRecievingSignalFromPixel, x, wavelength, image_vec, signal_vec, factor_vec)
        for(x=0; x<numPixels; x++) {
            travelTimeFromTransducerToPixel = travelTimes[detector*numPixels + x];
            firstRecievingSignalFromPixel = (int) ceil(travelTimeFromTransducerToPixel-circularInterpolationRadiusInSamples);
            lastRecievingSignalFromPixel = (int) floor(travelTimeFromTransducerToPixel+circularInterpolationRadiusInSamples);
            
            for(currentRecievingSignalFromPixel=firstRecievingSignalFromPixel; currentRecievingSignalFromPixel<=lastRecievingSignalFromPixel; currentRecievingSignalFromPixel++) {
                // Calculate scaling factor for interpolation & stream to SIMD vector
                factor = 1.0 - std::abs(0.0 + currentRecievingSignalFromPixel - travelTimeFromTransducerToPixel)/circularInterpolationRadiusInSamples;
                factor_vec = _mm256_set1_pd(factor);
               
                // use SIMD instructions to process 4 wavelengths simultaneously
                for(wavelength=0; wavelength<(numWavelengths-3); wavelength+=4) {
                    image_vec = _mm256_loadu_pd(&outputImage[x*numWavelengths + wavelength]);

                    signal_vec = _mm256_loadu_pd(&sigMat[detector*numWavelengths*numTimesteps + currentRecievingSignalFromPixel*numWavelengths + wavelength]);
                    image_vec = _mm256_add_pd(image_vec, _mm256_mul_pd(factor_vec,signal_vec));

                    _mm256_storeu_pd(&outputImage[x*numWavelengths + wavelength], image_vec);
                }

                // iterate over remaining wavelengths not dividable by 4
                for(wavelength; wavelength<numWavelengths; ++wavelength) {
                    outputImage[x*numWavelengths + wavelength] = outputImage[x*numWavelengths + wavelength] 
                            + factor*sigMat[detector*numWavelengths*numTimesteps + currentRecievingSignalFromPixel*numWavelengths + wavelength];
                }
            }
        }
    }
}
