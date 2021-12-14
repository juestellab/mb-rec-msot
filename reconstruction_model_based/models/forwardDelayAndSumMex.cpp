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
     * image: numWavelengths x numPixels double
     * travelTimes: numPixels x numDetectors double
     * numTimesteps: int32
     * circularInterpolationRadiusInSamples: double
     *
     * Output: wavelengths x numTimesteps x numDetectors  double
     */
    
    // Associate input variables
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","image must be type double.");
    }
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","travelTimes must be type double.");
    }
    if( !mxIsInt32(prhs[2])) {
        mexErrMsgIdAndTxt("forwardMex:notInt32","numTimesteps must be type int32.");
    }
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","circularInterpolationRadiusInSamples must be type double.");
    }
    
    double *image;
    image = mxGetDoubles(prhs[0]);
    
    double *travelTimes;
    travelTimes = mxGetDoubles(prhs[1]);
    
    int numTimesteps = mxGetScalar(prhs[2]);
    
    double circularInterpolationRadiusInSamples;
    circularInterpolationRadiusInSamples = mxGetScalar(prhs[3]);

    // Get required dimensions from input data (no consistency checking)
    int numPixels, numWavelengths, numDetectors;
    numPixels = mxGetN(prhs[0]);
    numWavelengths = mxGetM(prhs[0]);
    numDetectors = mxGetN(prhs[1]);
    
    /*mexPrintf(" numPixels: %i \n", numPixels);
    mexPrintf(" numWavelengths: %i \n", numWavelengths);
    mexPrintf(" numDetectors: %i \n", numDetectors);*/

   // Associate output variable
    double *outputSignal;
    const mwSize outputDims[3] = {(mwSize) numWavelengths,(mwSize) numTimesteps,(mwSize) numDetectors};
    plhs[0] = mxCreateNumericArray(3, outputDims, mxDOUBLE_CLASS, mxREAL);
    outputSignal = mxGetDoubles(plhs[0]);
    
    double travelTimeFromTransducerToPixel, factor; 
    int currentRecievingSignalFromPixel, firstRecievingSignalFromPixel, lastRecievingSignalFromPixel;
    int wavelength, detector, x;
    __m256d image_vec, factor_vec, signal_vec;
        
    #pragma omp parallel for private(travelTimeFromTransducerToPixel, factor, currentRecievingSignalFromPixel, firstRecievingSignalFromPixel, lastRecievingSignalFromPixel, detector, x, wavelength, image_vec, factor_vec, signal_vec)
    for(detector=0; detector<numDetectors; detector++) {

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
                    image_vec = _mm256_loadu_pd(&image[x*numWavelengths + wavelength]);

                    // Floor signals
                    signal_vec = _mm256_loadu_pd(&outputSignal[detector*numWavelengths*numTimesteps + currentRecievingSignalFromPixel*numWavelengths + wavelength]);
                    signal_vec = _mm256_add_pd(signal_vec, _mm256_mul_pd(factor_vec, image_vec));
                    _mm256_storeu_pd(&outputSignal[detector*numWavelengths*numTimesteps + currentRecievingSignalFromPixel*numWavelengths + wavelength], signal_vec);
                }

                // iterate over remaining wavelengths not dividable by 4
                for(wavelength; wavelength<numWavelengths; ++wavelength) {
                    outputSignal[detector*numWavelengths*numTimesteps + currentRecievingSignalFromPixel*numWavelengths + wavelength] 
                            += factor*image[x*numWavelengths + wavelength];
                }
            }
        }
    }
}
