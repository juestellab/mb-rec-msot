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
     * sirIndices: numPixels x numDetectors int32
     * travelTimes: numPixels x numDetectors int32
     * sir: numSirSteps (time) x numSir (corresponding to indices in sirIndices) double
     * sirTimestepsStartAndEnd: 2 x numSir (corresponding to indices in sirIndices) int32
     * circularInterpolationRadiusInSamples: double
     *
     * Output: wavelengths x numTimesteps x numDetectors  double
     */
    
    // Associate input variables
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","image must be type double.");
    }
    if( !mxIsInt32(prhs[1])) {
        mexErrMsgIdAndTxt("forwardMex:notInt32","sirIndices must be type int32.");
    }
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","travelTimes must be type double.");
    }
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","sir must be type double.");
    }
    if( !mxIsInt32(prhs[4])) {
        mexErrMsgIdAndTxt("forwardMex:notInt32","numTimesteps must be type int32.");
    }
    if( !mxIsInt32(prhs[5])) {
        mexErrMsgIdAndTxt("forwardMex:notInt32","sirTimestepsStartAndEnd must be type int32.");
    }
    if( !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("forwardMex:notDouble","circularInterpolationRadiusInSamples must be type double.");
    }
    
    double *image, *sir;
    image = mxGetDoubles(prhs[0]);
    sir = mxGetDoubles(prhs[3]);
    
    int *sirIndices, *sirStartTimesteps, *sirTimestepsStartAndEnd;
    double *travelTimes;
    sirIndices = mxGetInt32s(prhs[1]);
    travelTimes = mxGetDoubles(prhs[2]);
    sirTimestepsStartAndEnd = mxGetInt32s(prhs[5]);
    
    int numTimesteps = mxGetScalar(prhs[4]);
    
    double circularInterpolationRadiusInSamples;
    circularInterpolationRadiusInSamples = mxGetScalar(prhs[6]);
    
    // Get required dimensions from input data (no consistency checking)
    int numPixels, numWavelengths, numDetectors, numSirSteps, numSirs;
    numPixels = mxGetN(prhs[0]);
    numWavelengths = mxGetM(prhs[0]);
    numDetectors = mxGetN(prhs[1]);
    numSirSteps = mxGetM(prhs[3]);
    numSirs = mxGetN(prhs[3]);
    
    /*mexPrintf(" numPixels: %i \n", numPixels);
    mexPrintf(" numWavelengths: %i \n", numWavelengths);
    mexPrintf(" numDetectors: %i \n", numDetectors);
    mexPrintf(" numSirSteps: %i \n", numSirSteps);
    mexPrintf(" numTimesteps: %i \n", numTimesteps);*/

    // Associate output variable
    double *outputSignal;
    const mwSize outputDims[3] = {(mwSize) numWavelengths, (mwSize) numTimesteps, (mwSize) numDetectors};
    plhs[0] = mxCreateNumericArray(3, outputDims, mxDOUBLE_CLASS, mxREAL);
    outputSignal = mxGetDoubles(plhs[0]);
    
    // General variables for forward projection
    double travelTimeFromTransducerToPixel, factor; 
    int currentRecievingSignalFromPixel, firstRecievingSignalFromPixel, lastRecievingSignalFromPixel;
    int wavelength, detector, x;
    __m256d image_vec, factor_vec, signal_vec;
    
    // Variables specific to SIR correction
    int currentSirIndex, currentSignalIndex, currentSirTimestep;
    double currentSir;
    __m256d currentSir_vec;
	
	mwIndex ind_for_sir; // Quick fix to be able to index large SIR file using unsigned 64 bit int.
    
    #pragma omp parallel for private( \
        travelTimeFromTransducerToPixel,  \
        factor,  \
        currentRecievingSignalFromPixel, \
        firstRecievingSignalFromPixel, \
        lastRecievingSignalFromPixel, \
        wavelength, \
        detector, \
        x, \
        image_vec, \
        factor_vec, \
        signal_vec, \
        currentSirIndex, \
        currentSignalIndex, \
        currentSirTimestep, \
        currentSir, \
        currentSir_vec, \
		ind_for_sir)
    for(detector=0; detector<numDetectors; detector++) {
        
        for(x=0; x<numPixels; x++) {
            //Subtract one to convert from matlab to c++ Indices
            currentSirIndex = sirIndices[detector*numPixels + x]-1;
            
            travelTimeFromTransducerToPixel = travelTimes[detector*numPixels + x];
            firstRecievingSignalFromPixel = (int) ceil(travelTimeFromTransducerToPixel-circularInterpolationRadiusInSamples);
            lastRecievingSignalFromPixel = (int) floor(travelTimeFromTransducerToPixel+circularInterpolationRadiusInSamples);
            
            // iterate SIR of pixel
            for(currentSirTimestep=sirTimestepsStartAndEnd[2*currentSirIndex]-1; currentSirTimestep<=sirTimestepsStartAndEnd[2*currentSirIndex+1]-1; currentSirTimestep++) {
                ind_for_sir = (mwIndex) currentSirIndex*numSirSteps + currentSirTimestep;
				currentSir = sir[ind_for_sir];
                currentSir_vec = _mm256_set1_pd(currentSir);
                
                // iterate Recieving signals of pixel
                for(currentRecievingSignalFromPixel=firstRecievingSignalFromPixel; currentRecievingSignalFromPixel<=lastRecievingSignalFromPixel; currentRecievingSignalFromPixel++) {
                    currentSignalIndex = currentRecievingSignalFromPixel + currentSirTimestep - numSirSteps/2;

                    // Check if current sir-and-pixel-interpolation offset aims at a sample inside the sinogram
                    if(currentSignalIndex >=0 && currentSignalIndex<numTimesteps) {
                            factor = 1.0 - std::abs(0.0 + currentRecievingSignalFromPixel - travelTimeFromTransducerToPixel)/circularInterpolationRadiusInSamples;
                            factor_vec = _mm256_set1_pd(factor);

                        // use SIMD instructions to process 4 wavelengths simultaneously
                        for(wavelength=0; wavelength<(numWavelengths-3); wavelength+=4) {
                            image_vec = _mm256_loadu_pd(&image[x*numWavelengths + wavelength]);

                            signal_vec = _mm256_loadu_pd(&outputSignal[detector*numWavelengths*numTimesteps + currentSignalIndex*numWavelengths + wavelength]);
                            signal_vec = _mm256_add_pd(signal_vec, _mm256_mul_pd(image_vec, _mm256_mul_pd(currentSir_vec, factor_vec)));
                            _mm256_storeu_pd(&outputSignal[detector*numWavelengths*numTimesteps + currentSignalIndex*numWavelengths + wavelength], signal_vec);
                        }

                        // iterate over remaining wavelengths not dividable by 4
                        for(wavelength; wavelength<numWavelengths; wavelength++) {
                            outputSignal[detector*numWavelengths*numTimesteps + currentSignalIndex*numWavelengths + wavelength]
                                += factor * currentSir * image[x*numWavelengths + wavelength];
                        }
                    }
                }
            }
        }
    }
}
