# MSOT Model-based Reconstruction Toolbox
This is a toolbox to reconstruct initial pressure images in optoacoustic tomography, developed in the Jüstel Lab, Institute for Biological and Medical Imaging, Institute of Computational Biology, Helmholtz Zentrum München, and Chair of Biological Imaging, Technische Universität München.

## Environment setup
* Install the Matlab R2019b or higher with the toolboxes 'Signal Processing', 'Image Processing' and 'Parallel Computing'.
* Set up the Matlab instance for C/C++ Mex file compilation. For windows, it is sufficient to install MinGW. The Mex files in this toolbox require a x64_86 processor architecture and support OpenMP parallelization.
* Download Field II, version 3.30 from http://field-ii.dk/?downloading_2021.html for your operating system and extract the files in 'libraries/Field_II_ver_3_30_[os]/', where '[os]' is either 'windows' or 'linux'.
* Download ShearLab 3D, version 1.1 for Matlab from http://shearlab.math.lmu.de/software and save the extracted folder in 'libraries/ShearLab3Dv11/'.
* Download SpaRSA, version 2.0 from http://www.lx.it.pt/~mtf/SpaRSA and save the extracted folder in 'libraries/SpaRSA_2.0/'. Optional: Reduce the number of required model evaluations in 'SpaRSA.m' (and thus the reconstruction time of SpaRSA-based methods) by removing dispensable initializations and consistency checks at the beginning of the algorithm (lines 375, 383, 401, 427).

## Model-based reconstruction
This toolbox offers the following non-negative reconstruction methods:
* Non-negative reconstruction with L2 regularization (rec_nn_with_L2_reg.m)
* Non-negative reconstruction with Shearlet L1 regularization (rec_nn_with_Shearlet_reg.m)
* Non-negative reconstruction with L1 regularization (rec_nn_with_L1_eye_reg.m)
* Non-negative reconstruction with TV regularization (rec_nn_with_TV_reg.m)

### Example scripts
* Examples for all available reconstruction methods on a basic numerical phantom are provided in `./reconstruction_model_based/reconstruction_overview_script.m`.
* A template on how to apply the available preprocessing steps and reconstruction methods to real data is given in `reconstruction_example_script.m`. In order to load data from real-world measurements, replace the placeholder implementation in the function 'load_raw_data_for_reconstruction.m`' (currently the function returns simulated signals).

### Features
* The toolbox supports to include the total impulse response into the model of an imaging system, as introduced in [1] and [2].
* Imaging systems can be configured by specifying their properties with a new probeId in `reconstruction_model_based/models/Probe.m` and saving the characterized (individual) electrical impulse response (EIR) of the system in `reconstruction_model_based/data/EIRs`. This repository contains an example specification of a made-up imaging probe and an example EIR to illustrate which properties are required.
* During model initialization, the toolbox computes the travel times and the spatial impulse response (SIR) of a system (if SIR correction is enabled) and stores the computed values for reuse in upcoming reconstruction tasks in `reconstruction_model_based/data/travel_time_files` and `reconstruction_model_based/data/SPMRs`, respectively.

## Backprojection reconstruction
This toolbox also implements backprojection reconstruction according to the formula derived in [3] and [4], in the function 'backproject_waveeq'.

## Citation
If you use this toolbox for model-based reconstructions, please cite the papers [1] and [2].

## Contact
* Christoph Dehner (christoph.dehner@tum.de)
* Dominik Jüstel (dominik.juestel@helmholtz-muenchen.de)
* Maximilian Bader (maximilian.s.bader@tum.de)

## References
1. Chowdhury K.B., Prakash J., Karlas A., Jüstel D., and Ntziachristos V. A synthetic total impulse response characterization method for correction of hand-held optoacoustic images. IEEE Transactions on Medical Imaging, 39(10):3218-30, 2020.
2. Chowdhury K.B., Bader M., Dehner C., Jüstel D., and Ntziachristos V. Individual transducer impulse response characterization method to improve image quality of array-based handheld optoacoustic tomography. Optics Letters 46(1):1-4, 2021.
3. Kunyansky, L.A. Explicit inversion formulae for the spherical mean Radon transform. Inverse Problems, 23(1):373-383, 2007.
4. Kuchment, P. and Kunyansky, L. Mathematics of photoacoustic and thermoacoustic tomography. Handbook of Mathematical Methods in Imaging. O. Scherzer. New York, NY, Springer New York: 817-865, 2011.



