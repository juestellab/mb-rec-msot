% Compiling is tested using MinGW 6.3 on Windows 10 with Matlab2019b
mex -R2018a CXXFLAGS="$CXXFLAGS -fopenmp -mavx" LDFLAGS="$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3 -fwrapv -DNDEBUG" forwardMex.cpp
mex -R2018a CXXFLAGS="$CXXFLAGS -fopenmp -mavx" LDFLAGS="$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3 -fwrapv -DNDEBUG" transposeMex.cpp

mex -R2018a CXXFLAGS="$CXXFLAGS -fopenmp -mavx" LDFLAGS="$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3 -fwrapv -DNDEBUG" forwardDelayAndSumMex.cpp
mex -R2018a CXXFLAGS="$CXXFLAGS -fopenmp -mavx" LDFLAGS="$LDFLAGS -fopenmp" CXXOPTIMFLAGS="-O3 -fwrapv -DNDEBUG" transposeDelayAndSumMex.cpp
