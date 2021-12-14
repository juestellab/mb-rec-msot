% Add entire directory to MATLAB path
addpath(genpath('.'));

% Enable warnings
warning('on');

% Recompile mex-files
run(fullfile('.', 'reconstruction_model_based','models','compile_mex_files.m'));

if ispc
    run(fullfile('.', 'libraries','fast_marching','functions','compile_c_files_function.m'));
    run(fullfile('.', 'libraries','fast_marching','shortestpath','compile_c_files_shortestpath.m'));
else
    warning(['Mex files for fast-marching library are only available on windows. '...
    'Fast-marching algorithm will be anyways used to calculate the travel times for dual-sos-models without SIR correction but will be very slow.'])
end


