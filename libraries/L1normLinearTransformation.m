function L1norm = L1normLinearTransformation(x,linHandle)

% computes the L1 norm in a different basis specified by the handle linHandle

L1norm = sum(abs(vec(linHandle(x))));