function x_soft = softThresholdLinearTransformation(x,tau,linHandle,recHandle)

% performs a soft threshold on the coefficients in a different basis specified by the handle linHandle with inverse recHandle

x_soft = recHandle(soft(linHandle(x),tau));