function [y] = PhiForShearletReg(x,shearletTrafo)
% apply the Phi function on all wavelengths
    Phi_ow = @(inp) L1normLinearTransformation(inp, shearletTrafo);
	y = zeros(size(x,3), 1);
			
    parfor i = 1:size(x,3)
		y(i) = Phi_ow(x(:,:,i));
    end
    
    y = sum(y);
end
