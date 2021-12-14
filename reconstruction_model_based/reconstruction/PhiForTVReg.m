function [y] = PhiForTVReg(x)
% apply the Phi function on all wavelengths
	y = zeros(size(x,3), 1);
			
    parfor i = 1:size(x,3)
		y(i) = ITV(x(:,:,i));
    end
    
    y = sum(y);
end
