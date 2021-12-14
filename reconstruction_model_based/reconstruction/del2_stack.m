function d2x = del2_stack(x)

% compute 2D-Laplacian on stack of 2D images

d2x = zeros(size(x));
for j = 1:size(x,3)
    d2x(:,:,j) = del2(x(:,:,j));
end

end