function [f_post, f_sample, prefilt] = shift_inv_proj(f,L,phi)

%% Projection of a function onto a shift-invariant space
%
% [f_post, f_sample, prefilt] = shift_inv_proj(f,L,phi)
%
% Input: f - function to project (array of any dimension)
%        L - vector of integer lattice constants, length should equal the
%            dimension of f; should divide the number of pixels in the
%            respective dimension
%        phi - pixel shape model; same dimensions as f
%
% Output: f_post - projection of f onto the shift-invariant space defined
%                  by phi and L, evaluated at the same points as f
%         f_sample - coefficients of the shift-invariant function f_post
%         prefilt - prefilter for the shift-invariant space defined by phi
%                   and L
%
% Example: % shift-invariant projection of a disc
%          v = linspace(-1,1,400);
%          [x,y] = meshgrid(v);
%          f = double(sqrt(x.^2+y.^2) <= 0.5);
%          phi = exp(-(x.^2 + y.^2)/2/0.02^2);
%          L = [5,5];
%          [f_post, f_coeffs, prefilt] = shift_inv_proj(f,L,phi);
%          figure;
%          subplot(1,2,1); imagesc(v,v,f); axis equal;
%          title('Original function');
%          subplot(1,2,2); imagesc(v,v,f_post); axis equal;
%          title('Shift-invariant projection');
%          % note that the pixel model is too coarse to resolve edges
%
% Restrictions: If the integers in L don't divide the length of phi in
%               the respective dimension, you will get artifacts at the central axes of
%               the image due to overlap.

%% preliminaries

% handle exceptional 1D case
d = length(L);
phi_flipped = false;
f_flipped = false;
if d == 1
    L = [L 1];
    if(size(phi,1) == 1)
        phi = phi';
        phi_flipped = true;
    end
    if(size(f,1) ~= size(phi,1))
        f = f';
        f_flipped = true;
    end
end

% check, if integers in L divides divides length of phi in respective dimension
assert(sum(mod(size(phi),L))==0,'shift_inv_proj.m: Lattice constants L don''t divide length of phi in respective dimension.');

% define sampling operator
S = {};
for j = 1:d
    S = [S(:)' {1:L(j):size(f,j)}];
end

%% shift-invariant projection

% computation of the prefilter
phi_hat = fftn(phi);
autocorr_of_phi = ifftn(abs(phi_hat).^2,'symmetric');               % autocorr of phi (for circular conv.) 
autocorr_of_phi_on_sampling_grid = zeros(size(f));                  
autocorr_of_phi_on_sampling_grid(S{:}) = autocorr_of_phi(S{:});     % periodization of autocorr on reciprocal lattice (conv with comb of dirac impulses)
prefilt = conj(phi_hat ./ fftn(autocorr_of_phi_on_sampling_grid));  % conj, because phi is real & ifft(conj(phi_hat(omega))=conj(phi(-t))


% prefiltering
f_pre = ifftn(fftn(f).*prefilt,'symmetric');
% sampling
f_sample = zeros(size(f));
f_sample(S{:}) = f_pre(S{:});
% postfiltering
f_post = ifftn(fftn(f_sample).*phi_hat,'symmetric');

% extract coefficients
f_sample = f_sample(S{:});

% handle exceptional 1D case
if d ==1
    if phi_flipped
        prefilt = prefilt';
    end
    if f_flipped
        f_sample = f_sample';
        f_post = f_post';
    end
end