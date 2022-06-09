function [s,Ns] = freq_decomp_model(x,M,R0,F,R,lambda0,lambdaF,lambdaR)

% - x needs to be of size pixelNum*length(F).
% - R0 is the regularizer function handle for the reconstruction,
% - F and R are cell arrays of function handles of identical length,
%   containing the filters and regularizers.
% - No regularization/filter isindicated by an empty array.

% preliminaries
idxF = find(~cellfun('isempty',F));
idxR = find(~cellfun('isempty',R));
NF = length(idxF); % number of non-empty filters
NR = length(idxR); % number of non-empty regularizers
Nx = length(x)/length(F);

% precomputations
xd = reshape(x,Nx,[]);
md = M(xd);
fd = zeros(size(md));
for j = 1:NF
    fd(:,:,j) = F{idxF(j)}(md(:,:,idxF(j)));
end

% reconstruction
s = cell(1+NF+1+NR,1);
s{1} = vec(sum(md,3));

% filtering
for j = 1:NF
    s{1+j} = lambdaF(idxF(j))^2*vec(md(:,:,idxF(j))-fd(:,:,j));
end

% regularization
if ~isempty(R0)
    s{1+NF+1} = lambda0^2*vec(R0(sum(xd,2)));
    for j = 1:NR
        s{1+NF+1+j} = lambdaR(idxR(j))^2*vec(R{idxR(j)}(xd(:,idxR(j))));
    end
else
    for j = 1:NR
        s(1+NF+j) = lambdaR(idxR(j))^2*vec(R{idxR(j)}(xd(:,idxR(j))));
    end
end

Ns = cellfun('length',s);
s = cell2mat(s);