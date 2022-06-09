function x = freq_decomp_transp_model(s,MT,R0T,F,RT,lambda0,lambdaF,lambdaR,Ns)

% - Ns are the lengths of the components of s as returned from decomp_transp_model.
% - F and RT are array of function handles,
% - No regularization/filter isindicated by an empty array.

% preliminaries
idxF = find(~cellfun('isempty',F));
idxR = find(~cellfun('isempty',RT));
NF = length(idxF); % number of non-empty filters
%NR = length(idxR); % number of non-empty regularizers
Nc = length(Ns); % number of components of s

sd = cell(Nc,1);
Nsc = [0;cumsum(Ns)];
for j = 1:Nc
    sd{j} = s((Nsc(j)+1):Nsc(j+1));
end

% precomputations
mTs = cell(1+NF,1);
gd = cell(NF,1);
for j = 1:(NF+1)
    mTs{j} = MT(sd{j});
end
for j = 1:NF
    gd{j} = F{idxF(j)}(mTs{j+1});
end

% reconstruction
x = cell(length(F),1);

for j = 1:length(F)
    x{j} = vec(mTs{1});
    
    % filtering
    idx = find(idxF == j,1);
    if ~isempty(idx)
        x{j} = x{j}+lambdaF(j)^2*vec(mTs{idx+1}-gd{idx});
    end
    
    % regularization
    if ~isempty(R0T)
        x{j} = x{j}+lambda0^2*vec(R0T(sd{1+NF+1}));
    end
    idx = find(idxR == j, 1);
    if ~isempty(idx)
        x{j} = x{j}+lambdaR(j)^2*vec(RT{j}(sd{1+NF+1+idx}));
    end
    
end

x = cell2mat(x);