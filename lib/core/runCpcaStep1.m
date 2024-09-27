function OutStruct = runCpcaStep1(Z, G, H)

% X - rows are observations (Points), columns are features (variables)
% X should be zscored along the columns
% G - Predictor variable matrix along 1st dim -rows
% H - Predictor variable matrix along 2nd dim - cols

N = size(Z,1);
if ~all(size(G) == 1) && all(size(H) == 1)
    C = pinv(G'*G)*G'*Z; % Linear regression
    S = G*C; E = Z - S;
    OutStruct.C = C;
    OutStruct.B = NaN;
elseif all(size(G) == 1) && ~all(size(H) == 1)
    B = Z*H*pinv(H'*H); % Linear regression
    S = B*H'; E = Z - S;
    OutStruct.B = B;
    OutStruct.C = NaN;
elseif all(size(G) == 1) && all(size(H) == 1)
    S = Z; E = Z - S;
    OutStruct.C = NaN;
    OutStruct.B = NaN;
else
    return;
end

temp = corrcoef(zscore(S));
temp(isnan(temp)) = 0;
[~, OutStruct.EigValGH] = eig(temp);
OutStruct.EigValGH = sort(diag(OutStruct.EigValGH),'descend');
[OutStruct.UGH, OutStruct.DGH, OutStruct.VGH] = svd(S,'econ'); % SVD for GC or BH'

OutStruct.loadings_VD_sqrtN = OutStruct.VGH*diag(diag(OutStruct.DGH))/sqrt(N-1);
OutStruct.Z = Z;
OutStruct.G = G;
OutStruct.H = H;
OutStruct.E = E;
OutStruct.S = S;
OutStruct.subIdx = NaN;
OutStruct.nComp = NaN;

end