function varStruct = getVarianceTable(cpca_struct)
Z = cpca_struct.Z;
S = cpca_struct.S;
E = cpca_struct.E;

Scores = cpca_struct.scores_NU_GH;
Loadings = cpca_struct.loadings_VD_sqrtN_GH;

varStruct.varZ = estVar(Z); varStruct.perVarZ = 100;
varStruct.varS = estVar(S); varStruct.perVarS = varStruct.varS/varStruct.varZ*100;
varStruct.varE = estVar(E); varStruct.perVarE = varStruct.varE/varStruct.varZ*100;

varStruct.varEachZ = estEachVar(Z);
varStruct.varEachS = estEachVar(S);
varStruct.perVarEachS = varStruct.varEachS./varStruct.varEachZ*100;

varStruct.varScores = estVar(Scores*Loadings'); varStruct.perVarScores = varStruct.varScores/varStruct.varZ*100;

for ii = 1:size(Scores,2)
    varStruct.varScoresCompWise(ii) = estVar(Scores(:,ii)*Loadings(:,ii)');
    varStruct.perVarScoresCompWise(ii) = varStruct.varScoresCompWise(ii)./varStruct.varZ*100;
    varStruct.perVarScoresCompWiseBetween(ii) = varStruct.varScoresCompWise(ii)./varStruct.varS*100;
end
end

function varOut = estVar(X)
[~,D,~] = svd(X,"econ");
varOut = sum(D(:).^2);
end

function varOut = estEachVar(X)
varOut = NaN(size(X,2),1);
for ii = 1:size(X,2)
    [~,D,~] = svd(X(:,ii),"econ");
    varOut(ii) = sum(D(:).^2);
end
end

