function Res = SelectCompLoadings(cpca_struct, params)

Z_dim = size(cpca_struct.Z,2);
S = cpca_struct.S;
G = params.G_orig;

for ii = 1:Z_dim

    Rej_feat = ii;
    Sel_feat = setdiff(1:Z_dim, Rej_feat);

    N = size(S,1);
    G_CPCA = runCpcaStep1(S(:,Sel_feat), S(:,Rej_feat), 1);
    Res(ii).ResidueS = G_CPCA.E;
    Res(ii).ResidueScores = sqrt(N-1) .* Res(ii).ResidueS * cpca_struct.VGH(Sel_feat,:) / cpca_struct.DGH;

    if params.varimaxFlag
        Res(ii).ResidueScoresRot = Res(ii).ResidueScores(:,1:cpca_struct.nComp)*cpca_struct.T;
    else
        Res(ii).ResidueScoresRot = Res(ii).ResidueScores(:,1:cpca_struct.nComp);
    end

    Res(ii).ResiduePCorr = corr(G, Res(ii).ResidueScoresRot);
    Res(ii).meanPCorr_rel = sum(abs(Res(ii).ResiduePCorr).*cpca_struct.RelPredLoadMask, 1)./sum(cpca_struct.RelPredLoadMask, 1);
    Res(ii).meanPCorr_rel_corrected = sum(abs(Res(ii).ResiduePCorr).*cpca_struct.RelPredLoadMask_corrected, 1)./sum(cpca_struct.RelPredLoadMask_corrected, 1);

end