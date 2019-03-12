function MahalMat = fastMahalDist(Y)
    %Y is data where columns are observations
    %MahalMat is the squared Mahalanobis matrix or:
    %pdist2(Y',Y','mahalanobis').^2=fastMahalDist(Y)
    covMat =cov(Y);
    Y=Y';
    hx =chol(covMat);
    Ynew = hx'\Y;
    G=Ynew'*Ynew;
    dim=size(G,1);
    MahalMat=diag(G)*ones(1,dim)-2*G+ones(dim,1)*diag(G)';
end