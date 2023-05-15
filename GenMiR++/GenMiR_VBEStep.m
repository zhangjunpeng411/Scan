function Beta = GenMiR_VBEStep(C, Beta, Pi, Z, X, mu, Nu, Sigma, Omega, Phi, M, N, T)

[row col] = find(C);
K = nnz(C);
R = X-repmat(mu,N,1);

for kk=1:K

	ii = row(kk);
	jj = col(kk);
	zj = Z(jj,:)';
	Lk = Nu;
	Lk(jj)=[];
	lj = Nu(jj);

	if(~isempty(Lk))
		Zk = Z;
		Zk(jj,:)=[];
		Bk = Beta(ii,:);
		Bk(:,jj)=[];
		yi = 0.5*(Omega.^2+Phi)*(Bk*(Zk.*repmat(Lk,1,T)))';
		a = (Pi/(1-Pi))*exp(-lj*(zj'*(Sigma\(Omega*R(ii,:)' + yi + lj*(Omega.^2+Phi)*zj))));
	else
		a = (Pi/(1-Pi))*exp(-lj*(zj'*(Sigma\(Omega*R(ii,:)' + lj*(Omega.^2+Phi)*zj))));
	end;

	bij = a/(1+a);
	Beta(ii,jj) = bij;
end;
