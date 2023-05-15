function [Eval Enew dE] = GenMiR_evalE(E, Eval,C, Beta, Pi, Z, X, mu, Nu, alpha, Sigma, Omega, Phi, s, U, V, W, M, N, T, I)

B0 = Beta(C>0);
[row col] = find(C);
K = nnz(C);

Y0 = (Beta*(Z.*repmat(Nu,1,T)));
Y = Y0*Omega;
R = X-repmat(mu,N,1);
V0 = zeros(T,T);
for kk=1:K

	ii = row(kk);
	jj = col(kk);
	zj = Z(jj,:)';
	lj = Nu(jj);
	bij = Beta(ii,jj);
	V0 = V0 + (lj^2)*(2*bij-bij^2)*diag(zj.^2);	
end;

U = (1/N)*diag(sum(R.^2 + 2*R.*Y,1));
V = (1/N)*(Omega.^2+Phi)*V0;
W = (1/N)*(Omega.^2+Phi)*diag(sum(Y0.^2,1));

E_beta = sum(B0.*log((B0)/Pi)+(1-B0).*log((1-B0)./(1-Pi)));
E_gamma = (1/2)*(T*log(s) + (1/s)*trace((Omega-I).^2 + Phi) - sum(log(diag(Phi))));
E_lambda = sum(-log(Nu/alpha)+Nu/alpha);
E_data = (N/2)*sum(log(diag(Sigma))) + (N/2)*trace(Sigma\(U+V+W));

%Compute KL-divergence
Enew = E_beta + E_gamma + E_lambda + E_data;
dE = (E-Enew)/abs(E);
E = Enew;
Eval = [Eval Enew];

if(dE<0)
	warning('Objective function increased!');
end;
