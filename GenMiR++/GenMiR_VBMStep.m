function [Pi Nu alpha Omega Phi mu Sigma U V W] = GenMiR_VBMStep(C, Beta, Pi, Z, X, mu, Nu, alpha, Sigma, Omega, Phi, s, M, N, T, I)

                                               
V0 = zeros(T,T);
R = X-repmat(mu,N,1);

for jj=1:M
	zj = Z(jj,:)';
	Zk = Z;
	Zk(jj,:)=[];
	Bk = Beta;
	Bk(:,jj)=[];
	Lk = Nu;
	Lk(jj)=[];
	Bj = Beta(:,jj);
	Yj = (Bk*(Zk.*repmat(Lk,1,T)));
	Ej = (Sigma\(Omega.^2+Phi))*zj;

	aj = 2*sum(Bj)*zj'*Ej;
	bj = 1/alpha + (zj'*(Sigma\Omega)*R' + 0.5*Ej'*Yj')*Bj;
	cj = -1;

	Nu(jj) = max([(-bj + sqrt(bj^2-4*aj*cj))./(2*aj) (-bj - sqrt(bj^2-4*aj*cj))./(2*aj)]);
	V0 = V0 + (Nu(jj)^2)*diag(sum(2*Bj-Bj.^2).*(zj.^2));
end;

Y0 = (Beta*(Z.*repmat(Nu,1,T)));
A = Sigma\(V0 + diag(sum(Y0.^2,1)));
B = -Sigma\diag(sum(R.*Y0,1));
Omega = (B+(1/s)*I)/(A+(1/s)*I);
Phi = (A+(1/s)*I)\I;
                                               

Y0 = (Beta*(Z.*repmat(Nu,1,T)));
Y = Y0*Omega;
mu = mean(X+Y,1);

R = X-repmat(mu,N,1);
U = (1/N)*diag(sum(R.^2 + 2*R.*Y,1));
V = (1/N)*(Omega.^2+Phi)*V0;
W = (1/N)*(Omega.^2+Phi)*diag(sum(Y0.^2,1));
Sigma = diag(diag(U+V+W));

Pi = sum(sum(Beta))/nnz(C);
alpha = mean(Nu);
