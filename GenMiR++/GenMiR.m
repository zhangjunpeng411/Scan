% /filename GenMiR.m 
% /description Implementation of the GenMiR++ model & Variational Bayes algorithm for learning miRNA-target interactions from sequence and expression data
% /vars   X           NxT matrix of mRNA expression data
%         Z           MxT matrix of miRNA expression data
%         C           NxM matrix of putative regulatory relations bewteen miRNAs and mRNAs
%         Beta        MxN matrix of variational target selection parameters
%         Nu          Mx1 variational parameter vector for regulatory weights distribution
%		  Omega,Phi   TxT mean and variance diagonal parameter matrices for variational tissue scaling coefficients distribution
%		  a,s         Hyperparameters for regulatory weights and tissue-scaling coefficients
%         mu,Sigma    Tx1 background transcription rate and TxT diagonal covariance matrix
%         Pi          Base target selection probability  
%         glbl, mlbl  mRNA and miRNA labels, respectively
%
% /params display  Boolean flag for displaying data
%         num_iter Maximum number of iterations (default 100)
%         tol      Stopping criterion  (default 1e-5)
%
% /version 3.0, January 2007
% /author Jim Huang, PSI Group, University of Toronto
% $Id$

function Score = GenMiR(Z, X, C)

M = size(Z, 1);
N = size(X, 1);
T = size(Z, 2);

num_iter = 100;
tol = 1e-5;
display=0;    

h = waitbar(0,'Loading data...');
v = get(h, 'Position');
v(4) = 100;
set(h, 'Position', v, 'Name', 'GenMiR++');

[row col] = find(C);
K = length(row);

%Initialize parameters
I = eye(T);
s = 5e-2;
mu = mean(X,1);
alpha = 5e-4;
Nu = alpha*ones(M,1);
Pi = 0.5;
Beta = Pi*C;
Omega = I;
Phi = s*I;
Beta_mat = [];

%Compute D(q||p)
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
Sigma = diag(diag(U+V+W));

[Eval E dE] = GenMiR_evalE(1e30,[],C, Beta, Pi, Z, X, mu, Nu, alpha, Sigma, Omega, Phi, s, U, V, W, M, N, T, I);
t = 0;
n_iters = 1;

%Start GenMiR
waitbar(1,h,'Data loaded. Iterating GenMiR++ algorithm...');
disp([num2str(nnz(C)) ' predicted miRNA-target interactions involving ' num2str(nnz(sum(C,1))) ' miRNAs and ' num2str(nnz(sum(C,2))) ' target mRNAs']);
disp('Iterating GenMiR++ algorithm...');

tic
while (dE > tol & t < num_iter)
		
	%Variational Bayes E-step: Infer miRNA-target interactions	
	t=t+1;
 	waitbar(t/(num_iter),h,{'Variational Bayes E-step'; 'Inferring miRNA-target interactions from data'; ['D(q||p) = ' num2str(E)]; ['\Delta = ' num2str(dE)]});
	Beta = GenMiR_VBEStep(C, Beta, Pi, Z, X, mu, Nu, Sigma, Omega, Phi, M, N, T);

	%Variational Bayes M-step
 	waitbar(t/(num_iter),h,{'Variational Bayes M-step'; 'Marginalizing over parameters'; ['D(q||p) = ' num2str(E)] ; ['\Delta = ' num2str(dE)]});				
	[Pi Nu alpha Omega Phi mu Sigma U V W] = GenMiR_VBMStep(C, Beta, Pi, Z, X, mu, Nu, alpha, Sigma, Omega, Phi, s, M, N, T, I);
			
	%Evaluate bound on likelihood
	[Eval E dE] = GenMiR_evalE(E,Eval,C, Beta, Pi, Z, X, mu, Nu, alpha, Sigma, Omega, Phi, s, U, V, W, M, N, T, I);

	if(display)
		plot(Eval, 'LineWidth', 2);
		xlabel('Number of iterations','FontSize',24); ylabel('L(Q)','FontSize',24); drawnow;
		handle = gca;
		set(handle,'FontSize',24);
		box off;
	end;

end;

%Score miRNA-target interactions
waitbar(1,h,'Scoring putative interactions...');
disp('Scoring putative interactions...');
Score = zeros(nnz(C),1);
%Score = num2cell(Score);
for kk=1:K

	ii = row(kk);
	jj = col(kk);
	
	bij = Beta(ii,jj);
    %Score{kk,1} = mlbl{jj,1};
    %Score{kk,2} = glbl{ii,1};
	Score(kk) = log10((bij+eps)/(1-bij+eps));
end;

waitbar(1,h,['Algorithm terminated in ' num2str(t-1) ' iterations.']);
disp(['Algorithm terminated in ' num2str(t-1) ' iterations.']);
toc

set(h, 'Name', 'GenMiR++ - Done. Press any key to exit.')
waitbar(1,h,'Press any key to exit.');
disp('Press any key to exit.');
%pause
close(h);
