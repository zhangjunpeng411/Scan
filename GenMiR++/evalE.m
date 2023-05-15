function [E grad_u] = evalE(u,C,Beta,Pi,Z,X,mu,a,b,m,n,Sigma,G,K,T,I)

a_pos = u(1:K);
b_pos = u(K+1:2*K);
m_pos = u(2*K+1:2*K+T);
n_pos = u(2*K+1+T:end);

H_lambda = sum(gammaH(a_pos,b_pos));

Omega = diag(m_pos.*n_pos);
Phi = diag(m_pos.*(n_pos.^2));

Y0 = (Beta*(Z.*repmat(a_pos.*b_pos,1,T)));
Y = Y0*Omega;
R = X-repmat(mu,G,1);

U1 = (1/G)*diag(sum((R+Y).^2,1));
U2 = (1/G)*Phi*diag(sum(Y0.^2,1));

bv = sum(Beta,1)';
ll = bv.*a_pos.*(b_pos.^2);
U3 = (1/G)*diag(ll'*((Z.^2)*(Omega.^2+Phi)));

bv2 = sum(Beta-Beta.^2,1)';
ll = bv2.*(a_pos.^2).*(b_pos.^2);
U4 = (1/G)*diag(ll'*((Z.^2)*(Omega.^2+Phi)));

E_lambda = -sum((a-1)*Elog(a_pos,b_pos)  - a_pos.*b_pos/b - a*log(b) - gammaln(a));
E_data = (G/2)*sum(log(diag(Sigma))) + (G/2)*trace(Sigma\(U1+U2+U3+U4));

%Compute bound on likelihood
E = H_lambda + E_lambda + E_data;


%Compute gradient
g1 = ((Beta'*(R+Y)).*(repmat(b_pos,1,T).*Z))*diag(Omega/Sigma);
g2 = ((Beta'*Y0).*(repmat(b_pos,1,T).*Z))*diag(Phi/Sigma);
bv = sum(Beta,1)';
ll = bv.*(b_pos.^2);
g3 = 0.5*ll.*((Z.^2)*diag((Omega.^2+Phi)/Sigma));
bv2 = sum(Beta-Beta.^2,1)';
ll = repmat(bv2.*a_pos.*(b_pos.^2),1,T);
g4 = ll.*(Z.^2)*diag((Omega.^2+Phi)/Sigma);

grad_a = g1 + g2 + g3 + g4 - (a-1).*psi(1,a_pos) + b_pos/b - (a_pos-1).*psi(1,a_pos) + 1;
grad_b = zeros(K,1);
grad_m = zeros(T,1);
grad_n = zeros(T,1);

% grad_b = -a_pos/b + ;
% 
% [f, g] = Elog(m_pos,n_pos);
% 
% grad_m = -n_pos/n +  (1 + psi(m_pos) - m_pos.*psi(1,m_pos));
% grad_n = -m_pos/n;

grad_u = [grad_a; grad_b; grad_m; grad_n];

%Move this to main function
% dE = (E-Enew)/abs(E);
% E = Enew;
% Eval = [Eval Enew];
% 
% if(dE<0)
% 	warning('Objective function increased!');
% end;


function f = Elog(a,b)

f = log(b) + psi(0,a);


function H = gammaH(a,b)
H = a + log(b) + gammaln(a) + (1-a).*psi(0,a);

