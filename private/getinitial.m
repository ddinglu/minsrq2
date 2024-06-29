function [v0] = getinitial(A,B,C,alpha1,alpha2,gamma1,gamma2)
% function getinitial generate initial vectors for SCF for SRQ2
% minimization,using the eigenvector corresponding to the larger Rayleigh
% quotient.

n = size(A,1);

% find minimizer v1 of RQ1
G1 = A + alpha1*eye(n) + gamma1*C; 
[V1, E1] = eig(alpha1*eye(n) + gamma1*C, G1);
[~,idx1] = max(real(diag(E1)));
v1 = V1(:,idx1(1)); v1 = v1/norm(v1);
f1 = real(v1'*A*v1)/real(alpha1 + gamma1*v1'*C*v1);

% find minimizer v2 of RQ22
G2 = B + alpha2*eye(n) + gamma2*C; 
[V2, E2] = eig(alpha2*eye(n) + gamma2*C, G2);
[~,idx2] = max(real(diag(E2)));
v2 = V2(:,idx2(1)); v2 = v2/norm(v2);
f2 = real(v2'*A*v2)/real(alpha2 + gamma2*v2'*C*v2);

% select v0 from v1 and v2
if f1>f2
	v0 = v2;
else
	v0 = v1;
end

% END
return;