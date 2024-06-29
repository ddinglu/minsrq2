function [x, xcost, info] = runrtr(A, B, C, alpha1, alpha2, gamma1, gamma2, v00, tol)
% function runrtr applies Riemanian Trust Regions method to solve the
% optimization problem 
%	min (x'Ax)/(alpha1+gamma1*x'Cx) + (x'Bx)/(alpha2+gamma2*x'Cx) s.t. x'x=1
% where A, B, C are positive semidefinite definite.
%
% INPUT: 
% 	- {A, B, C, alpha1, gamma1, alpha2, gamma2}: coeffs for SRQ2 minimization 
%	- v00: starting vector for RTR (randn by default)
%	- tol: residual tolerance (1.0E-10 by defalut) 
%
% OUTPUT: 
%	- x: 		solution  
% 	- xcost: 	f(x)
%   - info: 	convergence history 
%

n = size(A,1);
isrealprob = isreal(A) * isreal(B) * isreal(C);
if nargin < 9, tol = 1.0E-10; end % residual tolerance
if nargin < 8
   	v00 = randn(n,1)+1i*randn(n,1)*(1-realprob); 
	v00 = v00/norm(v00);
end

% create the problem structure
if isrealprob
	manifold = spherefactory(n);
else
	manifold = spherecomplexfactory(n);
end
problem.M = manifold;
rx = @(x) real([x'*A*x, x'*B*x, x'*C*x]);
objfy = @(y) y(1)/(alpha1+gamma1*y(3)) + y(2)/(alpha2+gamma2*y(3));
problem.cost = @(x) objfy(rx(x)); 
problem.egrad = @(x) 2*A*x/real(alpha1+gamma1*x'*C*x) + 2*B*x/real(alpha2+gamma2*x'*C*x ) - 2*(gamma1*real(x'*A*x)/real(alpha1+gamma1*x'*C*x).^2 + gamma2*real(x'*B*x)/real(alpha2+gamma2*x'*C*x).^2)*C*x; % Euclidean gradient

% %Hessian: Not used for complex x.
%problem.ehess = @(x,u) 2*A*u/real(alpha1+gamma1*x'*C*x) + 2*B*u/real(alpha2+gamma2*x'*C*x ) -2*A*x/real(alpha1+gamma1*x'*C*x).^2*(2*real(gamma1*x'*C*u)) - 2*B*x/real(alpha2+gamma2*x'*C*x ).^2*(2*real(gamma2*x'*C*u)) - 2*(gamma1*C*u)*real(x'*A*x)/real(alpha1+gamma1*x'*C*x).^2 - 4*(gamma1*C*x)*real(x'*A*u)/real(alpha1+gamma1*x'*C*x).^2 + 8*(gamma1*C*x)*real(x'*A*x)/real(alpha1+gamma1*x'*C*x).^3*real(gamma1*x'*C*u) - 2*(gamma2*C*u)*real(x'*B*x)/real(alpha2+gamma2*x'*C*x).^2 - 4*(gamma2*C*x)*real(x'*B*u)/real(alpha2+gamma2*x'*C*x).^2 + 8*(gamma2*C*x)*real(x'*B*x)/real(alpha2+gamma2*x'*C*x).^3*real(gamma2*x'*C*u);

% %Numerically check gradient consistency (optional).
%checkgradient(problem);
%checkhessian(problem);
 
options.tolgradnorm = tol;
options.verbosity = 0;
[x, xcost, info] = trustregions(problem, v00, options);

%END
end