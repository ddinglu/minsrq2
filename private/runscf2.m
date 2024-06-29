function [v0, fv0, OBJFX, RRESD, RESD, neig] = runscf2(A, B, C, alpha1, alpha2, gamma1, gamma2, v00, tol)
% function runscf2 applies SCF iteration to solve NEPv from SRQ2 minimization 
% 	min (x'Ax)/(alpha1+gamma1*x'Cx) + (x'Bx)/(alpha2+gamma2*x'Cx) s.t. x'x=1 
% where A, B, C are positive semidefinite definite.
%
% Adaptive level-shifts (with replacement of solution) is used for stablizing
% convergence of SCF. 
%
% INPUT: 
% 	- {A, B, C, alpha1, gamma1, alpha2, gamma2}: coeffs for SRQ2 minimization 
%	- v00: starting vector for SCF (randn by default)
%	- tol: relative residual tolerance (1.0E-14 by defalut) 
%
% OUTPUT: 
%	- v0: 	solution  
% 	- fv0: 	f(v0)
%   - OBJFX: 	convergence history of objective values
% 	- RRESD: 	convergence history relative residual norms
% 	- RRES: 	convergence history of residual norms
% 	- neig:	number of calls to eig
%

n = size(A,1);
maxit = 300;

if nargin < 8
   	v00 = randn(n,1)+1i*randn(n,1); 
	v00 = v00/norm(v00);
end

if nargin < 9, tol = 1.0E-14; end

% define coefficient functions 
rx = @(x) real([x'*A*x, x'*B*x, x'*C*x]);
objfy = @(y) y(1)/(alpha1+gamma1*y(3)) + y(2)/(alpha2+gamma2*y(3));
objf = @(x) objfy(rx(x)); 
Hy = @(y) A/(alpha1+gamma1*y(3)) + B/(alpha2+gamma2*y(3)) - (gamma1*y(1)/(alpha1+gamma1*y(3)).^2 + gamma2*y(2)/(alpha2+gamma2*y(3)).^2)*C;
Hx = @(x) Hy(rx(x));

RRESD = []; % relative residual norm
RESD = []; 	% residual norm
OBJFX = []; % objective values

v0 = v00;

% main loop: SCF with level-shift
sigma = 0;
gamma = 2;
v0new = true;
neig = 0;
for ii = 1:maxit
	if v0new % v0 is updated in the previous iteration
		v0 = v0/norm(v0);
		HH = Hx(v0); 
		resdv = HH*v0 - v0*real(v0'*HH*v0);
		resd = norm(resdv);
	   	nh1 = norm(HH,1) + 1;
		rresd = resd/nh1;
		objfx = objf(v0);
		RESD = [RESD, resd*2]; % *2 so that RESD = gradnorm
		RRESD = [RRESD, rresd]; 
		OBJFX = [OBJFX, objfx]; 
	end

	if rresd < tol, break; end % residual tolerance

	% SCF
	HHs = HH -  sigma*v0*v0';
	[vv,ee] = eig(HHs); 
	[ee,idx] = sort(real(diag(ee)), 'ascend');
	neig = neig + 1;
	vv = vv(:,idx(1)); vv = vv/norm(vv);

	% backtracking line search for level-shift
	if sigma == 0, delta = real(ee(idx(2))-ee(idx(1))); end;
	if objf(vv) > objfx && norm(Hx(vv)*vv - vv*real(vv'*Hx(vv)*vv)) > resd 
	   	sigma = gamma*delta;
		gamma = gamma*2;
		v0new = false;
	else
		v0 = vv;
		v0new = true;
		sigma = 0;
		gamma = 2;
   	end
end

fv0 = objf(v0);

% END
end