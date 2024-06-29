function [G1, G2, G3, a1, a2, b1, b2] =buildSRQ2(Az, B, C, Pz, z, d, blkmode)
% function buildSRQ2 generates coefficient matrices for the SRQ2 minimization 
%	min x'*G1*x / (a1 + b1*x'G3x) + x'*G2*x / (a2 + b2*x'G3x)
% for the backward error analysis of the Rosenbrock system
%	[Az, B]
%	[C, Pz]
%
% INPUT:
%	Az (=A-zI), B, C, Pz (=P(z)): coefficient matrices
%	z: 	eigenvalues
%	d: 	degree of P(z)
%
%	blkmode:	string indicating perturbation type.
%				one block	: A, B, C, P
%				two blocks	: AB, AC, AP, BC, BP, CP
%				three blocks: ABC, ABP, ACP, BCP
%				four blocks	: ABCP
%				NOTE: non-trivial SRQ2 include AP, BC; ABC, ABP, ACP, BCP; ABCP
%
% OUTPUT:
%	Coefficient matrices and scalars; G2 = 0 indicates that only one Rayleigh 
%	quotient presents.
%
%
r = size(Az,1);
n = size(Pz,1);
k = n+r;

% ----------------------
% one-block perturbation
% ----------------------
if strlength(blkmode)==1
	if strfind(blkmode, 'A') 
		U = null([C, Pz]);
		G1 = [Az, B]*U; G1 = G1'*G1; 
		G2 = zeros(k, k);
		G3 = U(1:r,:); G3 = G3'*G3;
		a1 = 0; b1 = 1;
		a2 = 1; b2 = 0;

	elseif strfind(blkmode, 'B') 
		U = null([C, Pz]);
		G1 = [Az, B]*U; G1 = G1'*G1; 
		G2 = zeros(k, k);
		G3 = U(r+1:end,:); G3 = G3'*G3;
		a1 = 0; b1 = 1;
		a2 = 1; b2 = 0;

	elseif strfind(blkmode, 'C') 
		U = null([Az, B]);
		G1 = [C, Pz]*U; G1 = G1'*G1; 
		G2 = zeros(k, k);
		G3 = U(1:r,:); G3 = G3'*G3;
		a1 = 0; b1 = 1;
		a2 = 1; b2 = 0;

	elseif strfind(blkmode, 'P') 
		U = null([Az, B]);
		gamma = norm(z.^[0:1:d])^2;
		G1 = [C, Pz]*U; G1 = G1'*G1/gamma;
		G2 = zeros(k, k);
		G3 = U(r+1:end,:); G3 = G3'*G3;
		a1 = 0; b1 = 1;
		a2 = 1; b2 = 0;

	else 
		disp('Error: 1-block pattern not found');
	end

% ----------------------
% two-block perturbation
% ----------------------
elseif strlength(blkmode)==2
	if strfind(blkmode, 'AB')
		U = null([C, Pz]);
		G1 = [Az, B]*U; G1=G1'*G1;
		G2 = zeros(k, k);
		G3 = eye(r+n);
		a1 = 0; b1 = 1;
		a2 = 1; b2 = 0;

	elseif strfind(blkmode, 'AC')
		U = null([B', Pz']);
		G1 = [Az', C']*U; G1=G1'*G1;
		G2 = zeros(k, k);
		G3 = eye(r+n);
		a1 = 0; b1 = 1;
		a2 = 1; b2 = 0;

	elseif strfind(blkmode, 'AP') %**
		gamma = norm(z.^[0:1:d])^2;
		G1 = [Az, B]; G1 = G1'*G1;
		G2 = [C, Pz]; G2 = G2'*G2;
		G3 = blkdiag(eye(r), zeros(n,n));
		a1 = 0; b1 = 1;
		a2 = gamma; b2 = - gamma;

	elseif strfind(blkmode, 'BC') %**
		G1 = [Az, B]; G1 = G1'*G1;
		G2 = [C, Pz]; G2 = G2'*G2;
		G3 = blkdiag(eye(r), zeros(n,n));
		a1 = 1; b1 = -1;
		a2 = 0; b2 = 1;

	elseif strfind(blkmode, 'BP')
		U = null([Az', C']);
		gamma = norm(z.^[0:1:d])^2;
		Gamma = blkdiag(eye(r), gamma*eye(n));
		G1 = [B', Pz']*U/Gamma; G1 = G1'*G1;
		G2 = zeros(k, k);
		G3 = eye(r+n);
		a1 = 0; b1 = 1;
		a2 = 1; b2 = 0;

	elseif strfind(blkmode, 'CP')
		U = null([Az, B]);
		gamma = norm(z.^[0:1:d])^2;
		Gamma = blkdiag(eye(r), gamma*eye(n));
		G1 = [C, Pz]*U/Gamma; G1 = G1'*G1;
		G2 = zeros(k, k);
		G3 = eye(r+n);
		a1 = 0; b1 = 1;
		a2 = 1; b2 = 0;

	else 
		disp('Error: 2-block pattern not found');
	end

% ------------------------
% three-block perturbation
% ------------------------
elseif strlength(blkmode)==3
	if strfind(blkmode, 'ABC') %**
		G1 = [Az, B]; G1 = G1'*G1;
		G2 = [C, Pz]; G2 = G2'*G2;
		G3 = blkdiag(eye(r), zeros(n,n));
		a1 = 1; b1 = 0;
		a2 = 0; b2 = 1;

	elseif strfind(blkmode, 'ABP') %**
		gamma = norm(z.^[0:1:d])^2;
		G1 = [Az, B]; G1 = G1'*G1;
		G2 = [C, Pz]; G2 = G2'*G2;
		G3 = blkdiag(eye(r), zeros(n,n));
		a1 = 1; b1 = 0;
		a2 = gamma; b2 = -gamma;

	elseif strfind(blkmode, 'ACP') %**
		gamma = norm(z.^[0:1:d])^2;
		G1 = [Az, B]; G1 = G1'*G1;
		G2 = [C, Pz]; G2 = G2'*G2;
		G3 = blkdiag(eye(r), zeros(n,n));
		a1 = 0; b1 = 1;
		a2 = gamma; b2 = 1-gamma;

	elseif strfind(blkmode, 'BCP') %**
		gamma = norm(z.^[0:1:d])^2;
		G1 = [Az, B]; G1 = G1'*G1;
		G2 = [C, Pz]; G2 = G2'*G2;
		G3 = blkdiag(eye(r), zeros(n,n));
		a1 = 1; b1 = -1;
		a2 = gamma; b2 = 1-gamma;

	else 
		disp('Error: 3-block pattern not found');
	end

% ------------------------------
% four-block perturbation (full)
% ------------------------------
elseif strlength(blkmode)==4
	gamma = norm(z.^[0:1:d])^2;
	G1 = [Az, B]; G1 = G1'*G1;
	G2 = [C, Pz]; G2 = G2'*G2;
	G3 = blkdiag(eye(r), gamma*eye(n));
	a1 = 1; b1 = 0;
	a2 = 0; b2 = 1;

else
	disp('Error: Block pattern not found');
end

% symmetrize
G1 = (G1+G1')/2;
G2 = (G2+G2')/2;
G3 = (G3+G3')/2;
return