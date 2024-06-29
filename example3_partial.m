% Example 3: Backward error with partial block perturbation
% Demonstrate efficiency of NEPv approach against manifold
% optimization for computing eigenvalue backward error of Rosenbrock system
%	[ A - z*I, 	B 	]
%	[ C,		P(z)]
% with P(z)=P0+z*P1.

close all; clear;
warning('off', 'all');
set(0,'defaultTextInterpreter','latex'); 
rng(0); 		% for reproducing results
mtol = 1.0E-10;	% tolerance for relative residual norm 

% -------------
% testing mode 
% -------------
testmode = 1; 
% 1: 	r=10, n = [200, 400, 600, 800]; 
% 2:  	r = [200, 400, 600, 800], n=10
% 3:  	r = [200, 400, 600, 800], n=r
if testmode == 1
	% Case I: r=10, n = [200, 400, 600, 800];
    psize =  [ 	10, 200; 
    			10, 400; 
    			10, 600; 
    			10, 800];
elseif testmode == 2
    % Case II:  r = [200, 400, 600, 800], n=10
    psize =  [ 	200, 10; 
    			400, 10; 
    			600, 10; 
    			800, 10];
else
    % Case III:  r = [200, 400, 600, 800], n=r
    psize =  [ 	200, 200; 
    			400, 400; 
    			600, 600; 
    			800, 800];
end

TT1 = []; TT2 =[];
RR1 = []; RR2 =[];
FF1 = []; FF2 =[];
IITS1 = []; IITS2 = []; 
NNE = [];

Probs = ["AP", "BC", "ABC", "ABP", "ACP", "BCP", "ABCP"];
% --------------------------------
% main loop over each testing case 
% --------------------------------
for jj = 1:length(psize)
	r = psize(jj,1); n = psize(jj,2);
	% generate random data matrices for the Rosenbrock systems
	A = randn(r,r)+1i*randn(r,r); 
	B = randn(r,n)+1i*randn(r,n); 
	C = randn(n,r)+1i*randn(n,r); 
	P0 = randn(n,n)+1i*randn(n,n); 
	P1 = randn(n,n)+1i*randn(n,n); 
	d = 1; 		% degree of P

	% sample approximate eigenvalue z with controlled 'backward error'
	epsln = 1.0E-1;
	dA = A.*(1+randn(r,r)*epsln); 
	dB = B.*(1+randn(r,n)*epsln);
	dC = C.*(1+randn(n,r)*epsln);
	dP0 = P0.*(1 + randn(n,n)*epsln); 
	dP1 = P1.*(1 + randn(n,n)*epsln); 
	ee = eig([dA, dB; dC, dP0], blkdiag(eye(r),dP1));
	z = ee(end); 	

	Az = A-eye(r)*z;
	Pz = P0 - P1*z;

	% computation history
	T1 = []; T2 =[];
   	R1 = []; R2 =[];
	F1 = []; F2 =[];
	ITS1 = []; ITS2 = []; 
	NE = [];

	% loop over each individual testing case with different size of the problem 
	for prob = Probs
		% generate SRQ2 problem
		[G1, G2, G3, a1, a2, b1, b2] =buildSRQ2(Az, B, C, Pz, z, d, prob);

		% generate inital v00 using minimizers of RQ
		v00 = getinitial(G1,G2,G3,a1,a2,b1,b2);

		% run SCF
		tic;
		[v0, fv0, OBJFX, RRESD, RESD, neig] = runscf2(G1,G2,G3,a1,a2,b1,b2,v00,mtol);
		t1 = toc;
		its1 = length(OBJFX);

		% run RTR
        tol = RESD(end); % tolerance conversion: manopt uses residual norm as stopping criteria
		tic;
		[v0, fv0, info] = runrtr(G1,G2,G3,a1,a2,b1,b2,v00,tol);
		t2 = toc;
		RESDRTR = arrayfun(@(info) info.gradnorm, info);
		OBJRTR = arrayfun(@(info) info.cost, info);
		its2 = length(OBJRTR);

		% print results 
		disp([prob, num2str(r), num2str(n)])
		disp('Its: SCF (eig) / RTR');
		disp([its1, neig, its2]);
		disp('Timing: SCF / RTR / Speedup');
		disp([t1, t2, t2/t1]);
		disp('Resd: SCF / RTR');
		disp([RESD(end), RESDRTR(end)]);

		disp('Objective Value: SCF / RTR');
		disp([OBJFX(end), OBJRTR(end)]);

		% computation history 
		ITS1 = [ITS1, its1]; ITS2 = [ITS2, its2]; NE = [NE, neig]; 
		T1 = [T1, t1]; T2 = [T2, t2]; 
		R1 = [R1, RESD(end)];
	    R2 = [R2, RESDRTR(end)];
		F1 = [F1, OBJFX(end)];
	    F2 = [F2, OBJRTR(end)];
	end

	TT1 = [TT1; T1]; TT2 = [TT2; T2];
	RR1 = [RR1; R1]; RR2 = [RR2; R2];
	FF1 = [FF1; F1]; FF2 = [FF2; F2];
	IITS1 = [IITS1; ITS1]; IITS2 = [IITS2; ITS2];
	NNE = [NNE; NE]; 
end

% print results
disp('Timing: RTR/SCF')
disp([psize, round(TT2./TT1)])

disp('Fval: RTR/SCF')
disp([psize, round(FF2./FF1)])

disp('Resd: SCF:RTR')
disp(RR1)
disp(RR2)

disp('Its: SCF:RTR')
disp(IITS1)
disp(IITS2)

%END
return