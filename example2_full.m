% Example 2: Backward error with full perturbation
% Demonstrate efficiency of NEPv approach against manifold
% optimization for computing eigenvalue backward error of Rosenbrock system
%	[ A - z*I, 	B 	]
%	[ C,		P(z)]
% with P(z)=P0+z*P1.

close all; clear;
warning('off','all');
set(0,'defaultTextInterpreter','latex'); 
rng(0); 		% for reproducing results
mtol = 1.0E-10;	% relative residual tolerance 

% ----------------------
% generate data matrices 
% ----------------------
r = 10; n = 100; 
A = randn(r,r)+1i*randn(r,r); 
B = randn(r,n)+1i*randn(r,n);
C = randn(n,r)+1i*randn(n,r); 
P0 = randn(n,n)+1i*randn(n,n); 
P1 = randn(n,n)+1i*randn(n,n); 
d = 1; 		% degree of P

% sample approximate eigenvalue z with controlled 'backward error'
epsln = 1.0E-1;
dA = A + (randn(r,r)+randn(r,r)*1i)*epsln; 
dB = B + (randn(r,n)+randn(r,n)*1i)*epsln;
dC = C + (randn(n,r)+randn(n,r)*1i)*epsln;
dP0 = P0 + (randn(n,n)+randn(n,n)*1i)*epsln; 
dP1 = P1 + (randn(n,n)+randn(n,n)*1i)*epsln; 
ee = eig([dA, dB; dC, dP0], blkdiag(eye(r),dP1));
z = ee(end); 	

% generate SRQ2 coefficient matrices
Az = A-eye(r)*z;
Pz = P0 - P1*z;
prob = "ABCP";
[G1, G2, G3, a1, a2, b1, b2] =buildSRQ2(Az, B, C, Pz, z, d, prob);

% ---------------------------------------------------
% compare SCF and RTR with different starting vectors
% ---------------------------------------------------
nn = 20; 	% number of starting vectors
V00 = randn(n+r,nn)+1i*randn(n+r,nn); % random initial vectors
TT = []; % timing history
RR = []; % gradient history
FF = []; % function value history
INN = [];

for i = 1:nn
	% - starting vector
	v00 = V00(:,i); 
	v00 = v00/norm(v00);

	% - run SCF
	tic;
	[v0, fv0, OBJFX, RRESD, RESD, neig] = runscf2(G1,G2,G3,a1,a2,b1,b2,v00,mtol);
	t1 = toc;

	% - run Riemannian trustregions
	tol = RESD(end); % tolerance conversion: manopt uses residual norm 
	tic;
	[v1, fv1, info] = runrtr(G1,G2,G3,a1,a2,b1,b2,v00,tol);
	t2 = toc;
	RESDRTR = arrayfun(@(info) info.gradnorm, info);
	OBJRTR = arrayfun(@(info) info.cost, info);
    
	% - print results
	disp([prob, num2str(i)])
	disp('Timing: SCF / RTR / Speedup');
	disp([t1, t2, t2/t1]);
    disp('Inner: SCF / RTR');
    inn = [info.numinner];
    disp([neig, inn(end)]);
	disp('Resd: SCF / RTR');
	disp([RESD(end), RESDRTR(end)]);
	disp('Objective Value: SCF / RTR');
	disp([OBJFX(end), OBJRTR(end)]);

	% - save computation history
	TT = [TT; t1, t2];
    INN = [INN; neig, inn(end)];
	RR = [RR; RESD(end), RESDRTR(end)];
	FF = [FF; OBJFX(end), OBJRTR(end)];

    % - draw convergence history 
	figure(1); 
	semilogy(RESD,'--o'); hold on; 
	semilogy(RESDRTR,'--+'); 
	figure(2);
	plot(OBJFX, '--o'); hold on;
	semilogy(OBJRTR,'--+'); 
end

% label figures
figure(1);
legend('SCF', 'RTR');
xlabel('iteration');
ylabel('residual norms');
figure(2);
legend('SCF', 'RTR');
xlabel('iteration');
ylabel('objective values');

% print timing statistics  
tt = mean(TT);
vt = max(abs(TT-tt));

% --------------------------------
% visualization: global optimality 
% --------------------------------
% sample supporting points of joint numerical range 
[P, V00] = samplejnr3(G1,G2,G3);

% find minimal function value over sampled points
objf = @(x,y,z) x./(a1+b1*z) + y./(a2+b2*z);
objfp= arrayfun(objf, P(:,1), P(:,2), P(:,3));
objfpmin = min(objfp);

% print results 
disp('Objective Value: SCF / RTR / JNR');
disp([OBJFX(end), OBJRTR(end), objfpmin]);
disp('Backerr: SCF / RTR / JNR');
disp(sqrt([OBJFX(end), OBJRTR(end), objfpmin]));
disp('Timing: SCF / RTR / Speedup');
disp([tt(1), tt(2), tt(2)/tt(1)]);

% draw joint numerical range and contour surface
h = showjnr3(P);
axis square;
hold on; 
a = real(v0'*G1*v0); b = real(v0'*G2*v0); c = real(v0'*G3*v0);
showcontour(a1,b1,a2,b2, [a,b,c], [-200,800,1.0,1.01], h);

%END
return