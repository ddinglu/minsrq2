% Example 1: Maximization of SRQ2 v.s. Maximization over JNR, for a randomly
% generated SRQ2 maximization of size n=3:
%   max (x'*B1*x)/(x'*x) + (x'*B2*x)/(x'*B3*x)
close all; clear;
warning('off', 'map:removing:combntns');
set(0,'defaultTextInterpreter','latex'); 

%-----------------------
% test data and matrices 
%-----------------------
rng(0);
n = 3;
B1 = [  0.64   -0.15   -0.38;
       -0.15    0.60   -0.22;
       -0.38   -0.22    0.56];
B2 = [  0.73    0.24   -0.07;
        0.24    0.52   -0.04;
       -0.07   -0.04    0.38];
B3 = diag([0.53, 0.97, 0.38]);

% ---------------------------------------------------------------
% solve SRQ2 maximization by Manopt with sampled starting vectors
% ---------------------------------------------------------------
prob.M = spherefactory(n);
prob.cost = @(x) x'*B1*x + (x'*B2*x)/(x'*B3*x); 
prob.egrad= @(x) 2*B1*x + 1/(x'*B3*x)*(2*B2*x) - (x'*B2*x)/(x'*B3*x)^2*(2*B3*x); 
opts.verbosity = 0;

% sample random starting vectors and run RTR from each vector 
V00 = randn(n,100);
rr = []; vv1 = [];
for ii = 1:size(V00,2)
    % run RTR for a solution
    v00 = V00(:,ii);
    v00 = v00/norm(v00);
    [v1, xcost1, info1, opts] = trustregions(prob, v00, opts);

    % check if NEPv holds at the solution
    Hx = @(x) B1 + 1/(x'*B3*x)*B2 - (x'*B2*x)/(x'*B3*x)^2*B3;
    H0 = Hx(v1/norm(v1));
    rr = [rr, real(v1'*H0*v1)];
    vv1 = [vv1, v1];

end
% find non-indentical solution
[~, idxx] = unique(round(rr,4));

% ------------------------------------------
% draw contour surface y(1) + y(2)/y(3) = fc
% ------------------------------------------
% sample boundary points of JNR 
[P, V00] = samplejnr3(B1,B2,B3);
lb = min(P); ub = max(P);
lx = lb(1); ly = lb(2); lz = lb(3);
ux = ub(1); uy = ub(2); uz = ub(3);

% draw JNR with level surface for each solution found 
objf = @(x) x'*B1*x + (x'*B2*x)/(x'*B3*x); 

for ii = 1:length(idxx)
    % draw JNR
    hh = figure(ii);
    showjnr3(P, hh);

    % draw level surface at solution
    v0 = vv1(:, idxx(ii));
    fc = objf(v0);
    ffx = @(y,z) fc - y/z;
    ffy = @(y,z) y;
    ffz = @(y,z) z;
    fsurf(ffx, ffy, ffz, [ly-.2,uy,lz,uz], 'FaceAlpha', 1, 'EdgeColor', 'k');
    
    % draw solution point
    a = v0'*B1*v0; b = v0'*B2*v0; c = v0'*B3*v0;
    plot3(a, b, c,'o', 'color', 'red', 'Markersize', 10, 'linewidth', 2, 'MarkerFacecolor','w');
    axis([lx,ux,ly,uy,lz,uz]);

    % print residual of NEPv at solution
    Hxx = Hx(v0); Hxx = (Hxx + Hxx')/2;
    lmin = min(real(eig(Hxx)));
    resdvi = norm(Hxx*v0-lmin*v0);

    disp('y/x/nepv residual:')
    [a,b,c]
    v0
    resdvi

    % check definiteness of (real) Hessian over the complement of the solution:
    % Consider real representations of complex vector/matrix 
    BB1 = blkdiag(B1,B1); BB2 = blkdiag(B2,B2); BB3 = blkdiag(B3,B3); xx = [v0;0*v0]; 
    rhess = 2*BB1 + 2*BB2/real(xx'*BB3*xx) - 2*(BB3)*(xx'*BB2*xx)/(xx'*BB3*xx).^2 - 4/(xx'*BB3*xx).^2*((BB2*xx)*(xx'*BB3)) - 4/(xx'*BB3*xx).^2*((BB3*xx)*(xx'*BB2)) + 8*(xx'*BB2*xx)/(xx'*BB3*xx).^3*((BB3*xx)*(xx'*BB3));
    rhess_proj = (eye(2*n)-xx*xx')*rhess*(eye(2*n)-xx*xx');
    disp('eig of projected hessian:')
    sort(real(eig(rhess_proj)))

end

% END
return