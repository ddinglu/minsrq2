function [P, V00] = samplejnr3(A,B,C,nn)
% function samplejnr3 samples nn (800) boundry points for the
% joint numerical range of the Hermitian matrix tuple (A,B,C).

if nargin < 4, nn = 800; end 
n = ceil(sqrt(nn/2));

% - uniformly sample nn spherical coordinates 
THETA = linspace(0,2*pi, 2*n);
PHI = linspace(0, pi, n);

% - compute boundary points alone sampled search direction 
%	v = [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)]
%H2 = @(d) d(1)*A + d(2)*B + d(3)*C;
m = size(A,1);
ea = sort(real(eig(A))); ca = (ea(end)+ea(1))/2; sA = (A-ca*eye(m))/(ea(end)-ea(1));
eb = sort(real(eig(B))); cb = (eb(end)+eb(1))/2; sB = (B-cb*eye(m))/(eb(end)-eb(1));
ec = sort(real(eig(C))); cc = (ec(end)+ec(1))/2; sC = (C-cc*eye(m))/(ec(end)-ec(1));
H2 = @(d) d(1)*sA + d(2)*sB + d(3)*sC;

P = [];
for i = 1:length(THETA)
	theta = THETA(i);
	for j = 1:length(PHI)
		phi = PHI(j);
		v0 = [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)]; % normal dir

    	HH = H2(v0);
    	[vv,ee] = eig(HH+HH');
		[~,idx] = min(real(diag(ee))); x0 = vv(:,idx);
		V00{i,j} = x0;

		a = real(x0'*A*x0); 
		b = real(x0'*B*x0); 
		c = real(x0'*C*x0); 
		P = [P; [a,b,c]];
	end
end

return;
