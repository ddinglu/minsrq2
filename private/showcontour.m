function h = showcontour(a1,b1,a2,b2, y0, yz, h)
% function showcontour generate surface plot for contour g(x,y,z) = g(x0,y0,z0).
% g(x,y,z) = x/(a1+b1*z) + y/(a2+b1*z)

if nargin < 2 
	h = figure;
else
	figure(h);
end

%objf = @(x) real(x'*B1*x) + real(x'*B2*x)/real(x'*B3*x); % Hermitian matrices


g = @(y) y(1)./(a1+b1*y(3)) + y(2)./(a2+b2*y(3));

fc = g(y0);

ffx = @(y,z) (fc - y./(a2+b2*z) )*(a1+b1*z);
ffy = @(y,z) y;
ffz = @(y,z) z;


%fsurf(ffx,ffy,ffz, [ly-.2,uy,lz,uz],'FaceAlpha',0.9,'EdgeColor',[.7,.7,.7],'FaceColor',[.8,.8,.8]);
fsurf(ffx,ffy,ffz, yz,'FaceAlpha',1,'EdgeColor','k'); %,'FaceColor','y');

plot3(y0(1), y0(2), y0(3),'o', 'color', 'red', 'Markersize', 10, 'linewidth', 2, 'MarkerFacecolor','w');

%axis([lx,ux,ly,uy,lz,uz]);


