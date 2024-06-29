function h = showjnr3(P, h)
% function showjnr3 generate surface plot for the 3D joint numerical range using 
% sample points in P.

if nargin < 2 
	h = figure;
else
	figure(h);
end

% surface plot
KK = boundary(P(:,1), P(:,2), P(:,3));
trisurf(KK,P(:,1), P(:,2), P(:,3), 'Facecolor',[.9,.9,.9],'FaceAlpha',0.9, 'EdgeColor', 'k','linewidth', 1); hold on;
plot3(P(:,1),P(:,2),P(:,3), '.', 'Markersize', 6, 'color', [.1,.1,.1], 'linewidth', 1); hold on;
xlabel('$y_1$', 'Fontsize', 16);
ylabel('$y_2$', 'Fontsize', 16); 
zlabel('$y_3$', 'rotation', 0, 'Fontsize', 16);
grid on; axis equal;


% estimate bounding box for the JNR 
lbd = min(P); ubd = max(P);
axis( reshape([lbd;ubd], 6,1) );

return

