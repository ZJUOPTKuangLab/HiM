function plotPhasor(sB,gB,a,b)
%Spatial phasor plotting
%a,b: the slope and intercept of the projection line

yshift = 0.2;
sB = sB+yshift;
b = b+yshift;

mapResolution = 500;
PhasorMap = hist3([sB(:),gB(:)],'Edges', {1/mapResolution:1/mapResolution:1,1/mapResolution:1/mapResolution:1});
cmp = [zeros(3,33) jet(256)']';
cmp(1:33,1) = 1:-1/32:0;
cmp(1:33,2) = 1:-1/32:0;
cmp(1:33,3) = 1:-0.5/32:0.5;
figure;imagesc(PhasorMap(1:0.6*mapResolution,1:0.6*mapResolution)); colormap(cmp);
lx = (0.6-b)/a:0.001:-b/a;
hold on;plot(lx*mapResolution,(a*lx+b)*mapResolution,'--k','LineWidth',2);

hold on;quiver(-0.03*mapResolution,0.2*mapResolution,0.700*mapResolution,0,'color','k','LineWidth',2);
hold on;quiver(0,-0.03*mapResolution,0,0.7*mapResolution,'color','k','LineWidth',2);
text(-0.02*mapResolution,(yshift-0.04)*mapResolution,'0', 'FontSize', 16,'HorizontalAlignment','right');

xylabel = (-0.01:0.001:0.01)*mapResolution;
for ii= 0.2:0.2:0.4
    hold on;plot(ones(1,length(xylabel))*ii*mapResolution,xylabel+yshift*mapResolution,'-k','LineWidth',2);
    text(ii*mapResolution,(yshift-0.04)*mapResolution,num2str(ii), 'FontSize', 16,'HorizontalAlignment','center');
end
ii = 0.6;
text((ii-0.005)*mapResolution,(yshift-0.035)*mapResolution,'g', 'FontSize', 25,'HorizontalAlignment','right');

ii = 0;
hold on;plot(xylabel,ones(1,length(xylabel))*ii*mapResolution,'-k','LineWidth',2);
text((-0.02)*mapResolution,(ii+0.005)*mapResolution,num2str(ii-yshift), 'FontSize', 16,'HorizontalAlignment','right');
for ii = 0.4:0.2:0.4
    hold on;plot(xylabel,ones(1,length(xylabel))*ii*mapResolution,'-k','LineWidth',2);
    text((-0.02)*mapResolution,(ii+0.005)*mapResolution,num2str(ii-yshift), 'FontSize', 16,'HorizontalAlignment','right');
end
ii = 0.6;
text((-0.025)*mapResolution,(ii-0.015)*mapResolution,'s', 'FontSize', 25,'HorizontalAlignment','right');
axis xy; axis image;axis off;
end