function PSF_XY = calcPSF(lamda,psize,na,flagGpu)
Option.UseGpu = flagGpu;  % gpu flag
Option.Precision = 'single';
Option.FocusingMethod = 'czt';
Option.TimeIt = 0;

Option.BatchSize = 100;

xHalfScope = 1200;
yHalfScope = xHalfScope;

obj.NA = na;
obj.n = 1.515;

pupilRes = 100;
psfRes = xHalfScope*2/psize;

beam.PupilRes = pupilRes;
beam.wavelength = lamda;
load('Gauss.mat');
beam.amp = amp;
beam.phs = zeros(pupilRes);
beam.abr = zeros(pupilRes);
beam.plr = CircularPolarization(pupilRes,'l');

%% XY
scopeXY.xs = linspace(-xHalfScope+psize/2,xHalfScope-psize/2,psfRes);
scopeXY.ys = linspace(-yHalfScope+psize/2,yHalfScope-psize/2,psfRes);
scopeXY.zs = 0;

[Ex_XY, Ey_XY, Ez_XY] = singleobjectivepsf(obj,beam,scopeXY,Option);
PSF_XY = abs(Ex_XY).^2 + abs(Ey_XY).^2 + abs(Ez_XY).^2;

PSF_XY = PSF_XY-PSF_XY(length(PSF_XY)/2,1);
PSF_XY(PSF_XY<0) = 0;
PSF_XY = PSF_XY/max(PSF_XY(:));

%% figure
f = figure('Color','w');
imagesc(scopeXY.xs,scopeXY.ys,PSF_XY);
axis image xy;
xlabel('x');
ylabel('y');
colormap(gca,'hot');
colorbar;
title('XY');
end