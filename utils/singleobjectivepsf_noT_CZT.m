function [Ex,Ey,Ez] = singleobjectivepsf_noT_CZT(Obj,Beam,Scope,Option)
%SINGLEOBJECTIVEPSF_NOT_CZT calculates the PSF of single objective by 2D
%chirp-z transform (CZT) in Cartesian coordinate system.
%
% INPUT********************************************************************
% Obj.NA: scalar value, numerical aperture of objective lens
% Obj.n: scalar value, refractive index in focal space
% Beam.wavelength: scalar value, Beam.wavelength of light
% Beam.amp: M*M matrix, amplitude distribution on pupil plane
% Beam.phs: M*M matrix, phase distribution on pupil plane
% Beam.abr: M*M matrix, aberration distribution on pupil plane
% Beam.plr: M*M*3 matrix, polarization distribution on pupil plane
% Scope.xs: 1*N array, representing x axis of PSF
% Scope.ys: 1*N array, representing y axis of PSF
% Scope.zs: 1*N array, representing z axis of PSF
% Beam.PupilRes: scalar value, resolution of pupil
% Option.UseGpu: 0 or 1, option using GPU acceleration
% Option.Precision: 0 or 1, precision of numbers
%
% OUTPUT*******************************************************************
% [Ex,Ey,Ez]: electric field of focused field
%
% REFERENCES***************************************************************
% [1] M. Leutenegger, R. Rao, R. Leitgeb, T. Lasser, Fast focus field
%     calculations. Opt. Express 14, 11277-11291 (2006).
%
% [2] Y. Hu et al., Efficient full-path optical calculation of scalar and
%     vector diffraction using the Bluestein method. Light: Science &
%     Applications 9, 119 (2020).
%
% TIPS*********************************************************************
% (1). Ultral fast for xy plane calculation, but a bit slow in z direction
%
% TODO*********************************************************************
% 1. Use CZT algorithm on z direction.
%
% *************************************************************************
% LIU Xin
% liuxin2018@zju.edu.cn
% Apr.24, 2021

%% data initialization
if strcmp(Option.Precision,'single')
    Obj.NA = single(Obj.NA);
    Obj.n = single(Obj.n);
    Beam.wavelength = single(Beam.wavelength);
    Beam.amp = single(Beam.amp);
    Beam.phs = single(Beam.phs);
    Beam.abr = single(Beam.abr);
    Beam.plr = single(Beam.plr);
    Scope.xs = single(Scope.xs);
    Scope.ys = single(Scope.ys);
    Scope.zs = single(Scope.zs);
end

%% gpu data preparation
if Option.UseGpu == 1
    Obj.NA = gpuArray(Obj.NA);
    Obj.n = gpuArray(Obj.n);
    Beam.wavelength = gpuArray(Beam.wavelength);
    Beam.amp = gpuArray(Beam.amp);
    Beam.phs = gpuArray(Beam.phs);
    Beam.abr = gpuArray(Beam.abr);
    Beam.plr = gpuArray(Beam.plr);
    Scope.xs = gpuArray(Scope.xs);
    Scope.ys = gpuArray(Scope.ys);
    Scope.zs = gpuArray(Scope.zs);
end

%% generate pupil grid (theta, phi, rho)
r0 = 1;  % radius of pupil
[xx,yy] = meshgrid(linspace(-r0,r0,Beam.PupilRes));  % cartesian coordinate of pupil plane
[phi,rho] = cart2pol(xx,yy);  % polar coordinate of pupil plane

thetaMax = asin(Obj.NA/Obj.n);  % maximum convergent angle of objective
sinTheta = sin(thetaMax).*rho;
sinTheta(rho>r0) = 0;
theta = asin(sinTheta);  % convergance angle of each ray (theta)

%% interpolation of pupil
if ~(size(Beam.amp) == Beam.PupilRes)
    Beam.amp = pupilInterp(Beam.amp,Beam.PupilRes);
end
if ~(size(Beam.phs) == Beam.PupilRes)
    Beam.phs = pupilInterp(Beam.phs,Beam.PupilRes);
end
if ~(size(Beam.abr) == Beam.PupilRes)
    Beam.abr = pupilInterp(Beam.abr,Beam.PupilRes);
end

px = Beam.plr(:,:,1);
py = Beam.plr(:,:,2);
% pz = Beam.plr(:,:,3);

if ~(size(px) == Beam.PupilRes)
    px = pupilInterp(px,Beam.PupilRes);
end

if ~(size(py) == Beam.PupilRes)
    py = pupilInterp(py,Beam.PupilRes);
end

% if ~(size(pz) == Beam.PupilRes)
%     pz = pupilInterp(pz,Beam.PupilRes);
% end

%% pupil function
E0 = Beam.amp.*exp(1i.*(Beam.phs+Beam.abr));

% remove parts outside numerical aperture
E0(rho>r0) = 0;

%% check sampling condition
kc = Obj.NA/Beam.wavelength;  % cut-off frequency of objective (in k-space)
dk = 2*kc/(Beam.PupilRes-1);  % sampling period in k-space
fs_xy = 1/dk;  % sampling frequency in k-space

fx_max = max(abs(Scope.xs));
fy_max = max(abs(Scope.ys));
fz_max = max(abs(Scope.zs));

dkz_max = Obj.n/Beam.wavelength*...
    (sqrt(1-(sin(thetaMax)*(Beam.PupilRes-3)/(Beam.PupilRes-1))^2)-cos(thetaMax));
fs_z_min = 1/dkz_max;

% Nyquist–Shannon sampling theorem
sampling_condition_x = fs_xy > 2*fx_max;
sampling_condition_y = fs_xy > 2*fy_max;
sampling_condition_z1 = fs_z_min > 2*fz_max;

% sampling condition for kz
% the phase kz*z must not change by more than pi between neighboring
% sampling points in the pupil plane.
M0 = 2*Obj.NA^2/sqrt(Obj.n^2-Obj.NA^2)*(max(abs(Scope.zs))./Beam.wavelength);
M0 = round(M0);
sampling_condition_z2 = Beam.PupilRes > 2*M0;

if ~(sampling_condition_x && sampling_condition_y &&...
        sampling_condition_z1 && sampling_condition_z2)
    pupilResX = round((Beam.PupilRes-1)*2*fx_max/fs_xy) + 1;
    pupilResY = round((Beam.PupilRes-1)*2*fy_max/fs_xy) + 1;
    
    Chi = Beam.wavelength/Obj.n*(1/(2*fz_max))+cos(thetaMax);
    if fz_max == 0
        pupilResZ1 = 1;
    else
        pupilResZ1 = 2/(1-sqrt(1-Chi^2)/sin(thetaMax)) + 1;
    end
    pupilResZ2 = 2*M0;
    
    % the required least pupil resolution
    pupilResReq = max([pupilResX,pupilResY,pupilResZ1,pupilResZ2]);
    warning(['Pupil resolution should be larger than ',...
        num2str(pupilResReq), ' avoiding aliasing!']);
end

%% polarization in image space (exit pupil)
[Px,Py,Pz] = coortrans(px,py,0,theta,phi,'o2i');

%% amplitude apodization factor for energy conservation
% different objective may have different AF!!!
AAF = sqrt(cos(theta));  % only for objective obeying the sine condition

%% plane waves in image space
% amplitude projection factor from (theta, phi) to (kx, ky) coordinate for
% integral (FFT)
APF = 1./cos(theta);

%% amplitude transformation
E0 = E0.*AAF.*APF;

%% the coefficient in front of Debye integral
f = r0/sin(thetaMax);  % focal length
prefix = -1i*Obj.n*f/Beam.wavelength;
E0 = prefix*E0;

%% plane waves in image space (three polarization components)
Ex0 = E0.*Px;
Ey0 = E0.*Py;
Ez0 = E0.*Pz;

%% CZT (chirp-z transform)
lx = length(Scope.xs); ly = length(Scope.ys); lz = length(Scope.zs);
if lx ~= 1
    pixelSizeX = Scope.xs(2)-Scope.xs(1);
else
    pixelSizeX = 0;
end

if ly ~= 1
    pixelSizeY = Scope.ys(2)-Scope.ys(1);
else
    pixelSizeY = 0;
end

% start point in observation (Fourier) plane
fx1 = Scope.xs(1);
fy1 = Scope.ys(1);

% start point on unit circle (Equation S15)
Ax = exp( 1i*2*pi*fx1/fs_xy);
Ay = exp( 1i*2*pi*fy1/fs_xy);

% step size on unit circle (Equation S16)
Wx = exp(-1i*2*pi*pixelSizeX/fs_xy);
Wy = exp(-1i*2*pi*pixelSizeY/fs_xy);

% czt shift (Equation S18)
Mshiftx = -(Beam.PupilRes-1)/2;
Mshifty = -(Beam.PupilRes-1)/2;

% phase shift (Equation S19)
[xss,yss] = meshgrid(Scope.xs.*Mshiftx/fs_xy,Scope.ys.*Mshifty/fs_xy);
Pshift = exp(-1i*2*pi*(xss+yss));

% wave number in vacuum
k0 = 2*pi/Beam.wavelength;
kz = k0*Obj.n*cos(theta);

%% batch segmentation
% avoid Option.BatchSize = 0
if Option.BatchSize > lz
    Option.BatchSize = lz;
end

% number of blocks
batchNumber = floor(lz/Option.BatchSize);

% length of the last block
lastBatchSize = mod(lz,Option.BatchSize);

%% creat block matrix
batchVector = ones(1,batchNumber).*Option.BatchSize;
if lastBatchSize ~= 0
    batchVector = [batchVector,lastBatchSize];
end

totalBatches = length(batchVector);  % total cycles

zc = mat2cell(reshape(Scope.zs,1,1,lz),1,1,batchVector);

% memory preallocation
Ex = cell(1,totalBatches);
Ey = Ex;
Ez = Ex;

if Option.TimeIt == 1
    tstart = tic;
end

for ii = 1:totalBatches
    defocusTerm = exp(1i*kz.*zc{ii});
    
    Ewx = Ex0.*defocusTerm;
    Ewy = Ey0.*defocusTerm;
    Ewz = Ez0.*defocusTerm;
    
    EHold = myczt(Ewx,ly,Wy,Ay);
    EHold = permute(EHold,[2 1 3]);
    EHold = myczt(EHold,lx,Wx,Ax);
    Ex{ii} = permute(EHold,[2 1 3]).*Pshift;
    
    EHold = myczt(Ewy,ly,Wy,Ay);
    EHold = permute(EHold,[2 1 3]);
    EHold = myczt(EHold,lx,Wx,Ax);
    Ey{ii} = permute(EHold,[2 1 3]).*Pshift;
    
    EHold = myczt(Ewz,ly,Wy,Ay);
    EHold = permute(EHold,[2 1 3]);
    EHold = myczt(EHold,lx,Wx,Ax);
    Ez{ii} = permute(EHold,[2 1 3]).*Pshift;
    
    textwaitbar(ii, totalBatches, 'PSF in progress');
end
if Option.TimeIt == 1
    toc(tstart);
end

Ex = cat(3,Ex{:});
Ey = cat(3,Ey{:});
Ez = cat(3,Ez{:});

if Option.UseGpu == 1
    Ex = gather(Ex);
    Ey = gather(Ey);
    Ez = gather(Ez);
end

if strcmp(Option.Precision,'single')
    Ex = double(Ex);
    Ey = double(Ey);
    Ez = double(Ez);
end

%% energy efficiency calculation
%{
CC = AAF.*APF.*prefix;
pupilPixSize = 1/dk;
pupilMask = rho<=r0;
effAmp = Beam.amp(pupilMask).*CC(pupilMask);

EnergyIn = sum(abs(effAmp).^2,'all').*pupilPixSize^2;
EnergyOut = sum(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2,'all').*pixelSizeX.*pixelSizeY;
EnergyEfficiency = EnergyOut./EnergyIn;
%}
end