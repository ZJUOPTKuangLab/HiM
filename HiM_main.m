clear;clc;close all;
addpath('.\func');
addpath('.\utils');

%%
filePath = 'H:\HiM\data\';

flagGpu = 1;
k0=0.8;
na=1.45;
lamdaEX=640;
lamdaEM=690;
psize=30;
pinsize = 80*1000/(100*300/70*60/50)/2;
freqRound = 4;

%% load data and reassignment
imgAiry = readAiryData(filePath,1);
[yoffset,xoffset] = calcShift(imgAiry);
imgPR = pixelReassign(imgAiry,round(yoffset),round(xoffset));

%% deBKGD
rdets = zeros(1,19);
rdets(2:7) = 0.25;
rdets(8:2:18) = 0.75;
rdets(9:2:19) = 1;
imgPR = spatialPhasor(imgPR,rdets,4);
imgAiry = pixelReassign(imgPR,-yoffset,-xoffset);
Confocal = sum(imgAiry,3);
imgISM = sum(imgPR,3);
imgRes = size(imgAiry,1);
detNum = size(imgAiry,3);

%% PSF
[xx,yy]=meshgrid(-(imgRes-1)/2:(imgRes-1)/2,-(imgRes-1)/2:(imgRes-1)/2);
freqUnit = sqrt(xx.^2+yy.^2)/imgRes/psize;
PSFex = calcPSF(lamdaEX/1.3,psize,na,flagGpu);
padSize = (imgRes-length(PSFex))/2;
PSFex = padarray(PSFex,[padSize,padSize],0,'both');
PSFex = PSFex/sum(PSFex(:));
OTFex = fftshift(fft2(PSFex));
PSFem = calcPSF(lamdaEM,psize,na,flagGpu);
PSFem = padarray(PSFem,[padSize,padSize],0,'both');
Fpinhole = besselj(1,freqUnit*2*pi*pinsize)./freqUnit;
OTFem = fftshift(fft2(PSFem)).*Fpinhole;
OTFem = OTFem/max(OTFem(:));
PSFISM = PSFex.*abs(ifft2(ifftshift(OTFem)));
PSFISM = PSFISM/sum(PSFISM(:));
OTFISM = fftshift(fft2(PSFISM));

%% virtual modulation
kc = na*4*pi/lamdaEM;
pkc = na*2*psize*imgRes/lamdaEX;
ydet = zeros(imgRes,imgRes,detNum);
xdet = zeros(imgRes,imgRes,detNum);
for ii = 2 :detNum
    ydet(:,:,ii) = 2*yoffset(ii-1);
    xdet(:,:,ii) = 2*xoffset(ii-1);
end
k = k0*kc*psize;
[vsd,patterneff] = HiMmodulation(imgAiry,abs(OTFem),k,ydet,xdet,freqRound);

%% FP demodulation and sparse joint deconvolution
Iobj = HiMrecovery(vsd,abs(OTFex),abs(OTFISM),patterneff,imgISM,40);

%% save result
imwrite(uint16(Confocal),'.\result\Confocal.tif');
imwrite(uint16(imgISM),'.\result\ISM.tif');
imwrite(uint16(Iobj),'.\result\HiM.tif');



