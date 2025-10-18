function [imgPRIF,mask] = spatialPhasor_mask(imgPR,rdets,binning)
%Main function of spatial phasor analysis
%imgPR:  ISM data after pixel reassignment, 
%        a four-dimensional array of [xSize,ySize,detectorNumber,zSize]
%rdets:  normalized detector position
%binning:Sliding window radius in pixel size.

imgRes = size(imgPR,1);
detNum = size(imgPR,3);
n = size(imgPR,4);
imgRdets = repmat(reshape(rdets,1,1,detNum),imgRes,imgRes,1);

imgPR1 = zeros(imgRes,imgRes,detNum);
imgDark = zeros(imgRes,imgRes,n,'logical');
imgDark1 = zeros(imgRes,imgRes,n,'logical');
gB = zeros(imgRes,imgRes,n);
sB = zeros(imgRes,imgRes,n);
binning = binning+0.5;
for jj = 1:n
    kernel = fspecial('disk',binning);
    for ii = 1:19
        imgPR1(:,:,ii) = imfilter(imgPR(:,:,ii,jj),kernel);
    end
%     imgDark(:,:,jj) = max(imgPR1,[],3)>0.5;
    imgDark(:,:,jj) = max(imgPR1,[],3)>1;
    gB(:,:,jj) = sum(imgPR1.*cos(2*pi*imgRdets),3)./sum(imgPR1,3);
    sB(:,:,jj) = sum(imgPR1.*sin(2*pi*imgRdets),3)./sum(imgPR1,3);

    imgDark1(:,:,jj) = max(imgPR1,[],3)>1;
end

Coeffs = fit(sB(imgDark1),gB(imgDark1),'poly1');
a = 1/Coeffs.p1;
b = -Coeffs.p2/Coeffs.p1;
plotPhasor(sB,gB,a,b);

empmatr = ones(imgRes,imgRes,n);
x1 = 0;
y1 = a*x1+b;
xmatr = (gB+a*sB-a*b*empmatr)/(1+a^2);
ymatr = (a*gB+a^2*sB+b*empmatr)/(1+a^2);
dpie = sqrt((x1*empmatr-xmatr).^2+(y1*empmatr-ymatr).^2);
dpie = dpie-min(dpie(imgDark));
dpie(~imgDark) = max(dpie(imgDark));


[lowOut,highOut] = phasor2hist(dpie(imgDark));
% lowOut = 0.25;
% highOut = 0.4;

mask = zeros(imgRes,imgRes,detNum,n);
hOver = max(highOut-1,0);
for  jj = 1:n
mask1 = 1-imadjust(dpie(:,:,jj)-hOver,[lowOut,highOut]-hOver,[0,1]);
mask(:,:,:,jj) = repmat(filterGauss(mask1,imgRes/3/binning),1,1,detNum);
end
imgPRIF = imgPR.*mask;
end