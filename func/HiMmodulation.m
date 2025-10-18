function [vsd,patterneff] = HiMmodulation(imgAiry,OTF,k,ydet,xdet,freqRound)
imgRes = size(imgAiry,1);
detNum = size(imgAiry,3);
yscan = repmat((-(imgRes-1)/2:(imgRes-1)/2)',1,imgRes,detNum);
xscan = repmat((-(imgRes-1)/2:(imgRes-1)/2),imgRes,1,detNum);
patternNum = freqRound*(1+freqRound)/2*3*3+1;

num=1;
patterneff = ones(imgRes,imgRes,patternNum);
vsd = zeros(imgRes,imgRes,patternNum);
vsd(:,:,1) = sum(imgAiry,3);
for ii=1:freqRound
    freqAng = 60/ii/180*pi;
    kp = ii/freqRound*k;
    for theta=0:freqAng:pi-freqAng
        kx = kp*cos(theta);
        ky = kp*sin(theta);
        for phase = 0:2*pi/3:4*pi/3
            num = num+1;
            pattern = (cos(kx.*(-1*xdet+xscan)+ky.*(-1*ydet+yscan)+phase)); 
            vsd(:,:,num) = sum(imgAiry.*pattern,3);
            fpotf = fftshift(fft2(pattern(:,:,1))).*OTF;
            patterneff(:,:,num) = ifft2(ifftshift(fpotf));
        end
    end
end
end