function imgPR = pixelReassign(imgAiry,yshift,xshift)
%Pixel reassignment using shift vectors

imageRes = size(imgAiry,1);
det_num = size(imgAiry,3);
n = size(imgAiry,4);

ystart = -yshift.*(yshift < 0) + 1;
yend = imageRes - yshift.*(yshift >= 0);
xstart = -xshift.*(xshift < 0) + 1;
xend = imageRes - xshift.*(xshift >= 0);
imgPR = zeros(imageRes,imageRes,det_num,n);
imgPR(:,:,1,:) = imgAiry(:,:,1,:);
for k=1:det_num-1
    yIndex = ystart(k):yend(k);
    xIndex = xstart(k):xend(k);
    imgPR(yIndex,xIndex,k+1,:) = imgAiry(yIndex+yshift(k),xIndex+xshift(k),k+1,:);
end
end