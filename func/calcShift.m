function [yshift,xshift] = calcShift(imgAiry)
%Calculating shift vectors for ISM data 
%using correlation function

imgAiry1 = permute(imgAiry,[1,2,4,3]);
detNum = size(imgAiry1,4);
xshift = zeros(1,detNum-1);
yshift = zeros(1,detNum-1);

FI1 = fftn(imgAiry1(:,:,:,1));
FRI1I1 = FI1.*conj(FI1);
RI1 = ifftn(FRI1I1);
RI1 = fftshift(RI1);
numI1 = find(RI1==max(RI1(:)));
[iI1,jI1,~] = ind2sub(size(RI1), numI1);

for k=2:detNum
    FIn = fftn(imgAiry1(:,:,:,k));
    FRI1In = FI1.*conj(FIn);
    R = ifftn(FRI1In);
    R = fftshift(R);
    num = find(R==max(R(:)));
    [i,j,~] = ind2sub(size(R), num);
    if length(i) == 1
        yshift(k-1) = iI1-i;
        xshift(k-1) = jI1-j;
    end
end
end