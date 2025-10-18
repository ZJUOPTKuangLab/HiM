 function Iobj = HiMrecovery(raw_data,OTFex,OTFISM,patterneff,widefield,max_iter)

fidelity = 20;
constraint = 1;
iter = 40;

raw_data = gpuArray(single(raw_data));
OTFex = gpuArray(single(OTFex));
OTFISM = gpuArray(single(OTFISM));
patterneff = gpuArray(single(patterneff));
widefield = gpuArray(single(widefield));

num = size(raw_data,3);
imgRes = size(raw_data,1);
widefield = widefield/max(widefield(:));
Iobj = widefield;

n = 1;

while n <= max_iter
    tic;
    for i = 1 : num
        Iop = Iobj .* patterneff(:,:,i);
        fop = fftshift(fft2(Iop));
        fop1 = fop+conj(OTFex).* (fftshift(fft2(raw_data(:,:,i)))-OTFex.*fop);
        Iop1 = ifft2(ifftshift(fop1));
        Iobj1 = Iobj+patterneff(:,:,i).*(Iop1-Iop)./(max(max(patterneff(:,:,i)))^2);

        Iobj1(Iobj1<0) = 0;
        Iobj = abs(Iobj1);
        
    end
    
    figure(38);imshow(abs(Iobj),[]);
   
    IobjVik = Iobj;
    IobjEx = Iobj/max(Iobj(:));
    if n == 1
        IobjFliter = ifft2(ifftshift(fftshift(fft2(IobjEx)).*OTFISM));
        [~,ratioI] = matchInt(IobjEx,IobjFliter);
        imgISM = ratioI*widefield;
    else
        IobjEx = matchInt(IobjEx,IobjSpa);
    end
    
    IobjSpa = sparseDeconvCons(IobjEx,imgISM,OTFISM,fidelity,constraint*max((10-n+1),1),iter,1,1);
    IobjSpa(IobjSpa<0) = 0;
    figure(40);imshow(abs(IobjSpa),[]);
    
    if n < max_iter
    IobjSpa = filterGauss(IobjSpa,imgRes/2);
    Iobj = matchInt(IobjSpa,IobjVik);
    else
    Iobj = min(matchInt(IobjSpa,IobjVik),IobjVik);
    end
    
    ttime = toc ;
    disp(['  iter ' num2str(n) ' | ' num2str(max_iter) ', took ' num2str(ttime) ' secs']);
    n=n+1;
    
end

Iobj = gather(double(Iobj));

 end

 
 
 