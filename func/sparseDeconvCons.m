function Iout=sparseDeconvCons(Iini,Iblur,OTF,fidelity,constraint,iter,gpu,mu)

OTF = single(ifftshift(OTF));

if nargin < 6 || isempty(iter)
    iter = 1;
end
if nargin < 6 || isempty(gpu)
    gpu = cudaAvailable;
end
if nargin < 7 || isempty(mu)
    mu = 1;
end
Iini = single(Iini);
Iblur = single(Iblur);
b = zeros(size(Iini),'single');

if gpu == 1
    Iini = gpuArray(Iini);
    Iblur = gpuArray(Iblur);
    OTF = gpuArray(OTF);
    b = gpuArray(b);
end

norm = single((fidelity/mu)*OTF.^2+constraint/mu+1);
Ibb0 = (fidelity/mu)*fftn(Iblur).*OTF;
Iout = Iini;
for iter = 1:iter
    d = sign(Iout+b).*max(abs(Iout+b)-1/mu,0);
    b = b+(Iout-d);
    Ibb = Ibb0+(constraint/mu)*fftn(matchInt(Iini,Iout));
    Iout = real(ifftn((fftn(d-b)+Ibb)./norm));
end
Iout(Iout<0) = 0;
end
