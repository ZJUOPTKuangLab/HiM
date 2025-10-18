function [lowOut,highOut] = phasor2hist(imgPhasor)
%Calculation of in-focus and out-of-focus signal centers 
%from histogram of phasor points

phasorLine = imgPhasor(:);
phasorLine(phasorLine == 0) = [];
phasorLine(phasorLine == 1) = [];
Nbin = 200;
[phasorH,bins] = hist(phasorLine,Nbin);
options = fitoptions('gauss2', 'Lower', zeros(1,6));
phasorFit = fit(bins',phasorH','gauss2',options);

figure;plot(phasorFit);
hold on;histogram(phasorLine,Nbin);

if phasorFit.b1 < phasorFit.b2
    ifp = phasorFit.b1;
    ifs = phasorFit.c1/sqrt(2);
    ofp = phasorFit.b2;
    ofs = phasorFit.c2/sqrt(2);
else
    ifp = phasorFit.b2;
    ifs = phasorFit.c2/sqrt(2);
    ofp = phasorFit.b1;
    ofs = phasorFit.c1/sqrt(2);
end

if ofp-2*ofs > ifp+2*ifs
    lowOut = ifp+2*ifs;
else
    lowOut = ifp;
end


if ifp+2*ifs < ofp-2*ofs
    highOut = ofp-2*ofs;
elseif ifp+2*ifs < ofp
    highOut = ofp+2*ofs;
else
    highOut = ofp+2*ofs;
end

end