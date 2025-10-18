function[P] = CatmullRom_zoom(img,M)
[H,W,N] = size(img);
hc = linspace(0,H-1/M,H*M);
uh = [min(max(floor(hc)-1,0),H-1);min(max(floor(hc),0),H-1);min(max(floor(hc)+1,0),H-1);min(max(floor(hc)+1,0),H-1)];
uh = uh + 1;
Qh = (0.5.*[0,-1,2,-1;2,0,-5,3;0,1,4,-3;0,0,-1,1])*[ones(1,H*M);hc-floor(hc);(hc-floor(hc)).^2;(hc-floor(hc)).^3];
% Qh = (0.5.*[0,2,0,0;-1,0,1,0;2,-5,4,-1;-1,3,-3,1])*[ones(1,H*M);hc-floor(hc);(hc-floor(hc)).^2;(hc-floor(hc)).^3];
Ph = zeros(W,M*H,N);
for i = 1:N
    for j = 1:W
        for m = 1:M*H
            Ph_index(m,j) = img(uh(1,m),j,i).*Qh(1,m) + img(uh(2,m),j,i).*Qh(2,m) + img(uh(3,m),j,i).*Qh(3,m) + img(uh(4,m),j,i).*Qh(4,m);
            
        end
    end
    Ph(:,:,i) = Ph_index';
end


wc = linspace(0,W-1/M,W*M);
uw = [min(max(floor(wc)-1,0),W-1);min(max(floor(wc),0),W-1);min(max(floor(wc)+1,0),W-1);min(max(floor(wc)+2,0),W-1)];
uw = uw + 1;
Qw = (0.5.*[0,-1,2,-1;2,0,-5,3;0,1,4,-3;0,0,-1,1])*[ones(1,W*M);wc-floor(wc);(wc-floor(wc)).^2;(wc-floor(wc)).^3];
% Qw = (0.5.*[0,2,0,0;-1,0,1,0;2,-5,4,-1;-1,3,-3,1])*[ones(1,W*M);wc-floor(wc);(wc-floor(wc)).^2;(wc-floor(wc)).^3];
P = zeros(W*M,M*H,N);
for i = 1:N
    for j = 1:W*M
        for m = 1:M*H
            Pw_index(m,j) = Ph(uw(1,m),j,i).*Qw(1,m) + Ph(uw(2,m),j,i).*Qw(2,m) + Ph(uw(3,m),j,i).*Qw(3,m) + Ph(uw(4,m),j,i).*Qw(4,m);
            
        end
    end
    P(:,:,i) = Pw_index';
end
end