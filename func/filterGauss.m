% function Iout = filterGauss(Iin,sigma)
% %Frequency domain Gaussian filter
% 
% imgRes = size(Iin,1);
% [xx,yy]=meshgrid(-imgRes/2:imgRes/2-1,-imgRes/2:imgRes/2-1);
% rr = sqrt(xx.^2+yy.^2);
% fGauss = exp(-rr.^2/2/sigma^2);
% Iout = ifft2(ifftshift(fftshift(fft2(Iin)).*fGauss));
% Iout = abs(Iout);
% end


function Iout = filterGauss(Iin, sigma)
% 修正版频域高斯滤波
imgRes = size(Iin, 1);

% 生成正确的频率网格（确保零频率在中心）
[xx, yy] = meshgrid(-floor(imgRes/2):floor((imgRes-1)/2), ...
                    -floor(imgRes/2):floor((imgRes-1)/2));
rr = sqrt(xx.^2 + yy.^2);

% 生成频域高斯滤波器
fGauss = exp(-rr.^2 / (2 * sigma^2));

% 傅里叶变换并中心化
I_fft = fft2(Iin);
I_fft_shifted = fftshift(I_fft);  % 将零频率移到中心

% 应用滤波器
I_filtered_shifted = I_fft_shifted .* fGauss;

% 反中心化并逆变换
I_filtered = ifftshift(I_filtered_shifted);  % 移回原位
Iout = ifft2(I_filtered);

% 取绝对值确保输出为实数
Iout = abs(Iout);
end