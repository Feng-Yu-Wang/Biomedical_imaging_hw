% EE441000 ¤ý»ñ¦t 105081041 Hw2 12/7/2019
clc; close all;

% -----------Read data----------- %
% Real component
fid = fopen('hw2_r.dat','rb');
data_r = fread(fid,[256 256],'int32').';
fclose(fid);

% Imaginary component
fid = fopen('hw2_i.dat','rb');
data_i = fread(fid,[256 256],'int32').';
fclose(fid);

% k-space data
k_space_data = data_r+sqrt(-1)*data_i; 

% -----------(a)----------- %
% Real part
r_max = max(data_r);
r_min = min(data_r);
r_range = r_max-r_min;
r_max_row = max(data_r,[],2);
r_min_row = min(data_r,[],2);
r_range_row = r_max_row-r_min_row;
figure
subplot(1,2,1)
plot(r_range)
title('Data Range of Column')
xlabel('Pixel')
subplot(1,2,2)
plot(r_range_row)
title('Data Range of Row')
xlabel('Pixel')
r_ifft = 20*log10(abs(ifftshift(ifft2(data_r))));

figure
imagesc(r_ifft);
axis image
axis off
title('Image from Real Part')
colorbar
colormap(gray)

% Imaginary part
i_max = max(data_i);
i_min = min(data_i);
i_range = i_max-i_min;
i_max_row = max(data_i,[],2);
i_min_row = min(data_i,[],2);
i_range_row = i_max_row-i_min_row;
figure
subplot(1,2,1)
plot(i_range)
title('Data Range of Column')
xlabel('Pixel')
subplot(1,2,2)
plot(i_range_row)
title('Data Range of Row')
xlabel('Pixel')

i_ifft = 20*log10(abs(ifftshift(ifft2(data_i))));

figure
imagesc(i_ifft);
axis image
axis off
title('Image of Imaginary Part')
colorbar
colormap(gray)

% -----------(b)----------- %
k_space_image = 20*log10(abs(ifftshift(ifft2(k_space_data))));
figure
imagesc(k_space_image)
axis image
axis off
colorbar
colormap(gray)

% -----------(c)----------- %
r_mean = mean2(data_r);
i_mean = mean2(data_i);
data_r_rm = data_r-r_mean;
data_i_rm = data_i-i_mean;
data_r_rm(129,129) = 10*r_mean;
data_i_rm(129,129) = 10*i_mean;
k_space_data_rm = data_r_rm+sqrt(-1)*data_i_rm;
k_space_image_rm = 20*log10(abs(ifftshift(ifft2(k_space_data_rm))));
figure
imagesc(k_space_image_rm)
axis image
axis off
colorbar
colormap(gray)

% -----------(d)----------- %
k_space_data_new = zeros(256,256);
k_space_data_new(:,2:2:end) = k_space_data_rm(:,2:2:end);
k_space_image_new = 20*log10(abs(ifftshift(ifft2(k_space_data_new))));
figure
imagesc(k_space_image_new)
axis image
axis off
colormap(gray)
colorbar

% -----------(e)----------- %
HF(129:256,:) = k_space_data_rm(129:256,:);
HF(1:127,:) = conj(flip(flip(k_space_data_rm(130:256,:),1),2));
HF_data = 20*log10(abs(ifftshift(ifft2(HF))));
figure
imagesc(HF_data)
axis image
axis off
colormap(gray)
colorbar

% -----------(f)----------- %
% EMI at (160,160)
EMI_1 = k_space_data_rm;
EMI_1(160,160) = 10^6+sqrt(-1)*10^6;
EMI_1_image = 20*log10(abs(ifftshift(ifft2(EMI_1))));
figure
imagesc(EMI_1_image)
title('EMI at (160,160)')
axis image
axis off
colormap(gray)
colorbar

% EMI at (129,160)
EMI_2 = k_space_data_rm;
EMI_2(129,160) = 10^6;
EMI_2_image = 20*log10(abs(ifftshift(ifft2(EMI_2))));
figure
imagesc(EMI_2_image)
title('EMI at (129,160)')
axis image
axis off
colormap(gray)
colorbar

% EMI at (160,40)
EMI_3 = k_space_data_rm;
EMI_3(160,40) = 10^6;
EMI_3_image = 20*log10(abs(ifftshift(ifft2(EMI_3))));
figure
imagesc(EMI_3_image)
title('EMI at (160,40)')
axis image
axis off
colormap(gray)
colorbar