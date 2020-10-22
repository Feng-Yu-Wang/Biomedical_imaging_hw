% EE 441000 105081041 HW1 Part_I 10/31/2019
close all;

%--------------------------Parameters-----------------------------%
offset_time = 6.48; % offset time (us)
fs = 50;    % sampling rate (MHz)
aper_size = 6;  % aperture size (mm)
focal_pt = 5;   % focal point (mm)
dx = 0.05;    % distance between two successive scanning positions (mm)
v_sound = 1.54;  % speed of sound (mm/us)
%-----------------------------------------------------------------%

%--------------(a)--------------%
DR = 40;    % dynamic range (dB)
load points_rf_data;
data = points_rf_data; 
clear points_rf_data;
[m,n] = size(data);   % m: number of samples per scanline; n: number of scanlines
time = offset_time+(0:m-1)*(1/fs);   % time axis (us)
z_axis = time*v_sound/2;   % z axis (mm)
x_axis = -3+(0:n-1)*dx;    % x axis (mm)
env_data = abs(hilbert(data));   % envelope formation
norm_data = env_data/max(max(data)); % normalization with respect to maximum
log_data = 20*log10(norm_data+eps) + DR; % logarithm conversion
log_data(log_data < 0) = 0;

figure
image(x_axis,z_axis,log_data);
colormap(gray(DR))
colorbar;
title('Point targets')
xlabel('Lateral position (mm)')
ylabel('Depth (mm)')
axis image

%--------------(c)--------------%
scanline = 60;
figure
plot(z_axis,log_data(:,scanline));
xlabel('Depth (mm)')
ylabel('Signal amplitude (a.u.)');
xlim([5,20])
Title = sprintf('A-line Image of %dth Scanline ',scanline);
title(Title);

% axial resolution
scanned = log_data(:,scanline);
psf_1 = scanned(20:150); % PSF of point 1
psf_2 = scanned(190:330); % PSF of point 2
psf_3 = scanned(390:530); % PSF of point 3
psf_4 = scanned(590:730); % PSF of point 4
psf_5 = scanned(790:930); % PSF of point 5 

axial_distance = 66; % fix the axial distance equals to 1.0 mm
scatter_1 = [zeros(1,100),1,zeros(1,axial_distance),1,zeros(1,100)];
scatter_2 = [zeros(1,100),1,zeros(1,axial_distance),1,zeros(1,100)];
scatter_3 = [zeros(1,100),1,zeros(1,axial_distance),1,zeros(1,100)];
scatter_4 = [zeros(1,100),1,zeros(1,axial_distance),1,zeros(1,100)];
scatter_5 = [zeros(1,100),1,zeros(1,axial_distance),1,zeros(1,100)];

image_1 = conv(scatter_1, psf_1);
image_2 = conv(scatter_2, psf_2);
image_3 = conv(scatter_3, psf_3);
image_4 = conv(scatter_4, psf_4);
image_5 = conv(scatter_5, psf_5);

figure
plot(image_1)
%plot(image_2)
%plot(image_3)
%plot(image_4)
%plot(image_5)
xlabel('Axial position (mm)')
ylabel('Intensity (a.u.)')

% lateral resolution
scan_1 = log_data(85,:);   % pick five lateral lines
scan_2 = log_data(280,:);
scan_3 = log_data(480,:);
scan_4 = log_data(673,:);
scan_5 = log_data(871,:);
scatter_1_l = [zeros(1,100),1,zeros(1,70),1,zeros(1,100)];
scatter_2_l = [zeros(1,100),1,zeros(1,42),1,zeros(1,100)];
scatter_3_l = [zeros(1,100),1,zeros(1,45),1,zeros(1,100)];
scatter_4_l = [zeros(1,100),1,zeros(1,53),1,zeros(1,100)];
scatter_5_l = [zeros(1,100),1,zeros(1,68),1,zeros(1,100)];
image_1_l = conv(scatter_1_l,scan_1);
image_2_l = conv(scatter_2_l,scan_2);
image_3_l = conv(scatter_3_l,scan_3);
image_4_l = conv(scatter_4_l,scan_4);
image_5_l = conv(scatter_5_l,scan_5);

figure
plot(image_1_l)
plot(image_2_l)
plot(image_3_l)
plot(image_4_l)
plot(image_5_l)

%--------------(d)--------------%
%--------------Below is the calculation along the axial direction--------------%
x_1 = 1:1:length(psf_1);
q_1 = 1:0.1:length(psf_1);   % query point 1
ip_1 = interp1(x_1,psf_1,q_1);   % interpolation of point 1

figure
plot(x_1,psf_1,'.',q_1,ip_1,'r-');

half_1 = max(ip_1)-6;   % define the location of -6 dB
full_1 = max(ip_1)-20;   % define the location of -20 dB
first_6_db = find(ip_1 >= half_1,1);   % first point along the lateral line of -6 dB
last_6_db = find(ip_1 >= half_1,1,'last');   % last point along the lateral line of -6 dB
first_20_db = find(ip_1 >= full_1,1);
last_20_db = find(ip_1 >= full_1,1,'last');

% so the concept of calculating the width of specific dB is using last point to minus first point
width_6_db = (last_6_db-first_6_db+1)*(length(x_1)/length(q_1));   % calculate the width
width_20_db = (last_20_db-first_20_db+1)*(length(x_1)/length(q_1));

%--------------Below is the calculation along the lateral direction--------------%
la_x_1 = 1:1:length(scan_1);
la_q_1 = 1:0.1:length(scan_1);   % query point 1
la_ip_1 = interp1(la_x_1,scan_1,la_q_1);   % interpolation of point 1

%figure
%plot(la_x_1,scan_1,'.',la_q_1,la_ip_1,'r-');

la_half_1 = max(la_ip_1)-6;   % define the location of -6 dB
la_full_1 = max(la_ip_1)-20;   % define the location of -20 dB
la_first_6_db = find(la_ip_1 >= la_half_1,1);   % first point along the lateral line of -6 dB
la_last_6_db = find(la_ip_1 >= la_half_1,1,'last');   % last point along the lateral line of -6 dB
la_first_20_db = find(la_ip_1 >= la_full_1,1);   % first point along the lateral line of -20 dB
la_last_20_db = find(la_ip_1 >= la_full_1,1,'last');   % last point along the lateral line of -20 dB

% so the concept of calculating the width of specific dB is using last point to minus first point
la_width_6_db = (la_last_6_db-la_first_6_db+1)*(length(la_x_1)/length(la_q_1));   % calculate the width
la_width_20_db = (la_last_20_db-la_first_20_db+1)*(length(la_x_1)/length(la_q_1));
%------------------------------------------------------------------------------------%

%--------------(e)--------------%
x = data(400:550,60);
FFT = fft(x);
l = length(x);
P2 = abs(FFT/l);
P1 = P2(1:l/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1 = P1/max(P1);
f = fs*10^6*(0:(l/2))/l;
plot(f,P1);
title('Frequency-domain spectrum of the 60th A-line');
xlabel('f (Hz)');
ylabel('Amplitude (a.u.)');

%--------------(h)--------------%
PSF_f = data(429:528,40:80);   % define PSF at focal point 
PSF_5 = data(825:924,30:90);   % define PSF at the depth of 18 mm
N = 4500; % number of scatterers (assume there are 20 cells for each sample volume)
dx = 0.05; % define the volume of sample volume
dz = 0.02;
scale_x = 600;
scale_z = 1500;
scatterer_pos_x = rand(N,1)*scale_x; 
scatterer_pos_z = rand(N,1)*scale_z; 
Nz = 30/dz;               
Nx = 30/dx;
scatterer_dist = zeros(Nz,Nx); % spatial distribution of the scatterers, 
for iX = 1:Nx
	for iZ = 1:Nz
		if scatterer_pos_x(randi(N,1)) < scatterer_pos_z(randi(N,1))
			scatterer_dist(iZ, iX) = 1; % assign the same back-scattering coef.
		end
	end
end
Zcenter = Nz/2;
Xcenter = Nx/2;
Zrad = Nz/6;
Xrad = Nx/6;
for iiX = 1:Nx
	for iiZ = 1:Nz
            if (iiX-Xcenter)^2/Xrad^2  + (iiZ-Zcenter)^2/Zrad^2 <= 1
                scatterer_dist(iiZ, iiX) = 0;
            end
    end
end

RF_image_f = conv2(scatterer_dist, PSF_f);
RF_image_5 = conv2(scatterer_dist, PSF_5); 

env_data_f = abs(hilbert(RF_image_f));
log_data_dB_f = 20*log10(env_data_f/max(max(env_data_f))+eps)+40;
log_data_dB_f(log_data_dB_f<0) = 0;

figure
image([0 30],[0 30],log_data_dB_f);
colormap((gray(DR)));
colorbar;
title('Simulated ultrasound image of anechoic cyst at focal zone')
xlabel('Lateral position (mm)')
ylabel('Depth (mm)')
axis image

env_data_5 = abs(hilbert(RF_image_5));
log_data_dB_5 = 20*log10(env_data_5/max(max(env_data_5))+eps)+40;
log_data_dB_5(log_data_dB_5<0) = 0;

figure
image([0 30],[0 30],log_data_dB_5);
colormap((gray(DR)));
colorbar;
title('Simulated ultrasound image of anechoic cyst at depth of 18 mm')
xlabel('Lateral position (mm)')
ylabel('Depth (mm)')
axis image

figure
plot(log_data_dB_5(800,:));
hold on
plot(log_data_dB_f(800,:));
legend('18 mm','Focal point');
xlabel('Lateral position (mm)');
ylabel('Amplitude (a.u.)');
