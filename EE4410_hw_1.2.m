% EE 441000 105081041 HW1 Part_II 11/07/2019

close all;
USImageSim; 

%--------------(a)--------------%
figure
ROI_A = EnvelopeData(100:160,100:160);   % calculate the amplitude of a square region in the ROI 
histogram(ROI_A, 'FaceColor', '#31700F', 'EdgeColor', 'none');   % MATLAB recommends using the histogram function
xlabel('Amplitude (a.u.)');
ylabel('Number of pixels');
title('Probability Density Function of Amplitude in ROI');

figure
BKG_A = EnvelopeData(30:90,30:90);   % calculate the amplitude of a square region in the background
histogram(BKG_A, 'FaceColor', '#0F3C70', 'EdgeColor', 'none');
xlabel('Amplitude (a.u.)');
ylabel('Number of pixels');
title('Probability Density Function of Amplitude in background');

figure
ROI_I = (EnvelopeData(100:160,100:160)).^2;   % calculate the intensity of a square region in the ROI
histogram(ROI_I, 'FaceColor', '#31700F', 'EdgeColor', 'none');
xlabel('Intensity (a.u.)');
ylabel('Number of pixels');
title('Probability Density Function of Intensity in ROI');

figure
BKG_I = (EnvelopeData(30:90,30:90)).^2;   % calculate the intensity of a square region in the ROI
histogram(BKG_I, 'FaceColor', '#0F3C70', 'EdgeColor', 'none');
xlabel('Intensity (a.u.)');
ylabel('Number of pixels');
title('Probability Density Function of Intensity in background');

%--------------(b)--------------%
MEAN_ROI_A = mean2(ROI_A);   % calculate the mean value of rach region
MEAN_BKG_A = mean2(BKG_A);
MEAN_ROI_I = mean2(ROI_I);
MEAN_BKG_I = mean2(BKG_I);

STD_ROI_A = std2(ROI_A);   % calculate the standard deviation of each region
STD_BKG_A = std2(BKG_A);
STD_ROI_I = std2(ROI_I);
STD_BKG_I = std2(BKG_I);

SNR_ROI_A = MEAN_ROI_A ./ STD_ROI_A;
SNR_BKG_A = MEAN_BKG_A ./ STD_BKG_A;
SNR_ROI_I = MEAN_ROI_I ./ STD_ROI_I;
SNR_BKG_I = MEAN_BKG_I ./ STD_BKG_I;

%--------------(c)--------------%
dB_data_ROI = dBData(100:160,100:160);
dB_data_BKG = dBData(30:90,30:90);
STD_dB_data_ROI = std2(dB_data_ROI);
STD_dB_data_BKG = std2(dB_data_BKG);
