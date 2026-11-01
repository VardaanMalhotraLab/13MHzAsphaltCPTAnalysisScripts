%% Asphalt Electrical Properties
% Script to determine the electrical properties of asphalt sample using differential mixed mode parameters.
%Reference paper for formulas:
% "Dielectric measurement by open-ended coaxial line for hot-mix asphalt roads: From laboratory test to on-site investigation" by Fan et al.
% "Design and Calibration of a Large Open-Ended  Coaxial Probe for the  Measurement of  the Dielectric Properties of Concrete" by Filali et al.
% "Modeling, simulation and experimental validation of solid media in capacitive wireless power transfer" by Lecluyse et al.

clear
clc
close all


file_count = 10;

%Parameter definitions
Y0 = 1/100;
e0 = 8.854e-12; 
A_rexolite = pi*(72/1000)^2;
A_core1 = pi*(72/1000)^2;
A_core2 = pi*(71.9/1000)^2;
A_air =  A_rexolite; %the sample with the largest diameter to be set equal to A_air
d_total = 75/1000;
d_air = d_total;
d_rexolite = 60/1000;
d_core1 = 71.6/1000;
d_core2 = 72.3/1000;

e_air = 1.00059;



path_name = 'Data/';

f = importdata([path_name 'short1.csv']).data(:,1);
w = 2.*pi.*f;


%%%%%%   ---REXOLITE---   %%%%%%
rexolite_freq = [1.00e4, 1.00e7, 1.00e8, 3.00e8, 1.00e9, 3.00e9, 8.50e9, 1.40e10, 2.40e10];
rexoltie_er = [2.531, 2.531, 2.531, 2.53, 2.53, 2.529, 2.529, 2.529, 2.528];
rexolite_tan_delta = [0.00012, 0.00025, 0.00047, 0.00053, 0.00061, 0.00063, 0.00066, 0.00067, 0.00068];


rexolite_er_interpolated = pchip(rexolite_freq,rexoltie_er,f);
rexolite_tan_delta_interpolated = pchip(rexolite_freq,rexolite_tan_delta,f);

e_Rexolite = rexolite_er_interpolated;
e_dp_Rexolite = e_Rexolite .* rexolite_tan_delta_interpolated;

%%%%%%   ---EXPECTED REF. COEFF. ---   %%%%%%

% g1 = -1; %(Short)
cond = 5.96e7; %copper conductivity
a = 3.22e-4; %AWG 22 wire radius in meters
short_length = 0.1; %meters
mu_0 = 1.2566e-6;
skindepth_short = 1./(sqrt(pi.*f.*mu_0.*cond));
R_short = short_length./(cond.*pi.*(a^2 - (a - skindepth_short).^2));
L_short = round(2*short_length*100*(log((2*short_length*100)/(a*100)) - 1 + mu_0/4)); %output in nanohenries if input is in centimeters
X_short = 2.*pi.*f.*L_short.*1e-9; %L_short output is in nanohenries
Z_short = R_short + 1i*(X_short);
Z0 = 100;
g1 = (Z_short - Z0) ./ (Z_short + Z0); %(Short)
% g1= -0.9819 + 1i*0.1842; %(Short) @ 13.56 MHz
g2 = (Y0 - (1i.*w.*(e_air.*e0.*A_air./d_air))) ./ (Y0 + (1i.*w.*(e_air.*e0.*A_air./d_air)));  %Air
Y_rexolite = (w.*(e_dp_Rexolite.*e0.*A_rexolite./d_rexolite) + (1i.*w.*(e_Rexolite.*e0.*A_rexolite./d_rexolite)));
Y_rexolite = Y_rexolite + (1i.*w.*(e_air.*e0.*(A_air-A_rexolite)./(d_total)));
Y_rexolite_air = (1i.*w.*(e_air.*e0.*A_rexolite./(d_total - d_rexolite))); %Accounting for air gap between rexolite surface and capacitor plate
Y_rexolite_eq = (Y_rexolite.*Y_rexolite_air)./(Y_rexolite + Y_rexolite_air);
g3 = (Y0 - Y_rexolite_eq) ./ (Y0 + Y_rexolite_eq); %Lossy Rexolite


e_p_core1 = zeros(file_count, length(f));
e_dp_core1 = zeros(file_count, length(f));
loss_t_core1 = zeros(file_count, length(f));

e_p_core2 = zeros(file_count, length(f));
e_dp_core2 = zeros(file_count, length(f));
loss_t_core2 = zeros(file_count, length(f));

%Obtaining the S-parameters from each measurement dataset and converting it to an equivfalent sample admittance (and hence capacitance and resistance)
for num = 1:file_count

    p1 = importdata([path_name 'short' num2str(num) '.csv']).data(:,2) + 1i*importdata([path_name 'short' num2str(num) '.csv']).data(:,3);
    z1 = s2z(p1(357), 100);
    p2 = importdata([path_name 'air' num2str(num) '.csv']).data(:,2) + 1i*importdata([path_name 'air' num2str(num) '.csv']).data(:,3);
    z2 = s2z(p2(357), 100);
    p3 = importdata([path_name 'rexolite' num2str(num) '.csv']).data(:,2) + 1i*importdata([path_name 'rexolite' num2str(num) '.csv']).data(:,3);
    z3 = s2z(p3(357), 100);

    S11 = (g1.*g2.*p3.*(p1 - p2) + g1.*g3.*p2.*(p3 - p1) + g2.*g3.*p1.*(p2 - p3)) ./ (g1.*g2.*(p1 - p2) + g1.*g3.*(p3 - p1) + g2.*g3.*(p2 - p3));
    S22 = (g1.*(p2 - S11) + g2.*(S11 - p1)) ./ (g1.*g2.*(p2 - p1));
    S12_S21 = (p1 - S11).*(1 - S22.*g1) ./ (g1);

    pc1 = importdata([path_name 'core1_' num2str(num) '.csv']).data(:,2) + 1i*importdata([path_name 'core1_' num2str(num) '.csv']).data(:,3);
    zc1 = s2z(pc1(357), 100);

    pc2 = importdata([path_name 'core2_' num2str(num) '.csv']).data(:,2) + 1i*importdata([path_name 'core2_' num2str(num) '.csv']).data(:,3);
    zcq = s2z(pc2(357), 100);

    %%%%%%   ---LARGE SLAB CORE 1---   %%%%%%

    gc1 = (pc1 - S11) ./ (S22.*pc1 + S12_S21 - S11.*S22);
    Ys_eq1 = Y0*((1 - gc1)./(1 + gc1)); 
    
    Ycore1 = Ys_eq1 - 1i.*(w.*(e_air*e0.*(A_air-A_core1)./(d_total))); %removing parallel air "coating" around dielectric
    Ys_air1 = 1i.*(w.*(e_air*e0.*A_core1./(d_total - d_core1)));
    Ycore1_isolated = (Ys_air1.*Ycore1)./(Ys_air1 - Ycore1); %removing effect of series "air" dielectric between capacitor plate and surface of sample
    
    Cc1 = imag(Ycore1_isolated)./(w);
    Rc1 = 1./real(Ycore1_isolated);
    
    Cvacuum_core1 = e0.*A_core1/(d_core1);
    e_p_core1(num,:) = Cc1./Cvacuum_core1;
    e_dp_core1(num,:) = 1./(w.*Rc1.*Cvacuum_core1);
    loss_t_core1(num,:) = e_dp_core1(num,:)./e_p_core1(num,:);

    %%%%%%   ---LARGE SLAB CORE 2---   %%%%%%

    gc2 = (pc2 - S11) ./ (S22.*pc2 + S12_S21 - S11.*S22);
    Ys_eq2 = Y0*((1 - gc2)./(1 + gc2)); 
    
    Ycore2 = Ys_eq2 - 1i.*(w.*(e_air*e0.*(A_air-A_core2)./(d_total))); %removing parallel air "coating" around dielectric
    Ys_air2 = 1i.*(w.*(e_air*e0.*A_core2./(d_total - d_core2)));
    Ycore2_isolated = (Ys_air2.*Ycore2)./(Ys_air2 - Ycore2); %removing effect of series "air" dielectric between capacitor plate and surface of sample
    
    Cc2 = imag(Ycore2_isolated)./(w);
    Rc2 = 1./real(Ycore2_isolated);
    
    Cvacuum_core2 = e0.*A_core2/(d_core2);
    e_p_core2(num,:) = Cc2./Cvacuum_core2;
    e_dp_core2(num,:) = 1./(w.*Rc2.*Cvacuum_core2);
    loss_t_core2(num,:) = e_dp_core2(num,:)./e_p_core2(num,:);

end

%Calculating averages across datasets
e_p_core1_avg = mean(e_p_core1, 1);
e_dp_core1_avg = mean(e_dp_core1, 1);
loss_t_core1_avg = mean(loss_t_core1, 1);

e_p_core2_avg = mean(e_p_core2, 1);
e_dp_core2_avg = mean(e_dp_core2, 1);
loss_t_core2_avg = mean(loss_t_core2, 1);

%Calculating standard deviation across datasets to characterize noise/deviation from the VNA
e_p_core1_std = std(e_p_core1, 0, 1);
loss_t_core1_std = std(loss_t_core1, 0, 1);

e_p_core2_std = std(e_p_core2, 0, 1);
loss_t_core2_std = std(loss_t_core2, 0, 1);
