clear;
clc;
close all;

path_name = `Data/';

tx_distance_asphalt = importdata([path_name 'R11_1MHz.txt']).data(:,1);

%Importing Data
R11_1MHz = importdata([path_name 'R11_1MHz.txt']).data(:,2);
X11_1MHz = importdata([path_name 'X11_1MHz.txt']).data(:,2);
R12_1MHz = importdata([path_name 'R12_1MHz.txt']).data(:,2);
X12_1MHz = importdata([path_name 'X12_1MHz.txt']).data(:,2);
R21_1MHz = importdata([path_name 'R21_1MHz.txt']).data(:,2);
X21_1MHz = importdata([path_name 'X21_1MHz.txt']).data(:,2);
R22_1MHz = importdata([path_name 'R22_1MHz.txt']).data(:,2);
X22_1MHz = importdata([path_name 'X22_1MHz.txt']).data(:,2);

%Defining System Parameters
Vt = 1;

Qt = 450;
Qr = 450;

hRx_mm = 150;
hTx_mm = -tx_distance_asphalt * 10;
z_mm = 40;
d_mm = hRx_mm - hTx_mm;

z_d = z_mm ./ d_mm;
sepTx_mm = -(hTx_mm + z_mm);

%Defining Equivalent Two Port Network
Z11_1MHz = R11_1MHz + 1i * X11_1MHz;
Z12_1MHz = R12_1MHz + 1i * X12_1MHz;
Z21_1MHz = R21_1MHz + 1i * X21_1MHz;
Z22_1MHz = R22_1MHz + 1i * X22_1MHz;


Rs = 5; %Source resistance
Zt_1MHz = -imag(Z11_1MHz) .* (1i + 1 ./ Qt) + Rs;
Zr_1MHz = -imag(Z22_1MHz) .* (1i + 1 ./ Qr);

Rload_sweep = 1:1:10000;

Zin1_1MHz = zeros(length(Zr_1MHz), length(Rload_sweep));
Voc_1MHz = zeros(length(Zr_1MHz), length(Rload_sweep));
Zout2_1MHz = zeros(length(Zr_1MHz), length(Rload_sweep));

Ir_1MHz = zeros(length(Zr_1MHz), length(Rload_sweep));

Pr_1MHz = zeros(length(Zr_1MHz), length(Rload_sweep));

%Calculating optimal load resistance at 1 MHz
R_t_tot = (real(Zt_1MHz) + real(Z11_1MHz)); %real(Zt_13MHz) includes Rs
R_r_tot = (real(Zr_1MHz) + real(Z22_1MHz));
numerator_full = -real(Z12_1MHz).^2.*imag(Z12_1MHz).^2 + ...
    R_t_tot.^2.*R_r_tot.^2 - ...
    R_r_tot.*R_t_tot.*real(Z12_1MHz).^2 + ...
    R_r_tot.*R_t_tot.*imag(Z12_1MHz).^2;
root_numerator_full = sqrt(numerator_full);
R_load_optimal_1MHz = root_numerator_full./(R_t_tot);

%Calculating optimal load reactance at 1 MHz
x_num = (real(Zt_1MHz) + real(Z11_1MHz)).*(imag(Z22_1MHz)) - real(Z12_1MHz).*imag(Z12_1MHz);
x_den = real(Zt_1MHz) + real(Z11_1MHz);
X_load_rec_1MHz = -x_num./(x_den);

X_load_optimal_1MHz = X_load_rec_1MHz - imag(Zr_1MHz);

if(X_load_optimal_1MHz < 0)
    C_load_1MHz = -1/(2*pi*1e6*X_load_optimal_1MHz);
else
    L_load_1MHz = X_load_optimal_1MHz/(2*pi*1e6);
end

%Performing parameter sweep across different transmitter to asphalt seal distances
for index = 1:length(tx_distance_asphalt)

    Rload = R_load_optimal_1MHz(index);
    Zload = R_load_optimal_1MHz(index) + 1i*X_load_optimal_1MHz(index);

    Zin1_1MHz(:, index) = Zt_1MHz + Z11_1MHz - Z21_1MHz .* Z12_1MHz ./ (Zr_1MHz + Z22_1MHz + Zload);
    Voc_1MHz(:, index) = Vt .* Z21_1MHz ./ (Zt_1MHz + Z11_1MHz);
    Zout2_1MHz(:, index) = Zr_1MHz + Z22_1MHz - Z21_1MHz .* Z12_1MHz ./ (Zt_1MHz + Z11_1MHz);

    Ir_1MHz(:, index) = Voc_1MHz(:, index) ./ (Zout2_1MHz(:, index) + Zload);
    
    Pr_1MHz(:, index) = 0.5 * abs(Ir_1MHz(:, index)).^2 .* Rload;    

end

It_1MHz = Vt ./ Zin1_1MHz;

Pt_1MHz = 0.5 * real(Vt .* conj(It_1MHz));

eta_1MHz = Pr_1MHz ./ Pt_1MHz;

%Identifying optimal wireless link efficiency
[eta0_1MHz_asphalt, maxEtaIndices_1MHz] = max(eta_1MHz, [], 2);