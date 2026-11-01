clear;
clc;
close all;

path_name = `Data/';

tx_distance_asphalt = importdata([path_name 'X11_1MHz.txt']).data(:,1);

%Importing Two-Port Parameters
R11_1MHz = importdata([path_name 'R11_1MHz.txt']).data(:,2);
X11_1MHz = importdata([path_name 'X11_1MHz.txt']).data(:,2);
R12_1MHz = importdata([path_name 'R12_1MHz.txt']).data(:,2);
X12_1MHz = importdata([path_name 'X12_1MHz.txt']).data(:,2);
R21_1MHz = importdata([path_name 'R21_1MHz.txt']).data(:,2);
X21_1MHz = importdata([path_name 'X21_1MHz.txt']).data(:,2);
R22_1MHz = importdata([path_name 'R22_1MHz.txt']).data(:,2);
X22_1MHz = importdata([path_name 'X22_1MHz.txt']).data(:,2);

R11_6MHz = importdata([path_name 'R11_6MHz.txt']).data(:,2);
X11_6MHz = importdata([path_name 'X11_6MHz.txt']).data(:,2);
R12_6MHz = importdata([path_name 'R12_6MHz.txt']).data(:,2);
X12_6MHz = importdata([path_name 'X12_6MHz.txt']).data(:,2);
R21_6MHz = importdata([path_name 'R21_6MHz.txt']).data(:,2);
X21_6MHz = importdata([path_name 'X21_6MHz.txt']).data(:,2);
R22_6MHz = importdata([path_name 'R22_6MHz.txt']).data(:,2);
X22_6MHz = importdata([path_name 'X22_6MHz.txt']).data(:,2);

R11_13MHz = importdata([path_name 'R11_13MHz.txt']).data(:,2);
X11_13MHz = importdata([path_name 'X11_13MHz.txt']).data(:,2);
R12_13MHz = importdata([path_name 'R12_13MHz.txt']).data(:,2);
X12_13MHz = importdata([path_name 'X12_13MHz.txt']).data(:,2);
R21_13MHz = importdata([path_name 'R21_13MHz.txt']).data(:,2);
X21_13MHz = importdata([path_name 'X21_13MHz.txt']).data(:,2);
R22_13MHz = importdata([path_name 'R22_13MHz.txt']).data(:,2);
X22_13MHz = importdata([path_name 'X22_13MHz.txt']).data(:,2);

%Defining system parameters
Vt = 1;

Qt = 450;
Qr = 450;

hRx_mm = 150;
hTx_mm = -tx_distance_asphalt * 10;
z_mm = 40;
d_mm = hRx_mm - hTx_mm;

z_d = z_mm ./ d_mm;
sepTx_mm = -(hTx_mm + z_mm);

%Forming Equivalent Two-Port Parameter Impedances
Z11_1MHz = R11_1MHz + 1i * X11_1MHz;
Z12_1MHz = R12_1MHz + 1i * X12_1MHz;
Z21_1MHz = R21_1MHz + 1i * X21_1MHz;
Z22_1MHz = R22_1MHz + 1i * X22_1MHz;

Z11_6MHz = R11_6MHz + 1i * X11_6MHz;
Z12_6MHz = R12_6MHz + 1i * X12_6MHz;
Z21_6MHz = R21_6MHz + 1i * X21_6MHz;
Z22_6MHz = R22_6MHz + 1i * X22_6MHz;

Z11_13MHz = R11_13MHz + 1i * X11_13MHz;
Z12_13MHz = R12_13MHz + 1i * X12_13MHz;
Z21_13MHz = R21_13MHz + 1i * X21_13MHz;
Z22_13MHz = R22_13MHz + 1i * X22_13MHz;

Zt_1MHz = -imag(Z11_1MHz) .* (1i + 1 ./ Qt);
Zr_1MHz = -imag(Z22_1MHz) .* (1i + 1 ./ Qr);

Zt_6MHz = -imag(Z11_6MHz) .* (1i + 1 ./ Qt);
Zr_6MHz = -imag(Z22_6MHz) .* (1i + 1 ./ Qr);

Zt_13MHz = -imag(Z11_13MHz) .* (1i + 1 ./ Qt);
Zr_13MHz = -imag(Z22_13MHz) .* (1i + 1 ./ Qr);

Rload_sweep = 1:1:10000;

%Defining variables for FOR loop used to identify RL that leads to largest wireless link efficiency
Zin1_1MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Voc_1MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Zout2_1MHz = zeros(length(Zr_13MHz), length(Rload_sweep));

Zin1_6MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Voc_6MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Zout2_6MHz = zeros(length(Zr_13MHz), length(Rload_sweep));

Zin1_13MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Voc_13MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Zout2_13MHz = zeros(length(Zr_13MHz), length(Rload_sweep));

Ir_1MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Ir_6MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Ir_13MHz = zeros(length(Zr_13MHz), length(Rload_sweep));

Pr_1MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Pr_6MHz = zeros(length(Zr_13MHz), length(Rload_sweep));
Pr_13MHz = zeros(length(Zr_13MHz), length(Rload_sweep));

%Loop to identify RL that leads to largest wireless link efficiency
for index = 1:length(Rload_sweep)

    Rload = Rload_sweep(index);

    Zin1_1MHz(:, index) = Zt_1MHz + Z11_1MHz - Z21_1MHz .* Z12_1MHz ./ (Zr_1MHz + Z22_1MHz + Rload);
    Voc_1MHz(:, index) = Vt .* Z21_1MHz ./ (Zt_1MHz + Z11_1MHz);
    Zout2_1MHz(:, index) = Zr_1MHz + Z22_1MHz - Z21_1MHz .* Z12_1MHz ./ (Zt_1MHz + Z11_1MHz);
    
    Zin1_6MHz(:, index) = Zt_6MHz + Z11_6MHz - Z21_6MHz .* Z12_6MHz ./ (Zr_6MHz + Z22_6MHz + Rload);
    Voc_6MHz(:, index) = Vt .* Z21_6MHz ./ (Zt_6MHz + Z11_6MHz);
    Zout2_6MHz(:, index) = Zr_6MHz + Z22_6MHz - Z21_6MHz .* Z12_6MHz ./ (Zt_6MHz + Z11_6MHz);
    
    Zin1_13MHz(:, index) = Zt_13MHz + Z11_13MHz - Z21_13MHz .* Z12_13MHz ./ (Zr_13MHz + Z22_13MHz + Rload);
    Voc_13MHz(:, index) = Vt .* Z21_13MHz ./ (Zt_13MHz + Z11_13MHz);
    Zout2_13MHz(:, index) = Zr_13MHz + Z22_13MHz - Z21_13MHz .* Z12_13MHz ./ (Zt_13MHz + Z11_13MHz);
    
    Ir_1MHz(:, index) = Voc_1MHz(:, index) ./ (Zout2_1MHz(:, index) + Rload);
    Ir_6MHz(:, index) = Voc_6MHz(:, index) ./ (Zout2_6MHz(:, index) + Rload);
    Ir_13MHz(:, index) = Voc_13MHz(:, index) ./ (Zout2_13MHz(:, index) + Rload);
    
    Pr_1MHz(:, index) = 0.5 * abs(Ir_1MHz(:, index)).^2 .* Rload;
    Pr_6MHz(:, index) = 0.5 * abs(Ir_6MHz(:, index)).^2 .* Rload;
    Pr_13MHz(:, index) = 0.5 * abs(Ir_13MHz(:, index)).^2 .* Rload;    

end

It_1MHz = Vt ./ Zin1_1MHz;
It_6MHz = Vt ./ Zin1_6MHz;
It_13MHz = Vt ./ Zin1_13MHz;

Pt_1MHz = 0.5 * real(Vt .* conj(It_1MHz));
Pt_6MHz = 0.5 * real(Vt .* conj(It_6MHz));
Pt_13MHz = 0.5 * real(Vt .* conj(It_13MHz));

eta_1MHz = Pr_1MHz ./ Pt_1MHz;
eta_6MHz = Pr_6MHz ./ Pt_6MHz;
eta_13MHz = Pr_13MHz ./ Pt_13MHz;

%Identify max wireless link efficiency
[eta0_1MHz, maxEtaIndices_1MHz] = max(eta_1MHz, [], 2);
[eta0_6MHz_asphalt, maxEtaIndices_6MHz] = max(eta_6MHz, [], 2);
[eta0_13MHz, maxEtaIndices_13MHz] = max(eta_13MHz, [], 2);

%Approximation for RL,opt that was found to be valid for low-loss dielectrics
Rload0_1MHz_approx = real(Zr_1MHz + Z22_1MHz) .* sqrt(1 + imag(Z21_1MHz).^2 ./ (real(Zr_1MHz + Z22_1MHz).*real(Zt_1MHz + Z11_1MHz)));
Rload0_6MHz_approx = real(Zr_6MHz + Z22_6MHz) .* sqrt(1 + imag(Z21_6MHz).^2 ./ (real(Zr_6MHz + Z22_6MHz).*real(Zt_6MHz + Z11_6MHz)));
Rload0_13MHz_approx = real(Zr_13MHz + Z22_13MHz) .* sqrt(1 + imag(Z21_13MHz).^2 ./ (real(Zr_13MHz + Z22_13MHz).*real(Zt_13MHz + Z11_13MHz)));