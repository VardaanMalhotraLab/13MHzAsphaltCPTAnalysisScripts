%%
clear;
clc;
close all;

ifig = 0; %figure counter

Q = 3; %loaded quality factor
z_d_ratio_array = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1]; %ratio of lossy dielectric thickness to transmission distance

er = 5; %Lossy dielectric permitivity (real part) - aspproximate of asphalt over 0.5MHz to 15 MHz from the data
loss_t = logspace(-3,2,1000); %Loss tangent of lossy dielectric

ifig = ifig + 1;
figure(ifig)

for index=1:length(z_d_ratio_array)
    z_d_ratio = z_d_ratio_array(index);
    eta = 1./(1 + Q.*loss_t.*(z_d_ratio./(er.*(1+loss_t.^2).*(1-z_d_ratio) + z_d_ratio)));
    semilogx(loss_t, eta, 'DisplayName', sprintf("z/d = %.2f", z_d_ratio_array(index)), LineWidth=2)
    hold on
end
ax=gca;
ax.FontSize = 13; 
legend('show', fontsize = 13, Location='southwest')
xlabel('tan $\delta$','Interpreter','latex', 'FontSize', 20)
ylabel('$\eta$ (-)', 'Interpreter','latex', 'FontSize', 20)
ylim([0 1])
grid on