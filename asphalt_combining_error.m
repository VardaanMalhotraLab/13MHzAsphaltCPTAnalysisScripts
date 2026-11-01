%% Asphalt Slab Cores Measurement Uncertainties
clear
clc
close all

%Errors due to VNA
e_p_small_core_std = load("asphalt_small_core_ep_std.mat").e_p_core_std;
e_p_large_core1_std = load("asphalt_large_core1_ep_std.mat").e_p_core1_std;
e_p_large_core2_std = load("asphalt_large_core2_ep_std.mat").e_p_core2_std;

e_p_VNA_std = sqrt((e_p_small_core_std.^2 + e_p_large_core1_std.^2 + e_p_large_core2_std.^2)./3);

loss_t_small_core_std = load("asphalt_small_core_losst_std.mat").loss_t_core_std;
loss_t_large_core1_std = load("asphalt_large_core1_losst_std.mat").loss_t_core1_std;
loss_t_large_core2_std = load("asphalt_large_core2_losst_std.mat").loss_t_core2_std;

losst_VNA_std = sqrt((loss_t_small_core_std.^2 + loss_t_large_core1_std.^2 + loss_t_large_core2_std.^2)./3);

%Errors due to asphlt dimension variations


f = load("frequency.mat").f';

%Smaller Asphalt Core (~ 60 mm Thickness)
e_p_small_core_avg_r71p35_t57p68 = load("asphalt_small_core_er_r71p35_t57p68.mat").e_p_core_avg;
e_p_small_core_avg_r72p1_t57p8 = load("asphalt_small_core_er_r72p1_t57p8.mat").e_p_core_avg;
e_p_small_core_avg_r72p85_t57p92 = load("asphalt_small_core_er_r72p85_t57p92.mat").e_p_core_avg;

loss_t_small_core_avg_r71p35_t57p68 = load("asphalt_small_core_losst_r71p35_t57p68.mat").loss_t_core_avg;
loss_t_small_core_avg_r72p1_t57p8 = load("asphalt_small_core_losst_r72p1_t57p8.mat").loss_t_core_avg;
loss_t_small_core_avg_r72p85_t57p92 = load("asphalt_small_core_losst_r72p85_t57p92.mat").loss_t_core_avg;

%Larger Asphalt Core 1 (~ 75 mm Thickness)
e_p_large_core1_avg_r71p25_t71p48 = load("asphalt_large_core1_er_r71p25_t71p48.mat").e_p_core1_avg;
e_p_large_core1_avg_r72_t71p6 = load("asphalt_large_core1_er_r72_t71p6.mat").e_p_core1_avg;
e_p_large_core1_avg_r72p75_t71p72 = load("asphalt_large_core1_er_r72p75_t71p72.mat").e_p_core1_avg;

loss_t_large_core1_avg_r71p25_t71p48 = load("asphalt_large_core1_losst_r71p25_t71p48.mat").loss_t_core1_avg;
loss_t_large_core1_avg_r72_t71p6 = load("asphalt_large_core1_losst_r72_t71p6.mat").loss_t_core1_avg;
loss_t_large_core1_avg_r72p75_t71p72 = load("asphalt_large_core1_losst_r72p75_t71p72.mat").loss_t_core1_avg;


%Larger Asphalt Core 2 (~ 75 mm Thickness)
e_p_large_core2_avg_r71p15_t72p18= load("asphalt_large_core2_er_r71p15_t72p18.mat").e_p_core2_avg;
e_p_large_core2_avg_r71p9_t72p3 = load("asphalt_large_core2_er_r71p9_t72p3.mat").e_p_core2_avg;
e_p_large_core2_avg_r72p65_t72p42 = load("asphalt_large_core2_er_r72p65_t72p42.mat").e_p_core2_avg;

loss_t_large_core2_avg_r71p15_t72p18 = load("asphalt_large_core2_losst_r71p15_t72p18.mat").loss_t_core2_avg;
loss_t_large_core2_avg_r71p9_t72p3 = load("asphalt_large_core2_losst_r71p9_t72p3.mat").loss_t_core2_avg;
loss_t_large_core2_avg_r72p65_t72p42 = load("asphalt_large_core2_losst_r72p65_t72p42.mat").loss_t_core2_avg;



%Pooled standard deviations (relative permittivity)
min_ep_val_deviation = ( (e_p_small_core_avg_r72p85_t57p92(5:32:end) - e_p_small_core_avg_r72p1_t57p8(5:32:end)) + (e_p_large_core1_avg_r72p75_t71p72(5:32:end) - e_p_large_core1_avg_r72_t71p6(5:32:end)) + (e_p_large_core2_avg_r72p65_t72p42(5:32:end) - e_p_large_core2_avg_r71p9_t72p3(5:32:end)) ) ./3;
min_ep_val_deviation = sqrt(min_ep_val_deviation.^2 ...
    + (2*e_p_VNA_std(5:32:end)).^2 );

max_ep_val_deviation = ( (e_p_small_core_avg_r71p35_t57p68(5:32:end) - e_p_small_core_avg_r72p1_t57p8(5:32:end)) + (e_p_large_core1_avg_r71p25_t71p48(5:32:end) - e_p_large_core1_avg_r72_t71p6(5:32:end)) + (e_p_large_core2_avg_r71p15_t72p18(5:32:end) - e_p_large_core2_avg_r71p9_t72p3(5:32:end)) )./3;
max_ep_val_deviation = sqrt(max_ep_val_deviation.^2 ...
    + (2*e_p_VNA_std(5:32:end)).^2 );

%Pooled standard deviations (loss tangent)
min_losst_val_deviation = ( (loss_t_small_core_avg_r72p85_t57p92(5:32:end) - loss_t_small_core_avg_r72p1_t57p8(5:32:end)) + (loss_t_large_core1_avg_r72p75_t71p72(5:32:end) - loss_t_large_core1_avg_r72_t71p6(5:32:end)) + (loss_t_large_core2_avg_r72p65_t72p42(5:32:end) - loss_t_large_core2_avg_r71p9_t72p3(5:32:end)) )./3;
min_losst_val_deviation = sqrt(min_losst_val_deviation.^2 ...
    + (2*losst_VNA_std(5:32:end)).^2 );

max_losst_val_deviation = ( (loss_t_small_core_avg_r71p35_t57p68(5:32:end) - loss_t_small_core_avg_r72p1_t57p8(5:32:end)) + (loss_t_large_core1_avg_r71p25_t71p48(5:32:end) - loss_t_large_core1_avg_r72_t71p6(5:32:end)) + (loss_t_large_core2_avg_r71p15_t72p18(5:32:end) - loss_t_large_core2_avg_r71p9_t72p3(5:32:end)) )./3;
max_losst_val_deviation = sqrt(max_losst_val_deviation.^2 ...
    + (2*losst_VNA_std(5:32:end)).^2 );


