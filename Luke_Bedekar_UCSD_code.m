%% Multivariate time-dependent models for multiple events
% with corrected titration values for events with high AUC
% Rayanne Luke 2026
% Using data from the Pia Lab (UCSD)

set(groot, 'defaultLineLineWidth', 2, 'defaulttextinterpreter', 'latex', ...
    'defaultAxesFontSize', 20, 'defaultLegendInterpreter', 'latex')

data = readtable('../../Datasets/UCSD/sero_ped_days_vac.xlsm');

days_inf = data.days_af_infection;
days_vax1 = data.days_af_1st_vac;
days_vax2 = data.days_af_2nd_vac;
days_btwn_inf_vax2 = data.relative_event_days;
S_AUC = log2(data.spikeigg_auc);
cc = find(~S_AUC);
S_AUC(cc) = [];
days_inf(cc) = []; days_vax1(cc) = []; 
days_vax2(cc) = []; days_btwn_inf_vax2(cc) = [];
data(cc,:) = [];

inf_label = logical(strcmp(data.conf_nat_inf, 'inf_only'));
vax_label = logical(strcmp(data.conf_nat_inf, 'vac_only'));
inf_vax_label = logical(strcmp(data.conf_nat_inf, 'inf_vac'));
data_inf = data(inf_label, :);
data_vax = data(vax_label, :);
data_inf_vax = data(inf_vax_label, :);
nan_vax2 = ~isnan(data_inf_vax.days_af_2nd_vac);
data_inf_vax = data_inf_vax(nan_vax2, :);

G = findgroups(data.id);
G_inf = findgroups(data_inf.id);
G_vax = findgroups(data_vax.id);
G_inf_vax = findgroups(data_inf_vax.id);

G_cross = []; % group indices for crossover events
for i = 1:G(end)
    l = find(G == i);
    data_l = data(l,:);
    cr_i = sum(logical(strcmp(data_l.conf_nat_inf, 'inf_only')));
    cr_iv = sum(logical(strcmp(data_l.conf_nat_inf, 'inf_vac')));
    if (cr_i > 0 && cr_iv > 0)
        G_cross = [G_cross; i];
    end
end

I_days = days_inf(strcmp(data.conf_nat_inf, 'inf_only'));
I_AUC = S_AUC(strcmp(data.conf_nat_inf, 'inf_only'));

V2_days = days_vax2(strcmp(data.conf_nat_inf, 'vac_only'));
V_AUC = S_AUC(strcmp(data.conf_nat_inf, 'vac_only'));

IV_days = days_inf(strcmp(data.conf_nat_inf, 'inf_vac'));
IV_days_v2 = days_vax2(strcmp(data.conf_nat_inf, 'inf_vac'));
IV_c = find(isnan(IV_days_v2));
IV_AUC = S_AUC(strcmp(data.conf_nat_inf, 'inf_vac'));
IV_days(IV_c) = []; IV_days_v2(IV_c) = [];
IV_AUC(IV_c) = []; 

rel_days = data.relative_event_days;
time_btwn_I_V = rel_days(strcmp(data.conf_nat_inf, 'inf_vac'));
time_btwn_I_V(IV_c) = [];

%% Plotting raw data

h1t = tiledlayout('flow');
nexttile; hold on
splitapply(@myplot, I_days, I_AUC, G_inf)
ax = gca; ax.ColorOrderIndex = 2;
h12 = plot(I_days, I_AUC, 'x', 'markersize', 8);
xlabel('Days post infection', 'FontSize', 24)
legend(h12, 'Infection', 'fontsize', 20, 'location', 'southeast')
xlim([0, 400]), ylim([0, 16])
yticks([0, 4, 8, 12, 16])
hold off

nexttile
hold on
splitapply(@myplot, V2_days, V_AUC, G_vax)
ax = gca; ax.ColorOrderIndex = 1;
h11 = plot(V2_days, V_AUC, 'o', 'markersize', 8);
xlabel('Days post vaccination', 'FontSize', 24)
legend(h11, 'Vaccination', 'fontsize', 20, 'location', 'southeast')
xlim([0, 250]), ylim([0, 16])
yticks([0, 4, 8, 12, 16])
hold off

nexttile
hold on
splitapply(@myplot, IV_days_v2, IV_AUC, G_inf_vax)
ax = gca; ax.ColorOrderIndex = 5;
h13 = plot(IV_days_v2, IV_AUC, '^', 'markersize', 8);
xlabel('Days post vax after infection', 'FontSize', 24)
legend(h13, 'Infection then vaccination', 'fontsize', 20, 'location', 'southeast')
xlim([0, 250]), ylim([0, 16])
yticks([0, 4, 8, 12, 16])
hold off

ylabel(h1t, 'Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')

%% Titration-extrapolation

data_t = readtable('../../Datasets/UCSD/sero_ped_days_vac_titers.xlsx');
G_t = findgroups(data_t.id);
for i = 1:G_t(end)
    l = find(G_t == i);
    data_l = data_t(l,:);
    data_l = sortrows(data_l, [3, 4]);
    data_t(l,:) = data_l;
end
S_raw_t = data_t.spikeigg_auc;
S_AUC_t = log2(S_raw_t);
c_t = find(~S_AUC_t);
data_t(c_t, :) = [];
vax_label_t = logical(strcmp(data_t.conf_nat_inf, 'vac_only') | ...
    (strcmp(data_t.conf_nat_inf, 'inf_vac')));
inf_label_t = logical(strcmp(data_t.conf_nat_inf, 'inf_only'));
vax1_days_t = data_t.vax1_rel_days;
vax2_days_t = cell2mat(cellfun(@(x) str2double(x), data_t.vax2_rel_days, 'un', 0));
vax_auc = data_t.spikeigg_auc;

t_100 = data_t.Sp100;
t_300 = data_t.Sp300;
t_900 = data_t.Sp900;
t_2700 = data_t.Sp2700;
t_8100 = data_t.Sp8100;

t_matrix = [t_100 t_300 t_900 t_2700 t_8100];

t_label = ~isnan(t_matrix); t_label = t_label(:,1);
[r, c] = find(t_matrix < 0);
for i = 1:length(c)
    if c == 5
        t_matrix(r(i),c(i)) = -t_matrix(r(i), c(i)-1);
    else
       t_matrix(r(i), c(i):end) = -t_matrix(r(i), c(i)-1);
    end
end
OD_dec = find(all(diff(t_matrix')<0));
t_ph = ~ismember(1:length(t_label), OD_dec);
t_idx = OD_dec;
t_label(t_ph) = 0;

t_matrix_dec = t_matrix(OD_dec, :);
vax_auc_dec = vax_auc(OD_dec);

dil_vec = [100 300 900 2700 8100];
dil_vec2 = [dil_vec 3*8100]; % extend another one titration levels
levels = 1;

auc_test_calc = trapz(dil_vec, t_matrix_dec');
figure; plot(auc_test_calc, vax_auc_dec, 'o')
mdl = fitlm(auc_test_calc, vax_auc_dec);
hold on, plot(mdl), 
xlabel('AUC calculated', 'interpreter', 'latex'), 
ylabel('AUC reported', 'interpreter', 'latex'), title('')

shift_val = mdl.Coefficients.Estimate(1);
auc_calc = auc_test_calc + shift_val;
ll = find(auc_calc < 0); mm = find(vax2_days_t < 0);
lm = [ll'; mm];
auc_calc(lm) = []; vax2_days_t(lm) = []; t_matrix_dec(lm, :) = []; 
t_idx(lm) = []; 
lk = find(t_idx == 260 | t_idx == 263);
t_idx(lk) = []; auc_calc(lk) = []; t_matrix_dec(lk, :) = [];
% ^ don't inflate only 2/4 of a subject's meas.
t_label = zeros(length(t_label),1);
t_label(t_idx) = 1;

dil_vec = dil_vec'; t_matrix_dec = t_matrix_dec'; 

ext_vals = zeros(size(t_matrix_dec, levels), levels);
for i = 1:size(t_matrix_dec, 2)
    f = fit(log(dil_vec), log(t_matrix_dec(:,i)), 'poly2')';
    ext_vals(i, :) = exp(f(log(dil_vec2(end-(levels-1):end))));
end

tit_matrix_dec2 = [t_matrix_dec; ext_vals'];
figure;
loglog(dil_vec2(1:end-levels), tit_matrix_dec2(1:end-levels, :), '-o')
hold on, ax = gca; ax.ColorOrderIndex = 1;
loglog(dil_vec2(end-levels:end), tit_matrix_dec2(end-levels:end, :), '--o')
xlabel('Dilution ratio')
ylabel('Optical density'), legend('All samples')

auc_test_calc2 = trapz(dil_vec2, tit_matrix_dec2);
auc_calc2 = auc_test_calc2 + shift_val;

%% Read in data, plot personal trajectories (inf only, vax only, inf then vax)

neg_dist = [18.2, 0.152];

data = readtable('../../Datasets/UCSD/sero_ped_days_vac.xlsm');
G_d = findgroups(data.id);
for i = 1:G_d(end)
    l = find(G_d == i);
    data_l = data(l,:);
    data_l = sortrows(data_l, [3, 8]);
    data(l,:) = data_l;
end

days_inf = data.days_af_infection;
days_vax1 = data.days_af_1st_vac;
days_vax2 = data.days_af_2nd_vac;
days_btwn_inf_vax2 = data.relative_event_days;
S_raw = data.spikeigg_auc; S_raw(c_t) = [];
S_raw_old = data.spikeigg_auc; S_raw_old(c_t) = [];
S_raw(t_idx) = auc_calc2; % altered AUC values changed here!!!
S_AUC = log2(S_raw);
days_inf(c_t) = []; days_vax1(c_t) = []; 
days_vax2(c_t) = []; days_btwn_inf_vax2(c_t) = [];
data(c_t,:) = [];
data.S_raw_corr = S_raw;
data.S_AUC_log2 = S_AUC;

inf_label = logical(strcmp(data.conf_nat_inf, 'inf_only'));
vax_label = logical(strcmp(data.conf_nat_inf, 'vac_only'));
inf_vax_label = logical(strcmp(data.conf_nat_inf, 'inf_vac'));
data_inf = data(inf_label, :);
data_vax = data(vax_label, :);
data_inf_vax = data(inf_vax_label, :);
nan_vax2 = ~isnan(data_inf_vax.days_af_2nd_vac);
data_inf_vax = data_inf_vax(nan_vax2, :);

G = findgroups(data.id);
G_inf = findgroups(data_inf.id);
G_vax = findgroups(data_vax.id);
G_inf_vax = findgroups(data_inf_vax.id);

I_data_struct = struct([]); V_data_struct = struct([]);
for i = 1:G_inf(end)
    I_data_struct{i} = data_inf(G_inf == i, :);
end
for i = 1:G_vax(end)
    V_data_struct{i} = data_vax(G_vax == i, :);
end

G_cross = []; % group indices for crossover events
for i = 1:G(end)
    l = find(G == i);
    data_l = data(l,:);
    cr_i = sum(logical(strcmp(data_l.conf_nat_inf, 'inf_only')));
    cr_iv = sum(logical(strcmp(data_l.conf_nat_inf, 'inf_vac')));
    if (cr_i > 0 && cr_iv > 0)
        G_cross = [G_cross; i];
    end
end

I_days = days_inf(strcmp(data.conf_nat_inf, 'inf_only'));
I_AUC = S_AUC(strcmp(data.conf_nat_inf, 'inf_only'));

V2_days = days_vax2(strcmp(data.conf_nat_inf, 'vac_only'));
V_AUC = S_AUC(strcmp(data.conf_nat_inf, 'vac_only'));

IV_days = days_inf(strcmp(data.conf_nat_inf, 'inf_vac'));
IV_days_v2 = days_vax2(strcmp(data.conf_nat_inf, 'inf_vac'));
IV_c = find(isnan(IV_days_v2));
IV_AUC = S_AUC(strcmp(data.conf_nat_inf, 'inf_vac'));
IV_days(IV_c) = []; IV_days_v2(IV_c) = [];
IV_AUC(IV_c) = []; 

%remove entries where IV_days_v2 are negative
indices_neg = IV_days_v2<=0;
IV_days_v2(indices_neg) = [];
IV_AUC(indices_neg) = [];
G_inf_vax(indices_neg) = [];
IV_days(indices_neg) = [];


rel_days = data.relative_event_days;
time_btwn_I_V = rel_days(strcmp(data.conf_nat_inf, 'inf_vac'));
time_btwn_I_V(IV_c) = [];
time_btwn_I_V(indices_neg) = [];


idxs = find(ismember(G, G_cross));
cross_data = data(idxs, :);
cross_days = days_inf(idxs, :);
cross_S_AUC = S_AUC(idxs, :);
GGG = findgroups(cross_data.id);

I_days_c = cross_days(strcmp(cross_data.conf_nat_inf, 'inf_only'));
I_AUC_c = cross_S_AUC(strcmp(cross_data.conf_nat_inf, 'inf_only'));
cs = find(~I_AUC_c); I_AUC_c(cs) = []; I_days_c(cs) = [];

IV_AUC_c = cross_S_AUC(strcmp(cross_data.conf_nat_inf, 'inf_vac'));
IV_days_c = cross_days(strcmp(cross_data.conf_nat_inf, 'inf_vac'));

%% Plotting titration-extrapolated data

h1t = tiledlayout('flow');
nexttile; hold on
splitapply(@myplot, I_days, I_AUC, G_inf)
ax = gca; ax.ColorOrderIndex = 2;
h12 = plot(I_days, I_AUC, 'x', 'markersize', 8);
xlabel('Days post infection', 'FontSize', 24)
legend(h12, 'Infection', 'fontsize', 20, 'location', 'southeast')
xlim([0, 400]), ylim([0, 16])
yticks([0, 4, 8, 12, 16])
hold off

nexttile
hold on
splitapply(@myplot, V2_days, V_AUC, G_vax)
ax = gca; ax.ColorOrderIndex = 1;
h11 = plot(V2_days, V_AUC, 'o', 'markersize', 8);
xlabel('Days post vaccination', 'FontSize', 24)
legend(h11, 'Vaccination', 'fontsize', 20, 'location', 'southeast')
xlim([0, 250]), ylim([0, 16])
yticks([0, 4, 8, 12, 16])
hold off

nexttile
hold on
splitapply(@myplot, IV_days_v2, IV_AUC, G_inf_vax)
ax = gca; ax.ColorOrderIndex = 5;
h13 = plot(IV_days_v2, IV_AUC, '^', 'markersize', 8);
xlabel('Days post vax after infection', 'FontSize', 24)
legend(h13, 'Infection then vaccination', 'fontsize', 20, 'location', 'southeast')
xlim([0, 250]), ylim([0, 16])
yticks([0, 4, 8, 12, 16])
hold off

ylabel(h1t, 'Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')

%% Probability models

x = linspace(0, 4 + ceil(max([I_AUC; V_AUC; IV_AUC])), 1000);
t = linspace(0, max([I_days; V2_days; IV_days]) + 6, 1000);
[x, t] = meshgrid(x,t); x = x'; t = t';

t_delay = 150;
k1 = 1.15;
k2 = 1.2;

t_scale_fac = 100; % nondimensionalize t using 100 as characteristic time

obj_I = @(p) -sum(log(time_dep(I_AUC, I_days/t_scale_fac, [p(1), p(2), neg_dist], p(3))));
obj_V = @(p) -sum(log(time_dep(V_AUC, V2_days/t_scale_fac, [p(1), p(2), neg_dist], p(3))));
   
[opt_I, fval_I, exit_flagI] = fmincon(obj_I, [700, 33, 1.2], ...
    [], [], [], [], [sqrt(eps), sqrt(eps), 1 + sqrt(eps)]);
inf_fit = @(x,t) time_dep(x, t/t_scale_fac, [opt_I(1:2), neg_dist], opt_I(3));

[opt_V, fval_V, exit_flagV] = fmincon(obj_V, [700, 33, 1.2], ...
    [], [], [], [], [sqrt(eps), sqrt(eps), 1 + sqrt(eps)]);
vax_fit = @(x,t) time_dep(x, t/t_scale_fac, [opt_V(1:2), neg_dist], opt_V(3));

obj_IV = @(p) -sum(log(time_dep_two(IV_AUC, time_btwn_I_V/t_scale_fac, ...
    IV_days_v2/t_scale_fac, [opt_I(1:2), p(1), p(2), neg_dist], opt_I(3),  p(3))));

[opt_IV, fval_IV, exit_flagIV] = fmincon(obj_IV, [700, 33, 1.2], ...
    [], [], [], [], [sqrt(eps), sqrt(eps), 1 + sqrt(eps)]);

IV_fit = @(x, t) time_dep_two_plot(x, t/t_scale_fac, t_delay/t_scale_fac, ...
    [opt_I(1:2), opt_IV(1:2), neg_dist], opt_I(3), opt_IV(3));

contI = [0.4 0.3 0.2 0.1 0.01 0.001]; 
contV = contI;
        
        V_mat = vax_fit(x, t);
        I_mat = inf_fit(x, t);
        IV_mat = IV_fit(x, t);

prob_vals_I = zeros(length(I_data_struct), 5); % max # elements in a personal trajectory is 5
prob_vals_V = zeros(length(V_data_struct), 4); % max # elements in a personal trajectory is 4

for k = 1:length(I_data_struct)
    prob_vals_I(k, 1:size(I_data_struct{k},1)) ...
        = inf_fit(I_data_struct{k}.S_AUC_log2, I_data_struct{k}.days_af_infection)';
end
for k = 1:length(V_data_struct)
    prob_vals_V(k, 1:size(V_data_struct{k},1)) ...
        = vax_fit(V_data_struct{k}.S_AUC_log2, V_data_struct{k}.days_af_2nd_vac)';
end 

low_prob_vals_I = sum((prob_vals_I > 0) & (prob_vals_I < 0.001), 2);
high_prob_vals_I = sum((prob_vals_I > 0) & (prob_vals_I > 0.2), 2);
high_prob_I_perc = high_prob_vals_I./sum(prob_vals_I ~= 0, 2);
low_prob_vals_V = sum((prob_vals_V > 0) & (prob_vals_V < 0.001), 2);
high_prob_vals_V = sum((prob_vals_V > 0) & (prob_vals_V > 0.2), 2);
high_prob_V_perc = high_prob_vals_V./sum(prob_vals_V ~= 0, 2);

%% Plotting probability models with data and personal trajectories

rng(0)
N_AUC = random('Gamma', neg_dist(1), neg_dist(2), 250, 1);

figure; 
hh = tiledlayout('flow');


  nexttile 
histogram(N_AUC,'Normalization','pdf', 'facecolor', [211,211,211]/256)
hold on, plot(x, gampdf(x, neg_dist(1), neg_dist(2)), 'k')
xlabel('Log-scaled antibody meas.', 'FontSize', 18)
ylabel('PDF', 'Fontsize', 16)
legend('Na{\"i}ve', 'fontsize', 20)
xlim([0, 8])

 nexttile
 contourf(t, x, I_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.25);
 ax1 = gca;
 colormap(ax1, 'summer'), %set(gca,'ColorScale','log')
 clim([0, 0.3])
 hold on
  contour(t, x, I_mat, contI, 'linewidth', 1.5, ...
    'linecolor', [0.2,0.2,0.2])
  splitapply(@myplot, I_days, I_AUC, G_inf)
 ax = gca; ax.ColorOrderIndex = 2;
 hi = scatter(I_days, I_AUC, 25, 'x', 'MarkerFaceColor','flat', 'linewidth', 3);
    legend(hi,  'Infection', 'fontsize', 20)
    xlim([0, 400])

    nexttile 
contourf(t, x, V_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.2);
ax0 = gca;
colormap(ax0, 'parula'), %set(gca,'ColorScale','log')
clim([0, 0.3])
hold on
 contour(t, x, V_mat, contV, 'linewidth', 1.5, ...
    'linecolor', [0.2,0.2,0.2])
 splitapply(@myplot, V2_days, V_AUC, G_vax)
 ax0.ColorOrderIndex = 1;
 hv = scatter(V2_days, V_AUC, 30, 'd', 'MarkerFaceColor','flat');
 legend(hv, 'Vaccination', 'fontsize', 20)
 xlim([0, 250])

nexttile
 contourf(t, x, IV_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.2);
 hold on
contour(t, x, IV_mat, contI, 'linewidth', 1.5, ...
    'linecolor', [0.2,0.2,0.2])
ax = gca; 
colormap(ax, 'autumn'), %set(gca,'ColorScale','log')
clim([0, 0.3])
splitapply(@myplot, days_inf, S_AUC, G)
ax.ColorOrderIndex = 2;
hii =  scatter(I_days, I_AUC, 25, 'x', 'MarkerFaceColor','flat', 'linewidth', 3);
 ax.ColorOrderIndex = 5;
 hiv = scatter(IV_days, IV_AUC, 30, '^', 'MarkerFaceColor','flat');
  legend([hii, hiv], 'Infection', ...
      'Infection then vaccination', 'fontsize', 20, 'location', 'southeast')

xlabel(hh, 'Days post infection or vaccination', 'FontSize', 24, 'interpreter', 'latex')
ylabel(hh, 'Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')

%% Look at a few personal trajectories of interest of vaccination
% Representative trajectories and outliers

V_interest = [10, 13, 14, 32]; 
leg_V = cell(length(V_interest), 1);
clear hv
figure; hold on
contourf(t, x, V_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.2);
    contour(t, x, V_mat, contI, 'linewidth', 1.5, ...
        'linecolor', [0.2,0.2,0.2], 'ShowText','on')
for k = 1:length(V_interest)
    data_k = V_data_struct{V_interest(k)};
    days_vax2_k = data_k.days_af_2nd_vac;
    S_AUC_k = data_k.S_AUC_log2;
    ax = gca;  colormap(ax, 'parula'), clim([0,0.3])
    plot(days_vax2_k, S_AUC_k, 'k')
    hv(k) =  scatter(days_vax2_k, S_AUC_k, 70, 'd', 'MarkerFaceColor','flat');
    leg_V = {'Subject D2-7', 'Subject D2-8', 'Subject D2-9', 'Subject D2-10'};
    yticks([0, 4, 8, 12, 16]), ylim([0, 20])
    xlabel('Days post first dose', 'FontSize', 24, 'interpreter', 'latex')
    ylabel('Antibody meas. (AUC)', 'FontSize', 24, 'interpreter', 'latex')
    xlim([0,250])
end
legend(hv, leg_V, 'fontsize', 20, 'location', 'southeast')

%% Look at a few personal trajectories of interest of infection
% Representative trajectories and outliers

I_interest = [1, 5, 24, 27, 66, 105];
leg_I = cell(length(I_interest), 1);
clear hi
figure; hold on
contourf(t, x, I_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.25);
    contour(t, x, I_mat, contI, 'linewidth', 1.5, ...
        'linecolor', [0.596, 0.486, 0.463], 'ShowText','on')
for k = 1:length(I_interest)
    data_k = I_data_struct{I_interest(k)};
    days_inf_k = data_k.days_af_infection;
    S_AUC_k = data_k.S_AUC_log2;
    ax = gca;  colormap(ax, 'summer'), clim([0,0.7])
    plot(days_inf_k, S_AUC_k, 'k')
    hi(k) =  scatter(days_inf_k, S_AUC_k, 70, 'x', ...
        'MarkerFaceColor','flat', 'linewidth', 3);
    leg_I = {'Subject D2-1', 'Subject D2-2', 'Subject D2-3', 'Subject D2-4', 'Subject D2-5', 'Subject D2-6'};
    yticks([0, 4, 8, 12, 16])
    xlabel('Days post infection', 'FontSize', 24, 'interpreter', 'latex')
    ylabel('Antibody meas. (AUC)', 'FontSize', 24, 'interpreter', 'latex')
    xlim([0,400]), ylim([0, 20])
end
legend(hi, leg_I, 'fontsize', 20, 'location', 'southeast')

%% Identify crossover events

% Crossover events w/ personal trajectories

idxs = find(ismember(G, G_cross));
cross_data = data(idxs, :);
cross_days = days_inf(idxs, :);
cross_S_AUC = S_AUC(idxs, :);
GGG = findgroups(cross_data.id);
cross_data_struct = struct([]); cross_rel_days = zeros(GGG(end),1);
for i = 1:GGG(end)
    cross_data_struct{i} = cross_data(GGG == i, :);
    neg_vax2_day_idx = find(cross_data_struct{i}.days_af_2nd_vac < 0);
    if ~isempty(neg_vax2_day_idx)
        cross_data_struct{1,i}(neg_vax2_day_idx,:) = [];
    end
    cross_rel_days(i) = min(cross_data_struct{i}.relative_event_days);
end

[cross_rel_days_sort, idx] = sort(cross_rel_days);

I_days_c = cross_days(strcmp(cross_data.conf_nat_inf, 'inf_only'));
I_AUC_c = cross_S_AUC(strcmp(cross_data.conf_nat_inf, 'inf_only'));
cs = find(~I_AUC_c); I_AUC_c(cs) = []; I_days_c(cs) = [];

IV_AUC_c = cross_S_AUC(strcmp(cross_data.conf_nat_inf, 'inf_vac'));
IV_days_c = cross_days(strcmp(cross_data.conf_nat_inf, 'inf_vac'));
idxs_IV = find(ismember(IV_AUC, IV_AUC_c));

%% Look at a few crossover trajectories of interest

IV_interest = [25, 40, 42, 55,  61];
figure; 
hh = tiledlayout('flow'); 
for k = 1:length(IV_interest)
    data_k = cross_data_struct{idx(IV_interest(k))};
    days_inf_k = data_k.days_af_infection;
    days_vax1_k = data_k.days_af_1st_vac;
    days_vax2_k = data_k.days_af_2nd_vac;
    days_btwn_inf_vax2_k = data_k.relative_event_days;
    S_AUC_k = data_k.S_AUC_log2;
    
    inf_label_k = logical(strcmp(data_k.conf_nat_inf, 'inf_only'));
    vax_label_k = logical(strcmp(data_k.conf_nat_inf, 'vac_only'));
    inf_vax_label_k = logical(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    data_inf_k = data_k(inf_label_k, :);
    data_vax_k = data_k(vax_label_k, :);
    data_inf_vax_k = data_k(inf_vax_label_k, :);
    nan_vax2_k = ~isnan(data_inf_vax_k.days_af_2nd_vac);
    data_inf_vax_k = data_inf_vax_k(nan_vax2_k, :);

    I_days_k = days_inf_k(strcmp(data_k.conf_nat_inf, 'inf_only'));
    I_AUC_k = S_AUC_k(strcmp(data_k.conf_nat_inf, 'inf_only'));
    
    V_AUC_k = S_AUC_k(strcmp(data_k.conf_nat_inf, 'vac_only'));
    
    IV_days_k = days_inf_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    IV_days_v2_k = days_vax2_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    IV_days_v1_k = days_vax1_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    IV_c_k = find(isnan(IV_days_v1_k));
    IV_AUC_k = S_AUC_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    IV_days_k(IV_c_k) = []; 
    
    %remove entries where IV_days_v2 are negative
    indices_neg_k = IV_days_v2_k < 0;
    if ~isempty(find(indices_neg_k, 1))
        s_idx = find(S_AUC_k == IV_AUC_k(indices_neg_k));
        S_AUC_k(s_idx) = [];
        days_inf_k(s_idx) = [];
        IV_AUC_k(indices_neg_k) = [];
        IV_days_k(indices_neg_k) = [];
    end
    
    rel_days_k = data_k.relative_event_days;
    time_btwn_I_V_k = rel_days_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    time_btwn_I_V_k(IV_c_k) = [];
    time_btwn_I_V_k(indices_neg_k) = [];
    time_btwn_I_V_k = time_btwn_I_V_k(1);


        IV_fit_individ = @(x, t) time_dep_two_plot(x, t/t_scale_fac, ...
            time_btwn_I_V_k/t_scale_fac, [opt_I(1:2), opt_IV(1:2), neg_dist], opt_I(3), opt_IV(3));

        IV_mat_individ = IV_fit_individ(x, t);

    nexttile, hold on
    contourf(t, x, IV_mat_individ, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.25)
    contour(t, x, IV_mat_individ, contI, 'linewidth', 1.5, ...
        'linecolor', [0.2, 0.2, 0.2], 'showtext', 'on')
    ax = gca; colormap(ax, 'autumn'), clim([0, 0.7])
    plot(days_inf_k, S_AUC_k, 'k')
    plot(ones(100,1)*time_btwn_I_V_k, linspace(0, 20, 100), 'w--')
    text(time_btwn_I_V_k-30, 14, 'Vaccination', 'color', 'w', 'fontsize', 16, 'Rotation', 90)
    ax.ColorOrderIndex = k+3;
    h2  = scatter(I_days_k, I_AUC_k, 50, 'x', 'markerfacecolor', 'flat', 'linewidth', 3);
    ax = gca; ax.ColorOrderIndex = k+3;
    leg_IV = {'Subject D2-11', 'Subject D2-12', 'Subject D2-13', 'Subject D2-14', 'Subject D2-15'};
    h3 = scatter(IV_days_k, IV_AUC_k, 50, '^', 'markerfacecolor', 'flat');
        legend(h2, leg_IV, 'fontsize', 20, ...
        'location', 'southeast');
    yticks([0, 4, 8, 12, 16]) 
end
 xlabel(hh, 'Days post infection', 'FontSize', 24, 'interpreter', 'latex')
    ylabel(hh, 'Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')

%% Individual models for crossover events

file_name = 'UCSD_two_event_model.gif';

prob_vals = zeros(length(cross_data_struct), 6); % max # elements in a personal trajectory is 6
for k = 1:length(cross_data_struct)
    zoom off
    drawnow
    data_k = cross_data_struct{idx(k)};
    days_inf_k = data_k.days_af_infection;
    days_vax1_k = data_k.days_af_1st_vac;
    days_vax2_k = data_k.days_af_2nd_vac;
    days_btwn_inf_vax2_k = data_k.relative_event_days;
    S_AUC_k = data_k.S_AUC_log2;
    
    inf_label_k = logical(strcmp(data_k.conf_nat_inf, 'inf_only'));
    vax_label_k = logical(strcmp(data_k.conf_nat_inf, 'vac_only'));
    inf_vax_label_k = logical(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    data_inf_k = data_k(inf_label_k, :);
    data_vax_k = data_k(vax_label_k, :);
    data_inf_vax_k = data_k(inf_vax_label_k, :);
    nan_vax2_k = ~isnan(data_inf_vax_k.days_af_2nd_vac);
    data_inf_vax_k = data_inf_vax_k(nan_vax2_k, :);

    I_days_k = days_inf_k(strcmp(data_k.conf_nat_inf, 'inf_only'));
    I_AUC_k = S_AUC_k(strcmp(data_k.conf_nat_inf, 'inf_only'));
    
    V_AUC_k = S_AUC_k(strcmp(data_k.conf_nat_inf, 'vac_only'));
    
    IV_days_k = days_inf_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    IV_days_v2_k = days_vax2_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    IV_days_v1_k = days_vax1_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    IV_c_k = find(isnan(IV_days_v1_k));
    IV_AUC_k = S_AUC_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    IV_days_k(IV_c_k) = []; 
    
    %remove entries where IV_days_v2 are negative
    indices_neg_k = IV_days_v2_k < 0;
    if ~isempty(find(indices_neg_k, 1))
        s_idx = find(S_AUC_k == IV_AUC_k(indices_neg_k));
        S_AUC_k(s_idx) = [];
        days_inf_k(s_idx) = [];
        IV_AUC_k(indices_neg_k) = [];
        IV_days_k(indices_neg_k) = [];
    end
    
    rel_days_k = data_k.relative_event_days;
    time_btwn_I_V_k = rel_days_k(strcmp(data_k.conf_nat_inf, 'inf_vac'));
    time_btwn_I_V_k(IV_c_k) = [];
    time_btwn_I_V_k(indices_neg_k) = [];
    time_btwn_I_V_k = time_btwn_I_V_k(1);

        IV_fit_individ = @(x, t) time_dep_two_plot(x, t/t_scale_fac, ...
            time_btwn_I_V_k/t_scale_fac, [opt_I(1:2), opt_IV(1:2), neg_dist], opt_I(3), opt_IV(3));

        IV_mat_individ = IV_fit_individ(x, t);

    prob_vals(k, 1:length([I_AUC_k; IV_AUC_k])) = IV_fit_individ([I_AUC_k; IV_AUC_k], [I_days_k; IV_days_k]);    

    figure;
    contourf(t, x, IV_mat_individ, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.2);
    hold on
    contour(t, x, IV_mat_individ, contI, 'linewidth', 1.5, ...
        'linecolor', [0.596, 0.486, 0.463])
    ax = gca; ax.ColorOrderIndex = 2;  colormap(ax, 'autumn')
      scatter(I_days_k, I_AUC_k, 10, 'MarkerFaceColor','flat')
     ax.ColorOrderIndex = 5;
      scatter(IV_days_k, IV_AUC_k, 10, 'MarkerFaceColor','flat')
      plot(ones(100,1)*time_btwn_I_V_k, linspace(0, 20, 100), 'w--')
    text(time_btwn_I_V_k-20, 14, 'Vaccination', 'color', 'w', 'fontsize', 18, 'Rotation', 90)
    plot(days_inf_k, S_AUC_k, 'k')
    ax = gca; ax.ColorOrderIndex = k+3;
    h2  = scatter(I_days_k, I_AUC_k, 50, 'x', 'markerfacecolor', 'flat', 'linewidth', 3);
    ax = gca; ax.ColorOrderIndex = k+3;
    h3 = scatter(IV_days_k, IV_AUC_k, 50, '^', 'markerfacecolor', 'flat');
    yticks([0, 4, 8, 12, 16])
    xlabel('Days post infection', 'FontSize', 24, 'interpreter', 'latex')
    ylabel('Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')
    exportgraphics(gcf,file_name,'Append',true);
end

%% Functions

function myplot(dm,d)
plot(dm,d,'k', 'linewidth', 0.5);
end

function shapetpos = shapetpos(t, params, expt) 
%gives the shape parameter for positive distribution (R) at timeperiod t.

shape = params(3);
theta = params([1 2]);
shapetpos = (theta(1).*t)./(1+(theta(2).*t.^expt)) + shape;

end

function td_fit = time_dep(x, t, p, expt)

a = shapetpos(t, p(1:3), expt);
td_fit = gampdf(x, a, p(4));

end

function W = time_dep_two(r, t_inf_vac, t_vax_meas, p, exptI, exptIV)
% two-event distribution

shape_1 = shapetpos(t_inf_vac, [p(1), p(2), p(5)], exptI);
shape_2 = shapetpos(t_vax_meas, [p(3), p(4), p(5)], exptIV);
shape = shape_1 + shape_2 - p(5);
scale = p(6);
W = gampdf(r, shape, scale);

end

function W = time_dep_two_plot(r, t, t_delay, p, exptI, exptIV)
% two-event distribution
%code by Joseph Nakao!
shape_1 = shapetpos(t, [p(1), p(2), p(5)], exptI);

t_since = t - t_delay;
shape_2 = shapetpos(t_since, [p(3), p(4), p(5)], exptIV);
shape_2(t < t_delay) = 0;
shape = shape_1 + shape_2;
shape(t >= t_delay) = shape(t >= t_delay) - p(5);

scale = p(6);
W = gampdf(r, shape, scale);

end
