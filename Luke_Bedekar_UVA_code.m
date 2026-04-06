%% Multivariate time-dependent models for multiple events
% Rayanne Luke 2026
% Using data from the Long COVID Recovery Cohort and Vaccine Cohort (UVA)

set(groot, 'defaultLineLineWidth', 2, 'defaulttextinterpreter', 'latex', ...
    'defaultAxesFontSize', 24, 'defaultLegendInterpreter', 'latex')

%% Read in and pre-process data
% LC = long COVID cohort
% VC = vaccine cohort

data_LC_all = readtable('../../Datasets/UVA/UVA_antibody_data.xlsx', ...
    'sheet', 'Timedep_meas_ONLY');
data_VC = readtable('../../Datasets/UVA/UVA_antibody_data_vax_cohort_days_since.xlsx');

data_LC = data_LC_all(~isnan(data_LC_all.Days_Start_Inf), :);

inf_label_LC = logical(isnan(data_LC.Days_LastVax) ...
    & ~isnan(data_LC.Spike_IgG) & ~isnan(data_LC.Days_Start_Inf));
inf_vax_label_LC = logical(~isnan(data_LC.Days_Vax_1) &  ...
    isnan(data_LC.Days_Vax_2) & ~strcmp(data_LC.Days_Start_Inf, 'NA') ...
    & ~isnan(data_LC.Spike_IgG) & (data_LC.Days_Start_Inf > data_LC.Days_LastVax));
inf_vax_vax_label_LC = logical(~isnan(data_LC.Days_Vax_2) & ...
     isnan(data_LC.Days_Vax_3) & ~strcmp(data_LC.Days_Start_Inf, 'NA') ...
     & ~isnan(data_LC.Spike_IgG));

inf_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Inf') ...
    & ~isnan(data_VC.Spike_IgG) & data_VC.Days_Inf > 0); 
inf_vax_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Inf_Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax2_dose_2));
inf_vax_vax_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Inf_Vax_Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax3_booster));
inf_first_label_VC = inf_label_VC | inf_vax_label_VC | inf_vax_vax_label_VC;

vax_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax2_dose_2));
vax_inf_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax_Inf') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Inf));
vax_vax_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax_Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax3_booster));
vax_vax_inf_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax_Vax_Inf') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Inf));
vax_only_label_VC = vax_label_VC | vax_vax_label_VC;
vax_first_label_VC = vax_label_VC | vax_vax_label_VC | vax_inf_label_VC ...
    | vax_vax_inf_label_VC;

days_since_inf_LC = data_LC.Days_Start_Inf;

days_inf_LC = data_LC.Days_Start_Inf(inf_label_LC);
days_inf_vax_LC = data_LC.Days_Vax_1(inf_vax_label_LC); 
days_inf_vax_vax_LC = data_LC.Days_Vax_2(inf_vax_vax_label_LC); 

days_since_VC = data_VC.Days_since_first_event;
days_since_VC(vax_first_label_VC) = data_VC.Days_Vax2_dose_2(vax_first_label_VC);

days_since_inf_VC = data_VC.Days_Inf;
days_inf_VC = data_VC.Days_Inf(inf_label_VC);
days_inf_vax_VC = data_VC.Days_Vax2_dose_2(inf_vax_label_VC);
days_inf_vax_vax_VC = data_VC.Days_Vax3_booster(inf_vax_vax_label_VC);

days_since_vax_VC = data_VC.Days_Vax2_dose_2; 
days_vax_VC = data_VC.Days_Vax2_dose_2(vax_label_VC);
days_vax_vax_VC = data_VC.Days_Vax3_booster(vax_vax_label_VC);
days_vax_vax_inf_VC = data_VC.Days_Inf(vax_vax_inf_label_VC);

rel_days_inf_vax_LC = data_LC.Days_Start_Inf(inf_vax_label_LC) - days_inf_vax_LC;
rel_days_inf_vax_LC_IVV = data_LC.Days_Start_Inf(inf_vax_vax_label_LC)...
    - data_LC.Days_Start_Vax(inf_vax_vax_label_LC);
rel_days_inf_vax_vax_LC = data_LC.Days_Start_Vax(inf_vax_vax_label_LC) ...
    - days_inf_vax_vax_LC;

rel_days_vax_vax_VC = rmmissing(data_VC.Days_Vax2_dose_2(vax_vax_label_VC))...
    - days_vax_vax_VC;
ng_VC = find(rel_days_vax_vax_VC < 0 | isnan(rel_days_vax_vax_VC));
rel_days_vax_vax_VC(ng_VC) = [];
days_vax_vax_VC(ng_VC) = [];
rel_days_vax_vax_VC_VVI = data_VC.Days_Vax2_dose_2(vax_vax_inf_label_VC)...
    - data_VC.Days_Vax3_booster(vax_vax_inf_label_VC);
rel_days_vax_vax_inf_VC = data_VC.Days_Vax3_booster(vax_vax_inf_label_VC) ...
    - days_vax_vax_inf_VC;
rel_days_inf_vax_VC = data_VC.Days_since_first_event(inf_vax_label_VC) ...
    - days_inf_vax_VC;
rel_days_inf_vax_VC_IVV = data_VC.Days_since_first_event(inf_vax_vax_label_VC) ...
    - data_VC.Days_Vax2_dose_2(inf_vax_vax_label_VC);
ng_VC2 = find(rel_days_inf_vax_VC < 0 | isnan(rel_days_inf_vax_VC));
rel_days_inf_vax_VC(ng_VC2) = [];
days_inf_vax_VC(ng_VC2) = [];
rel_days_inf_vax_vax_VC = data_VC.Days_Vax1_dose_1(inf_vax_vax_label_VC) ...
    - days_inf_vax_vax_VC;

days_vax_only_VC = days_since_vax_VC(vax_only_label_VC);

S_LC = log2(data_LC.Spike_IgG + 2) - 1; 
% ^ log-scale and shift to avoid negative values
S_neg_LC = data_LC_all.Spike_IgG(strcmp(data_LC_all.Visit_ID, 'Base'));
S_inf_LC = S_LC(inf_label_LC);
S_inf_vax_LC = S_LC(inf_vax_label_LC);
S_inf_vax_vax_LC = S_LC(inf_vax_vax_label_LC); 

S_VC = log2(data_VC.Spike_IgG + 2) - 1; 
S_nans = find(isnan(S_VC));
S_vax_only_VC = S_VC(vax_only_label_VC);

S_inf_VC = S_VC(inf_label_VC);
S_vax_VC = S_VC(vax_label_VC); 
S_vax_vax_VC = S_VC(vax_vax_label_VC);
S_vax_vax_VC(ng_VC) = [];
S_inf_vax_VC = S_VC(inf_vax_label_VC);
S_inf_vax_VC(ng_VC2) = [];
S_inf_vax_vax_VC = S_VC(inf_vax_vax_label_VC);
S_vax_vax_inf_VC = S_VC(vax_vax_inf_label_VC);

S_neg_VC = data_VC.Spike_IgG(data_VC.Days_Vax1_dose_1 <= 0 ...
    & isnan(data_VC.Days_Inf));
S_neg = [S_neg_LC; S_neg_VC];
S_neg = log2(S_neg + 2) - 1;

G_LC = findgroups(data_LC.Subject_ID);
data_vax_LC = data_LC(~inf_label_LC, :);
G_vax_LC = findgroups(data_vax_LC.Subject_ID);
data_inf_LC = data_LC(inf_label_LC, :);
G_inf_LC = findgroups(data_inf_LC.Subject_ID);

G_VC = findgroups(data_VC.Subject_ID);
data_inf_VC = data_VC(inf_first_label_VC, :);
data_vax_VC = data_VC(vax_first_label_VC, :);
data_vax1_VC = data_VC(vax_label_VC, :);
data_vax_only_VC = data_VC(vax_only_label_VC, :);
G_vax_VC = findgroups(data_vax_VC.Subject_ID);
G_vax1_VC = findgroups(data_vax1_VC.Subject_ID);
G_vax_only_VC = findgroups(data_vax_only_VC.Subject_ID);
G_inf_VC = findgroups(data_inf_VC.Subject_ID);

data_inf_VC_days = data_inf_VC.Days_since_first_event;
data_inf_VC_S = log2(data_inf_VC.Spike_IgG+2)-1;

data_vax_VC_days = data_vax_VC.Days_Vax2_dose_2;
data_vax_VC_S = log2(data_vax_VC.Spike_IgG+2)-1;

data_struct_LC = struct([]); 
for i = 1:G_LC(end)
    data_struct_LC{i} = data_LC(G_LC == i, :);
end

I_data_struct_LC = struct([]); 
for i = 1:G_inf_LC(end)
    I_data_struct_LC{i} = data_inf_LC(G_inf_LC == i, :);
end

IV_data_struct_LC = struct([]);
data_inf_vax_LC = [];
k = 1;
G_inf_vax_LC = []; % group indices for crossover events
for i = 1:G_LC(end)
    l = find(G_LC == i);
    data_l = data_LC(l,:);
    cr_i = sum(logical(strcmp(data_l.Vaxed, 'n')));
    cr_iv = sum(logical(strcmp(data_l.Vaxed, 'y')));
    cr_reld = diff(data_l.Days_Start_Inf);
    if (cr_i > 0 && cr_iv > 0 && all(cr_reld > 0))
        data_l = data_l((strcmp(data_l.Vaxed, 'n') ...
            | strcmp(data_l.Vaxed, 'y')) & data_l.Days_Start_Inf > 0 ...
           & ~isnan(data_l.Spike_IgG) & isnan(data_l.Days_Vax_2), :);
        data_inf_vax_LC = [data_inf_vax_LC; data_l];
        IV_data_struct_LC{k} = data_l;
        k = k + 1;
    end
end

 G_inf_vax_LC = findgroups(data_inf_vax_LC.Subject_ID);

IV_data_struct_VC = struct([]);
data_inf_vax_VC = [];
G_inf_vax_VC = []; % group indices for crossover events
k = 1;
for i = 1:G_VC(end)
    l = find(G_VC == i);
    data_l = data_VC(l,:);
    cr_i = sum(logical(strcmp(data_l.Inf_Vax_traj, 'Inf')));
    cr_iv = sum(logical(strcmp(data_l.Inf_Vax_traj, 'Inf_Vax')));
    cr_reld = diff(data_l.Days_since_first_event);
    if ((cr_i > 0 || cr_iv > 0))
        data_l = data_l(((strcmp(data_l.Inf_Vax_traj, 'Inf') ...
            | strcmp(data_l.Inf_Vax_traj, 'Inf_Vax')) & ...
            ~isnan(data_l.Spike_IgG) & data_l.Days_Inf > 0 ...
            & ~isnan(data_l.Days_Vax2_dose_2)), :);
        if ~isempty(data_l)
            data_inf_vax_VC = [data_inf_vax_VC; data_l];
            IV_data_struct_VC{k} = data_l;
            k = k + 1;
        end
    end
end

G_inf_vax_VC = findgroups(data_inf_vax_VC.Subject_ID);

V_data_struct_VC = struct([]); 
for i = 1:G_vax_VC(end)
    V_data_struct_VC{i} = data_vax_VC(G_vax_VC == i, :);
end

VV_data_struct_VC = struct([]);
data_vax_vax_VC = [];
G_vax_vax_VC = []; % group indices for crossover events
k = 1;
for i = 1:G_VC(end)
    l = find(G_VC == i);
    data_l = data_VC(l,:);
    cr_i = sum(logical(strcmp(data_l.Inf_Vax_traj, 'Vax')));
    cr_iv = sum(logical(strcmp(data_l.Inf_Vax_traj, 'Vax_Vax')));
    cr_reld = diff(data_l.Days_since_first_event);
    if ((cr_i > 0 && cr_iv > 0))
        data_l = data_l(((strcmp(data_l.Inf_Vax_traj, 'Vax') ...
            | strcmp(data_l.Inf_Vax_traj, 'Vax_Vax')) & ...
            ~isnan(data_l.Spike_IgG) ...
            & ~isnan(data_l.Days_Vax2_dose_2)), :);
        if size(data_l,1) > 1
            data_vax_vax_VC = [data_vax_vax_VC; data_l];
            VV_data_struct_VC{k} = data_l;
            k = k + 1;
        end
    end
end

VVI_data_struct_VC = struct([]);
k = 1;
data_vax_vax_inf = [];
for i = 1:G_vax_VC(end)
    temp = data_vax_VC(G_vax_VC == i, :);
    ct = length(find(strcmp(temp.Inf_Vax_traj, 'Vax_Vax_Inf')));
    if ct > 0
        data_vax_vax_inf = [data_vax_vax_inf; temp];
        VVI_data_struct_VC{k} = temp;
        k = k + 1;
    end
end

G_vax_inf = findgroups(data_vax_vax_inf.Subject_ID);

V_only_data_struct_VC = struct([]); 
for i = 1:G_vax_only_VC(end)
    V_only_data_struct_VC{i} = data_vax_only_VC(G_vax_only_VC == i, :);
end

I_data_struct_VC = struct([]); 
for i = 1:G_inf_VC(end)
    I_data_struct_VC{i} = data_inf_VC(G_inf_VC == i, :);
end

%% Probability models for infection and subsequent events

exptI = 2;
exptIV = 1.25;
exptIVV = 1.7;

neg_dist = fitdist(S_neg, 'Gamma');

days_inf = [days_inf_LC; days_inf_VC];
days_inf_vax = [days_inf_vax_LC; days_inf_vax_VC];
days_inf_vax_vax = [days_inf_vax_vax_LC; days_inf_vax_vax_VC];
S_inf = [S_inf_LC; S_inf_VC];
S_inf_vax = [S_inf_vax_LC; S_inf_vax_VC];
S_inf_vax_vax = [S_inf_vax_vax_LC; S_inf_vax_vax_VC];
rel_days_inf_vax = [rel_days_inf_vax_LC; rel_days_inf_vax_VC];
rel_days_inf_vax_IVV = [rel_days_inf_vax_LC_IVV; rel_days_inf_vax_VC_IVV];
rel_days_inf_vax_vax = [rel_days_inf_vax_vax_LC; rel_days_inf_vax_vax_VC];

x = linspace(0, 12, 1000);
t = linspace(0, 700, 1000);
[x, t] = meshgrid(x,t); x = x'; t = t';

t_scale_fac = 100; % nondimensionalize t using 100 as characteristic time

obj_I = @(p) -sum(log(time_dep(S_inf, days_inf/t_scale_fac, [p(1), p(2), ...
    neg_dist.a, neg_dist.b], p(3))));
    
options_I = optimoptions('fmincon', 'Algorithm','sqp');
[opt_I, I_fval, exit_flag_I] = fmincon(obj_I, [5, 5, 0.5], [], [], ...
    [], [], [sqrt(eps), sqrt(eps), 1 + sqrt(eps)],[], [], options_I);
inf_fit = @(x,t) time_dep(x, t/t_scale_fac, ...
    [opt_I(1:2), neg_dist.a, neg_dist.b], opt_I(end));

obj_IV = @(p) -sum(log(time_dep_two(S_inf_vax, rel_days_inf_vax/t_scale_fac,...
    days_inf_vax/t_scale_fac, [opt_I(1:2), p(1), p(2), neg_dist.a, ...
    neg_dist.b], opt_I(end), p(3))));
% ^ both exptI, exptIV vary as parameters

[opt_IV, IV_fval, exit_flag_IV] = fmincon(obj_IV, [5, 5, 0.5], [], [], ...
    [], [], [sqrt(eps), sqrt(eps), 1 + sqrt(eps)], [], [], options_I);
t_delay = 100;

IV_fit = @(x, t) time_dep_two_plot(x, t/t_scale_fac, t_delay/t_scale_fac, ...
    [opt_I(1:2), opt_IV(1:2), neg_dist.a, neg_dist.b], opt_I(end), opt_IV(end));
IV_fit_2t = @(x, t, t_delay) time_dep_two_plot(x, t/t_scale_fac, t_delay/t_scale_fac, ...
    [opt_I(1:2), opt_IV(1:2), neg_dist.a, neg_dist.b], opt_I(end), opt_IV(end));

rel_days_inf_vax_vax([2,4,5]) = [];
days_inf_vax_vax([2,4,5]) = [];
rel_days_inf_vax_IVV([2,4,5]) = [];
S_inf_vax_vax([2,4,5]) = [];

obj_IVV = @(p) -sum(log(time_dep_three(S_inf_vax_vax, ...
    rel_days_inf_vax_IVV/t_scale_fac, rel_days_inf_vax_vax/t_scale_fac, ...
    days_inf_vax_vax/t_scale_fac, [opt_I(1:2), opt_IV(1:2), ...
    p(1), p(2), neg_dist.a, neg_dist.b], opt_I(end), opt_IV(end), p(3))));
% ^ exptI, exptIV, exptIVV vary as parameters

options_IVV = optimset('Algorithm', 'sqp');
init_IVV = [0.05, 0.1, exptIVV];
[opt_IVV, IVV_fval, exit_flag1] = fmincon(obj_IVV, init_IVV, ...
   [], [], [], [], [sqrt(eps), sqrt(eps), 1+sqrt(eps)], ...
   [Inf, Inf, Inf], [], options_IVV);
t_delay2 = 250;

IVV_fit = @(x, t) time_dep_three_plot(x, t/t_scale_fac, t_delay/t_scale_fac, ...
    t_delay2/t_scale_fac, [opt_I(1:2), opt_IV(1:2), opt_IVV(1:2), ...
    neg_dist.a, neg_dist.b], opt_I(end), opt_IV(end), opt_IVV(end));

contI = [0.6 0.5 0.4 0.3 0.2 0.1 0.01, 0.001]; 

        I_mat = inf_fit(x, t);
        IV_mat = IV_fit(x, t);
        IVV_mat = IVV_fit(x, t);

%% Probability models for vaccination and boosting

exptV = 2;
exptVV = 2;

neg_dist2 = fitdist(S_neg, 'Gamma');

x2 = linspace(0, 12, 1000);
t2 = linspace(0, 750, 1000); 

[x2, t2] = meshgrid(x2,t2); x2 = x2'; t2 = t2';

obj_V = @(p) -sum(log(time_dep(S_vax_VC, days_vax_VC/t_scale_fac, ...
    [p(1), p(2), neg_dist2.a, neg_dist2.b], p(3))));
    
[opt_V, V_fval, exit_flagV] = fmincon(obj_V, [5, 5, exptV], [], [], [], [], ...
    [sqrt(eps), sqrt(eps), 1 + sqrt(eps)], [], [], options_I);
vax_fit = @(x,t) time_dep(x, t/t_scale_fac, ...
    [opt_V(1:2), neg_dist2.a, neg_dist2.b], opt_V(end));

obj_VV = @(p) -sum(log(time_dep_two(S_vax_vax_VC, ...
    rel_days_vax_vax_VC/t_scale_fac, days_vax_vax_VC/t_scale_fac, ...
    [opt_V(1:2), p(1), p(2), neg_dist2.a, neg_dist2.b], opt_V(end), p(3))));
% ^ both exptI, exptIV vary as parameters

[opt_VV, VV_fval, exit_flagVV] = fmincon(obj_VV, [5, 5, exptVV], [], [], [], [], ...
    [sqrt(eps), sqrt(eps), 1 + sqrt(eps)], [], [], options_I);
t_delay_VV = 275;
VV_fit = @(x, t) time_dep_two_plot(x2, t2/t_scale_fac, t_delay_VV/t_scale_fac, ...
    [opt_V(1:2), opt_VV(1:2), neg_dist2.a, neg_dist2.b], opt_V(end), opt_VV(end));
VV_fit_2t = @(x, t, t_delay) time_dep_two_plot(x2, t2/t_scale_fac, t_delay/t_scale_fac, ...
    [opt_V(1:2), opt_VV(1:2), neg_dist2.a, neg_dist2.b], opt_V(end), opt_VV(end));

contI = [0.6 0.5 0.4 0.3 0.2 0.1 0.01, 0.001]; 

        V_mat = vax_fit(x2, t2);
        VV_mat = VV_fit(x2, t2);

prob_vals_I_LC = zeros(length(I_data_struct_LC), 7); % max # elements in a personal trajectory is 7
prob_vals_I_VC = zeros(length(I_data_struct_VC), 7); % max # elements in a personal trajectory is 7
prob_vals_V = zeros(length(V_data_struct_VC), 17); % max # elements in a personal trajectory is 17
prob_vals_IV_LC = zeros(length(IV_data_struct_LC), 7); % max # elements in a personal trajectory is 7
prob_vals_IV_VC = zeros(length(IV_data_struct_VC), 3); % max # elements in a personal trajectory is 3
prob_vals_VV = zeros(length(VV_data_struct_VC), 20); % max # elements in a personal trajectory is ?16


for k = 1:length(I_data_struct_LC)
    prob_vals_I_LC(k, 1:size(I_data_struct_LC{k},1)) ...
        = inf_fit(log2(I_data_struct_LC{k}.Spike_IgG+2)-1, I_data_struct_LC{k}.Days_Start_Inf)';
end
for k = 1:length(I_data_struct_VC)
    prob_vals_I_VC(k, 1:size(I_data_struct_VC{k},1)) ...
        = inf_fit(log2(I_data_struct_VC{k}.Spike_IgG+2)-1, I_data_struct_VC{k}.Days_Inf)';
end
for k = 1:length(V_data_struct_VC)
    prob_vals_V(k, 1:size(V_data_struct_VC{k},1)) ...
        = vax_fit(log2(V_data_struct_VC{k}.Spike_IgG+2)-1, V_data_struct_VC{k}.Days_Vax2_dose_2)';
end

low_prob_vals_I_LC = sum((prob_vals_I_LC > 0) & (prob_vals_I_LC < 0.001), 2);
low_prob_vals_I_VC = sum((prob_vals_I_VC > 0) & (prob_vals_I_VC < 0.001), 2);
high_prob_vals_I_LC = sum((prob_vals_I_LC > 0) & (prob_vals_I_LC > 0.2), 2);
high_prob_vals_I_VC = sum((prob_vals_I_VC > 0) & (prob_vals_I_VC > 0.2), 2);
high_prob_I_perc_LC = high_prob_vals_I_LC./sum(prob_vals_I_LC ~= 0, 2);
high_prob_I_perc_VC = high_prob_vals_I_VC./sum(prob_vals_I_VC ~= 0, 2);
low_prob_vals_V = sum((prob_vals_V > 0) & (prob_vals_V < 0.001), 2);
high_prob_vals_V = sum((prob_vals_V > 0) & (prob_vals_V > 0.2), 2);
high_prob_V_perc = high_prob_vals_V./sum(prob_vals_V ~= 0, 2);

%% Plotting infection, infection + vaccination models with personal trajectories

figure; 
hh = tiledlayout('flow');

nexttile 
histogram(S_neg, 10, 'Normalization','pdf', 'facecolor', [211,211,211]/256)
hold on, y = pdf(neg_dist, x); plot(x,y, 'k')
xlabel('Log-scaled antibody meas.', 'FontSize', 24)
ylabel('PDF', 'Fontsize', 24)
legend('Na{\"i}ve', 'fontsize', 24)
xlim([0, ceil(max(S_neg))]), ylim([0, 3])
hold off

nexttile 
hold on
contourf(t, x, I_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.25);
contour(t, x, I_mat,  contI, 'linewidth', 1.5, ...
    'linecolor', [0.2,0.2,0.2])
  splitapply(@myplot, days_inf_LC, S_inf_LC, G_inf_LC)
ax = gca; ax.ColorOrderIndex = 2;
h1 = scatter(days_inf, S_inf, 35, 'x', 'linewidth', 3, 'MarkerFaceColor','flat');
xlabel('Days post infection', 'fontsize', 24)
ylabel('Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')
legend(h1, 'Infection', 'fontsize', 24)
ylim([0, max(max(x))]) 
 colormap(ax, 'autumn'), clim([0, 0.7])
hold off

nexttile
  hold on
contourf(t, x, IV_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.2);
contour(t, x, IV_mat,  contI, 'linewidth', 1.5, ...
    'linecolor', [0.2,0.2,0.2])
  splitapply(@myplot, days_inf_LC, S_inf_LC, G_inf_LC)
splitapply(@myplot, data_inf_vax_VC.Days_Inf, ...
    log2(data_inf_vax_VC.Spike_IgG + 2) - 1, G_inf_vax_VC)
splitapply(@myplot, data_inf_vax_LC.Days_Start_Inf, ...
    log2(data_inf_vax_LC.Spike_IgG + 2) - 1, G_inf_vax_LC)
ax1 = gca; ax1.ColorOrderIndex = 2;
h1 = scatter([days_since_inf_LC(inf_label_LC); days_since_inf_VC(inf_label_VC)], ...
    [S_inf_LC; S_inf_VC], 35, 'x', ...
    'MarkerFaceColor','flat', 'linewidth', 3);
ax1.ColorOrderIndex = 5;
h2 = scatter([days_since_inf_LC(inf_vax_label_LC); ...
    days_since_inf_VC(inf_vax_label_VC)], [S_inf_vax_LC; ...
    S_inf_vax_VC], ...
    35, '^', 'MarkerFaceColor','flat', 'linewidth', 2.5);
xlabel('Days post infection', 'fontsize', 24)
legend( h2,'Infection then vaccination', 'fontsize', 24)
ylim([0,12]), xlim([0, 700]) 
colormap(ax1, 'summer'), clim([0, 0.7])
        
%% Plotting vaccination, vaccination^2 ONLY
% with personal trajectories

figure; 
hh = tiledlayout('flow');
nexttile 
hold on
contourf(t2, x2, V_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.2);
contour(t2, x2, V_mat,  contI, 'linewidth', 1.5, ...
    'linecolor', [0.2,0.2,0.2])
splitapply(@myplot, days_vax_VC, S_vax_VC, G_vax1_VC)
ax = gca; ax.ColorOrderIndex = 1;
h1 = scatter(days_vax_VC, S_vax_VC, 35, 'd', 'MarkerFaceColor','flat');
ylabel('Log-scaled antibody meas.', 'FontSize', 28, 'interpreter', 'latex')
legend(h1, 'Vaccination', 'fontsize', 24)
ylim([0,max(max(x2))]), xlim([0, 650])
 colormap(ax, 'parula'), clim([0, 0.7])
hold off

 nexttile
   hold on
 contourf(t2, x2, VV_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.15);
 contour(t2, x2, VV_mat,  contI, 'linewidth', 1.5, ...
     'linecolor', [0.2,0.2,0.2])
 splitapply(@myplot,days_vax_only_VC, S_vax_only_VC, G_vax_only_VC )
 ax1 = gca; ax1.ColorOrderIndex = 1;
 h1 = scatter(days_vax_VC, S_vax_VC, 35, 'd', 'MarkerFaceColor','flat');
ax1.ColorOrderIndex = 3;
 h2 = scatter(days_since_vax_VC(vax_vax_label_VC), ...
    S_VC(vax_vax_label_VC), 43, ...
     'MarkerFaceColor','flat', 'MarkerEdgeColor','none');
 legend(h2, 'Booster', 'fontsize', 24)
 ylim([0,max(max(x2))]), xlim([0, 650])
 colormap(ax1, 'hot'), clim([0, 0.7])
 hold off

 xlabel(hh, 'Days post vaccination', 'fontsize', 28, 'interpreter', 'latex')

%% Look at a few personal trajectories of interest of infection
% Representative trajectories and outliers

I_interest = [8, 22, 28, 37, 42];
leg_I = cell(length(I_interest), 1);
clear hi
figure; hold on
contourf(t, x, I_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.25);
    contour(t, x, I_mat, contI, 'linewidth', 1.5, ...
        'linecolor', [0.596, 0.486, 0.463], 'ShowText','on')
    colormap(gca, 'summer'), clim([0, 0.7])
for k = 1:length(I_interest)
    data_k = I_data_struct_LC{I_interest(k)};
    days_inf_k = data_k.Days_Start_Inf;
    S_k = log2(data_k.Spike_IgG + 2) - 1;
    ax = gca;
    plot(days_inf_k, S_k, 'k')
    hi(k) =  scatter(days_inf_k, S_k, 50, 'x', 'MarkerFaceColor','flat', ...
        'linewidth', 4);
    leg_I = {'Subject D1-1', 'Subject D1-2', 'Subject D1-3', 'Subject D1-4', 'Subject D1-5'};
    yticks([0, 4, 8, 12, 16])
    xlabel('Days post infection', 'FontSize', 24, 'interpreter', 'latex')
    ylabel('Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')
    xlim([0,500]), ylim([0, 12])
end
legend(hi, leg_I, 'fontsize', 20, 'location', 'northeast')

%% Look at a few personal trajectories of interest of infection + vax (LC)
% Representative trajectories and outliers

I_interest = [27,  32, 57];
clear hi
figure; 
hh = tiledlayout('flow'); 
for k = 1:length(I_interest)
    nexttile, hold on
    data_k = data_struct_LC{I_interest(k)};
    days_inf_k = data_k.Days_Start_Inf;
    days_I_k = data_k.Days_Start_Inf(strcmp(data_k.Vaxed, 'n'));
    days_IV_k = days_inf_k(strcmp(data_k.Vaxed, 'y'));
    S_k = log2(data_k.Spike_IgG + 2) - 1;
    S_I_k = S_k(strcmp(data_k.Vaxed, 'n'));
    S_IV_k = S_k(strcmp(data_k.Vaxed, 'y'));
    rel_days_k = rmmissing(data_k.Days_Start_Inf - data_k.Days_Vax_1);
    rel_days_k = rel_days_k(1);
    IV_fit_pers = @(x, t) time_dep_two_plot(x, t/t_scale_fac, rel_days_k/t_scale_fac, ...
    [opt_I(1:2), opt_IV(1:2), neg_dist.a, neg_dist.b], opt_I(end), opt_IV(end));
    IV_mat_pers = IV_fit_pers(x, t);
    contourf(t, x, IV_mat_pers, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.25);
    contour(t, x, IV_mat_pers, contI, 'linewidth', 1.5, ...
        'linecolor', [0.2,0.2,0.2], 'ShowText','on')
    colormap(gca, 'autumn'), clim([0, 0.7])
    plot(ones(100,1)*rel_days_k, linspace(0, 12, 100), 'w--')
    text(rel_days_k-20, 8, 'Vaccination', 'color', 'w', 'fontsize', 20, 'Rotation', 90)
    ax = gca; ax.ColorOrderIndex = k+3;
    plot(days_inf_k, S_k, 'k')
    hi(k) =  scatter(days_I_k, S_I_k, 50, 'x', 'MarkerFaceColor','flat', ...
        'linewidth', 4);
    ax.ColorOrderIndex = k+3;
    hy(k) = scatter(days_IV_k, S_IV_k, 50, '^', 'MarkerFaceColor','flat', 'linewidth', 3);
    yticks([0, 4, 8, 12, 16])
    if k == 1
        ylabel('Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')
    end
    xlim([0,700]), ylim([0, 12])
end
leg_IV = {'Subject D1-12', 'Subject D1-13', 'Subject D1-14'};
    
legend(hi, leg_IV, 'fontsize', 20, 'location', 'northeast')
 xlabel(hh, 'Days post infection', 'FontSize', 24, 'interpreter', 'latex')

%% Look at a few personal trajectories of interest of vaccination
% Representative trajectories and outliers

V_interest = [5, 7, 9, 111, 135, 149];
leg_V = cell(length(V_interest), 1);
clear hi
figure; hold on
contourf(t, x, V_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.25);
    contour(t, x, V_mat, contI, 'linewidth', 1.5, ...
        'linecolor', [0.2,0.2,0.2], 'ShowText','on')
    colormap(gca, 'parula'), clim([0, 0.3])
for k = 1:length(V_interest)
    data_k = V_data_struct_VC{V_interest(k)};
    data_k = data_k(logical(strcmp(data_k.Inf_Vax_traj, 'Vax') ...
    & ~isnan(data_k.Spike_IgG) & ~isnan(data_k.Days_Vax2_dose_2)),:);
    days_vax_k = data_k.Days_Vax2_dose_2;
    S_k = log2(data_k.Spike_IgG + 2) - 1;
    ax = gca;  
    plot(days_vax_k, S_k, 'k')
    hi(k) =  scatter(days_vax_k, S_k, 50, 'd', 'MarkerFaceColor','flat');
    leg_V = {'Subject D1-6', 'Subject D1-7', 'Subject D1-8', 'Subject D1-9', 'Subject D1-10', 'Subject D1-11'};
    yticks([0, 4, 8, 12, 16])
    xlabel('Days post vaccination', 'FontSize', 24, 'interpreter', 'latex')
    ylabel('Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')
    xlim([0,500]), ylim([0, 12])
end
legend(hi, leg_V, 'fontsize', 20, 'location', 'northeast')

%% Look at a few personal trajectories of interest of vaccination + boost
% Representative trajectories and outliers

V_interest = [4, 6, 16, 38, 63, 73];
leg_V = cell(length(V_interest), 1);
clear hi
figure; 
hh = tiledlayout('flow'); 
for k = 1:length(V_interest)
    nexttile, hold on
    data_k = VV_data_struct_VC{V_interest(k)};
    days_vax_k = data_k.Days_Vax2_dose_2;
    S_k = log2(data_k.Spike_IgG + 2) - 1;
    S_vax2_k = S_k(strcmp(data_k.Inf_Vax_traj, 'Vax'));
    S_boost_k = S_k(strcmp(data_k.Inf_Vax_traj, 'Vax_Vax'));
    days_vax2_k = data_k.Days_Vax2_dose_2(strcmp(data_k.Inf_Vax_traj, 'Vax'));
    days_boost_k = data_k.Days_Vax2_dose_2(strcmp(data_k.Inf_Vax_traj, 'Vax_Vax'));
    rel_days_k = rmmissing(data_k.Days_Vax2_dose_2 - data_k.Days_Vax3_booster);
    rel_days_k = rel_days_k(1);
    VV_fit_pers = @(x, t) time_dep_two_plot(x, t/t_scale_fac, rel_days_k/t_scale_fac, ...
    [opt_V(1:2), opt_VV(1:2), neg_dist.a, neg_dist.b], opt_V(end), opt_VV(end));
    VV_mat_pers = VV_fit_pers(x, t);
    contourf(t, x, VV_mat_pers, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.15);
    contour(t, x, VV_mat_pers, contI, 'linewidth', 1.5, ...
        'linecolor', [0.2,0.2,0.2], 'ShowText','on')
    colormap(gca, 'hot'), clim([0, 0.7])
    plot(ones(100,1)*rel_days_k, linspace(0, 12, 100), 'w--')
    text(rel_days_k-20, 8, 'Booster', 'color', 'w', 'fontsize', 20, 'Rotation', 90)
    ax = gca;  ax.ColorOrderIndex = k+3;
    plot(days_vax_k, S_k, 'k')
    hi(k) =  scatter(days_vax2_k, S_vax2_k, 50, 'd', 'MarkerFaceColor','flat');
     ax.ColorOrderIndex = k+3;
    hy(k) =  scatter(days_boost_k, S_boost_k, 50, 'o', 'MarkerFaceColor','flat');
    yticks([0, 4, 8, 12, 16])
    xlim([0,700]), ylim([0, 12])
end
leg_VV = {'Subject D1-15', 'Subject D1-16', 'Subject D1-17', 'Subject D1-18', 'Subject D1-19', 'Subject D1-20'};
    
ylabel(hh, 'Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')
 xlabel(hh, 'Days post vaccination', 'FontSize', 24, 'interpreter', 'latex')
 legend(hi, leg_VV, 'fontsize', 20, 'location', 'northeast')

 %% Read in and pre-process data
% LC = long COVID cohort
% VC = vaccine cohort

data_LC_all = readtable('../../Datasets/UVA/UVA_antibody_data.xlsx', ...
    'sheet', 'Timedep_meas_ONLY');
data_VC = readtable('../../Datasets/UVA/UVA_antibody_data_vax_cohort_days_since.xlsx');

data_LC = data_LC_all(~isnan(data_LC_all.Days_Start_Inf), :);

inf_label_LC = logical(isnan(data_LC.Days_LastVax) ...
    & ~isnan(data_LC.Spike_IgG) & ~isnan(data_LC.Days_Start_Inf));
inf_vax_label_LC = logical(~isnan(data_LC.Days_Vax_1) &  ...
    isnan(data_LC.Days_Vax_2) & ~strcmp(data_LC.Days_Start_Inf, 'NA') ...
    & ~isnan(data_LC.Spike_IgG) & (data_LC.Days_Start_Inf > data_LC.Days_LastVax));
inf_vax_vax_label_LC = logical(~isnan(data_LC.Days_Vax_2) & ...
     isnan(data_LC.Days_Vax_3) & ~strcmp(data_LC.Days_Start_Inf, 'NA') ...
     & ~isnan(data_LC.Spike_IgG));

inf_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Inf') ...
    & ~isnan(data_VC.Spike_IgG) & data_VC.Days_Inf > 0); 
inf_vax_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Inf_Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax1_dose_1));
inf_vax_vax_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Inf_Vax_Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax3_booster));
inf_first_label_VC = inf_label_VC | inf_vax_label_VC | inf_vax_vax_label_VC;

vax_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax1_dose_1));
% ^ dose 1 and dose 2
vax1_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax') ...
    & ~isnan(data_VC.Spike_IgG) & isnan(data_VC.Days_Vax2_dose_2));
vax2_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax2_dose_2));
vax_inf_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax_Inf') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Inf));
vax_vax_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax_Vax') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Vax3_booster));
vax_vax_inf_label_VC = logical(strcmp(data_VC.Inf_Vax_traj, 'Vax_Vax_Inf') ...
    & ~isnan(data_VC.Spike_IgG) & ~isnan(data_VC.Days_Inf));
vax_only_label_VC = vax_label_VC | vax_vax_label_VC;
vax_first_label_VC = vax_label_VC | vax_vax_label_VC | vax_inf_label_VC ...
    | vax_vax_inf_label_VC;

days_since_inf_LC = data_LC.Days_Start_Inf;

days_inf_LC = data_LC.Days_Start_Inf(inf_label_LC);
days_inf_vax_LC = data_LC.Days_Vax_1(inf_vax_label_LC); 
days_inf_vax_vax_LC = data_LC.Days_Vax_2(inf_vax_vax_label_LC); 

days_since_VC = data_VC.Days_since_first_event;
days_since_VC(vax_first_label_VC) = data_VC.Days_Vax1_dose_1(vax_first_label_VC);
% ^ after dose 1

days_since_inf_VC = data_VC.Days_Inf;
days_inf_VC = data_VC.Days_Inf(inf_label_VC);
days_inf_vax_VC = data_VC.Days_Vax1_dose_1(inf_vax_label_VC);
% ^ after dose 1
days_inf_vax_vax_VC = data_VC.Days_Vax3_booster(inf_vax_vax_label_VC);

days_since_vax_VC = data_VC.Days_Vax1_dose_1; 
days_vax_VC = data_VC.Days_Vax1_dose_1(vax_label_VC);
% ^ after dose 1
days_vax_vax_VC = data_VC.Days_Vax3_booster(vax_vax_label_VC);
days_vax_vax_inf_VC = data_VC.Days_Inf(vax_vax_inf_label_VC);

rel_days_inf_vax_LC = data_LC.Days_Start_Inf(inf_vax_label_LC) - days_inf_vax_LC;
rel_days_inf_vax_LC_IVV = data_LC.Days_Start_Inf(inf_vax_vax_label_LC)...
    - data_LC.Days_Start_Vax(inf_vax_vax_label_LC);
rel_days_inf_vax_vax_LC = data_LC.Days_Start_Vax(inf_vax_vax_label_LC) ...
    - days_inf_vax_vax_LC;

rel_days_vax_vax_VC = rmmissing(data_VC.Days_Vax1_dose_1(vax_vax_label_VC))...
    - days_vax_vax_VC;
% ^ after dose 1
ng_VC = find(rel_days_vax_vax_VC < 0 | isnan(rel_days_vax_vax_VC));
rel_days_vax_vax_VC(ng_VC) = [];
days_vax_vax_VC(ng_VC) = [];
rel_days_vax_vax_VC_VVI = data_VC.Days_Vax1_dose_1(vax_vax_inf_label_VC)...
    - data_VC.Days_Vax3_booster(vax_vax_inf_label_VC);
% ^ after dose 1
rel_days_vax_vax_inf_VC = data_VC.Days_Vax3_booster(vax_vax_inf_label_VC) ...
    - days_vax_vax_inf_VC;
rel_days_inf_vax_VC = data_VC.Days_since_first_event(inf_vax_label_VC) ...
    - days_inf_vax_VC;
rel_days_inf_vax_VC_IVV = data_VC.Days_since_first_event(inf_vax_vax_label_VC) ...
    - data_VC.Days_Vax1_dose_1(inf_vax_vax_label_VC);
ng_VC2 = find(rel_days_inf_vax_VC < 0 | isnan(rel_days_inf_vax_VC));
rel_days_inf_vax_VC(ng_VC2) = [];
days_inf_vax_VC(ng_VC2) = [];
rel_days_inf_vax_vax_VC = data_VC.Days_Vax1_dose_1(inf_vax_vax_label_VC) ...
    - days_inf_vax_vax_VC;

days_vax_only_VC = days_since_vax_VC(vax_only_label_VC);

S_LC = log2(data_LC.Spike_IgG + 2) - 1; 
% ^ log-scale and shift to avoid negative values
%S_LC(isnan(S_LC)) = [];
S_neg_LC = data_LC_all.Spike_IgG(strcmp(data_LC_all.Visit_ID, 'Base'));
S_inf_LC = S_LC(inf_label_LC);
S_inf_vax_LC = S_LC(inf_vax_label_LC);
S_inf_vax_vax_LC = S_LC(inf_vax_vax_label_LC); 

S_VC = log2(data_VC.Spike_IgG + 2) - 1; 
S_nans = find(isnan(S_VC));
S_vax_only_VC = S_VC(vax_only_label_VC);

S_inf_VC = S_VC(inf_label_VC);
S_vax_VC = S_VC(vax_label_VC); 
S_vax_vax_VC = S_VC(vax_vax_label_VC);
S_vax_vax_VC(ng_VC) = [];
S_inf_vax_VC = S_VC(inf_vax_label_VC);
S_inf_vax_VC(ng_VC2) = [];
S_inf_vax_vax_VC = S_VC(inf_vax_vax_label_VC);
S_vax_vax_inf_VC = S_VC(vax_vax_inf_label_VC);

S_neg_VC = data_VC.Spike_IgG(data_VC.Days_Vax1_dose_1 <= 0 ...
    & isnan(data_VC.Days_Inf));
S_neg = [S_neg_LC; S_neg_VC];
S_neg = log2(S_neg + 2) - 1;

G_LC = findgroups(data_LC.Subject_ID);
data_vax_LC = data_LC(~inf_label_LC, :);
G_vax_LC = findgroups(data_vax_LC.Subject_ID);
data_inf_LC = data_LC(inf_label_LC, :);
G_inf_LC = findgroups(data_inf_LC.Subject_ID);

G_VC = findgroups(data_VC.Subject_ID);
data_inf_VC = data_VC(inf_first_label_VC, :);
data_vax_VC = data_VC(vax_first_label_VC, :);
data_vax1_VC = data_VC(vax_label_VC, :);
data_vax_only_VC = data_VC(vax_only_label_VC, :);
G_vax_VC = findgroups(data_vax_VC.Subject_ID);
G_vax1_VC = findgroups(data_vax1_VC.Subject_ID);
G_vax_only_VC = findgroups(data_vax_only_VC.Subject_ID);
G_inf_VC = findgroups(data_inf_VC.Subject_ID);

data_inf_VC_days = data_inf_VC.Days_since_first_event;
data_inf_VC_S = log2(data_inf_VC.Spike_IgG+2)-1;

data_vax_VC_days = data_vax_VC.Days_Vax1_dose_1;
data_vax_VC_S = log2(data_vax_VC.Spike_IgG+2)-1;

I_data_struct_LC = struct([]); 
for i = 1:G_inf_LC(end)
    I_data_struct_LC{i} = data_inf_LC(G_inf_LC == i, :);
end

IV_data_struct_LC = struct([]);
data_inf_vax_LC = [];
k = 1;
G_inf_vax_LC = []; % group indices for crossover events
for i = 1:G_LC(end)
    l = find(G_LC == i);
    data_l = data_LC(l,:);
    cr_i = sum(logical(strcmp(data_l.Vaxed, 'n')));
    cr_iv = sum(logical(strcmp(data_l.Vaxed, 'y')));
    cr_reld = diff(data_l.Days_Start_Inf);
    if (cr_i > 0 && cr_iv > 0 && all(cr_reld > 0))
        data_l = data_l((strcmp(data_l.Vaxed, 'n') ...
            | strcmp(data_l.Vaxed, 'y')) & data_l.Days_Start_Inf > 0 ...
           & ~isnan(data_l.Spike_IgG) & isnan(data_l.Days_Vax_2), :);
        data_inf_vax_LC = [data_inf_vax_LC; data_l];
        IV_data_struct_LC{k} = data_l;
        k = k + 1;
    end
end

 G_inf_vax_LC = findgroups(data_inf_vax_LC.Subject_ID);

IV_data_struct_VC = struct([]);
data_inf_vax_VC = [];
G_inf_vax_VC = []; % group indices for crossover events
k = 1;
for i = 1:G_VC(end)
    l = find(G_VC == i);
    data_l = data_VC(l,:);
    cr_i = sum(logical(strcmp(data_l.Inf_Vax_traj, 'Inf')));
    cr_iv = sum(logical(strcmp(data_l.Inf_Vax_traj, 'Inf_Vax')));
    cr_reld = diff(data_l.Days_since_first_event);
    if ((cr_i > 0 || cr_iv > 0))
        data_l = data_l(((strcmp(data_l.Inf_Vax_traj, 'Inf') ...
            | strcmp(data_l.Inf_Vax_traj, 'Inf_Vax')) & ...
            ~isnan(data_l.Spike_IgG) & data_l.Days_Inf > 0 ...
            & ~isnan(data_l.Days_Vax1_dose_1)), :);
        if ~isempty(data_l)
            data_inf_vax_VC = [data_inf_vax_VC; data_l];
            IV_data_struct_VC{k} = data_l;
            k = k + 1;
        end
    end
end

G_inf_vax_VC = findgroups(data_inf_vax_VC.Subject_ID);

V_data_struct_VC = struct([]); 
for i = 1:G_vax_VC(end)
    V_data_struct_VC{i} = data_vax_VC(G_vax_VC == i, :);
end

VVI_data_struct_VC = struct([]);
k = 1;
data_vax_vax_inf = [];
for i = 1:G_vax_VC(end)
    temp = data_vax_VC(G_vax_VC == i, :);
    ct = length(find(strcmp(temp.Inf_Vax_traj, 'Vax_Vax_Inf')));
    if ct > 0
        data_vax_vax_inf = [data_vax_vax_inf; temp];
        VVI_data_struct_VC{k} = temp;
        k = k + 1;
    end
end

G_vax_inf = findgroups(data_vax_vax_inf.Subject_ID);

V_only_data_struct_VC = struct([]); 
for i = 1:G_vax_only_VC(end)
    V_only_data_struct_VC{i} = data_vax_only_VC(G_vax_only_VC == i, :);
end

I_data_struct_VC = struct([]); 
for i = 1:G_inf_VC(end)
    I_data_struct_VC{i} = data_inf_VC(G_inf_VC == i, :);
end

 %% Plotting vaccination, vaccination^2, and vaccination, + infection events
% with personal trajectories

figure; 
hh = tiledlayout('flow');
nexttile 
hold on
contourf(t2, x2, V_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.2);
contour(t2, x2, V_mat,  contI, 'linewidth', 1.5, ...
    'linecolor', [0.2,0.2,0.2])
splitapply(@myplot, days_vax_VC, S_vax_VC, G_vax1_VC)
ax = gca; ax.ColorOrderIndex = 1;
h1 = scatter(days_since_vax_VC(vax1_label_VC), S_VC(vax1_label_VC), ...
    35, 'd',  'MarkerFaceColor','flat', 'linewidth', 2);
ax.ColorOrderIndex = 4;
h2 = scatter(days_since_vax_VC(vax2_label_VC), S_VC(vax2_label_VC), ...
    35, 'd',  'MarkerFaceColor','none', 'linewidth', 2);
ylabel('Log-scaled antibody meas.', 'FontSize', 28, 'interpreter', 'latex')
legend([h1, h2], 'Vax Dose 1', 'Vax Dose 2', 'fontsize', 24)
ylim([0,max(max(x2))]), xlim([0, 650])
 colormap(ax, 'parula'), clim([0, 0.3])
hold off

 nexttile
   hold on
 contourf(t2, x2, VV_mat, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.15);
 contour(t2, x2, VV_mat,  contI, 'linewidth', 1.5, ...
     'linecolor', [0.2,0.2,0.2])
 splitapply(@myplot,days_vax_only_VC, S_vax_only_VC, G_vax_only_VC )
 ax1 = gca; ax1.ColorOrderIndex = 1;
h1 = scatter(days_since_vax_VC(vax1_label_VC), S_VC(vax1_label_VC), ...
    35, 'd',  'MarkerFaceColor','flat', 'linewidth', 2);
ax1.ColorOrderIndex = 4;
h3 = scatter(days_since_vax_VC(vax2_label_VC), S_VC(vax2_label_VC), ...
    35, 'd',  'MarkerFaceColor','none', 'linewidth', 2);
ax1.ColorOrderIndex = 3;
 h2 = scatter(days_since_vax_VC(vax_vax_label_VC), ...
    S_VC(vax_vax_label_VC), 43, ...
     'MarkerFaceColor','flat', 'MarkerEdgeColor','none');
 legend(h2, 'Booster', 'fontsize', 24)
 ylim([0,max(max(x2))]), xlim([0, 650])
 colormap(ax1, 'hot'), clim([0, 0.3])
 hold off

 xlabel(hh, 'Days since vaccination', 'fontsize', 28, 'interpreter', 'latex')

%% Movie for booster with personal trajectories one at a time

% first, find and sort relative days between vaccination and booster
rel_days_to_sort = zeros(length(VV_data_struct_VC), 1);
for k = 1:length(VV_data_struct_VC)
    data_k = VV_data_struct_VC{k};
    days_vax_k = data_k.Days_Vax2_dose_2;
    days_vax2_k = data_k.Days_Vax2_dose_2(strcmp(data_k.Inf_Vax_traj, 'Vax'));
    days_boost_k = data_k.Days_Vax2_dose_2(strcmp(data_k.Inf_Vax_traj, 'Vax_Vax'));
    rel_days_k = rmmissing(data_k.Days_Vax2_dose_2 - data_k.Days_Vax3_booster);
    rel_days_to_sort(k) = rel_days_k(1);
end
[rel_days_sorted, rel_days_sort_idx] = sort(rel_days_to_sort);

file_name = 'UVA_two_event_model_25_vax_boost_individ.gif';
for k = 1:length(VV_data_struct_VC)
    zoom off
    drawnow

     figure; hold on
    data_k = VV_data_struct_VC{rel_days_sort_idx(k)};
    days_vax_k = data_k.Days_Vax2_dose_2;
    S_k = log2(data_k.Spike_IgG + 2) - 1;
    S_vax2_k = S_k(strcmp(data_k.Inf_Vax_traj, 'Vax'));
    S_boost_k = S_k(strcmp(data_k.Inf_Vax_traj, 'Vax_Vax'));
    days_vax2_k = data_k.Days_Vax2_dose_2(strcmp(data_k.Inf_Vax_traj, 'Vax'));
    days_boost_k = data_k.Days_Vax2_dose_2(strcmp(data_k.Inf_Vax_traj, 'Vax_Vax'));
    rel_days_k = rmmissing(data_k.Days_Vax2_dose_2 - data_k.Days_Vax3_booster);
    rel_days_k = rel_days_k(1);
    VV_fit_pers = @(x, t) time_dep_two_plot(x, t/t_scale_fac, rel_days_k/t_scale_fac, ...
    [opt_V(1:2), opt_VV(1:2), neg_dist.a, neg_dist.b], opt_V(end), opt_VV(end));
    VV_mat_pers = VV_fit_pers(x, t);
    contourf(t, x, VV_mat_pers, 100, 'EdgeColor', 'None', 'FaceAlpha', 0.15);
    contour(t, x, VV_mat_pers, contI, 'linewidth', 1.5, 'linecolor', [0.2,0.2,0.2])
    colormap(gca, 'hot'), clim([0, 0.7])
    plot(ones(100,1)*rel_days_k, linspace(0, 12, 100), 'w--')
    text(rel_days_k-20, 8, 'Booster', 'color', 'w', 'fontsize', 20, 'Rotation', 90)
    ax = gca;  ax.ColorOrderIndex = k+3;
    plot(days_vax_k, S_k, 'k')
    hi(k) =  scatter(days_vax2_k, S_vax2_k, 50, 'd', 'MarkerFaceColor','flat');
     ax.ColorOrderIndex = k+3;
    hy(k) =  scatter(days_boost_k, S_boost_k, 50, 'o', 'MarkerFaceColor','flat');
    leg_V{k} = ['Subject ', num2str(data_k.Subject_ID(1))];
    yticks([0, 4, 8, 12, 16])
        ylabel('Log-scaled antibody meas.', 'FontSize', 24, 'interpreter', 'latex')
    xlim([0,700]), ylim([0, 12])
     xlabel('Days post vaccination', 'FontSize', 24, 'interpreter', 'latex')

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

function W = time_dep_two(r, t_btwn_events, t_event2_meas, p, expt1, expt2)
% two-event distribution

shape_1 = shapetpos(t_btwn_events, [p(1), p(2), p(5)], expt1);
shape_2 = shapetpos(t_event2_meas, [p(3), p(4), p(5)], expt2);
shape = shape_1 + shape_2 - p(5);
scale = p(6);
W = gampdf(r, shape, scale);

end

function W = time_dep_two_plot(r, t, t_delay, p, expt1, expt2)
% two-event distribution
%code by Joseph Nakao!
shape_1 = shapetpos(t, [p(1), p(2), p(5)], expt1);

t_since = t - t_delay;
shape_2 = shapetpos(t_since, [p(3), p(4), p(5)], expt2);
shape_2(t < t_delay) = 0;
shape = shape_1 + shape_2;
shape(t >= t_delay) = shape(t >= t_delay) - p(5);

scale = p(6);
W = gampdf(r, shape, scale);

end
