% load GAMM fitting resutls,
ft_perc_data = table2array(readtable(['../01_growth_curve/saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']));
pt_perc_data = table2array(readtable(['../01_growth_curve/saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']));
% run it.
%inter_code_2
ft_pma = ft_perc_data(:,2);
ft_pert = ft_perc_data(:,14:18);

pt_pma = pt_perc_data(:,2);
pt_pert = pt_perc_data(:,14:18);

fill([pt_pma; flipud(pt_pma)], ...
    [pt_pert(:,1); flipud(pt_pert(:,5))], ...
    c_deeppink, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); hold on;

% Shaded area for preterm: 25–75%
fill([pt_pma; flipud(pt_pma)], ...
    [pt_pert(:,2); flipud(pt_pert(:,4))], ...
    c_deeppink, 'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on;

% Shaded area for term: 5–95%
fill([ft_pma; flipud(ft_pma)], ...
    [ft_pert(:,1); flipud(ft_pert(:,5))], ...
    c_deepblue, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); hold on;

% Shaded area for term: 25–75%
fill([ft_pma; flipud(ft_pma)], ...
    [ft_pert(:,2); flipud(ft_pert(:,4))], ...
    c_deepblue, 'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on;

for ind = 1:5
    plot(ft_pma,ft_pert(:,ind),'color',c_deepblue,'linewidth',0.8); hold on;
end
plot(ft_pma,ft_pert(:,3),'color',c_deepblue,'linewidth',1.5); hold on;

for ind = 1:5
    plot(pt_pma,pt_pert(:,ind),'color',c_deeppink,'linewidth',0.8); hold on;
end
plot(pt_pma,pt_pert(:,3),'color',c_deeppink,'linewidth',1.5); hold on;
grid
max_value = max(ft_pert(:));
ax = gca; ax.YAxis.Exponent = floor(log10(abs(max_value)));
set(findall(gcf,'-property','fontweight'),'fontweight','bold');
title(cur_region_name);
