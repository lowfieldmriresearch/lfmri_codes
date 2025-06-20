% Extract GAMM infos:
PMA = conf_data.PMA;
Group = conf_data.Group;
brain_volume = conf_data.brain_volume;
lower_5 = conf_data.lower_5;
upper_95 = conf_data.upper_95;
lower = conf_data.lower;
upper = conf_data.upper;
raw_x = conf_data.raw_x;
raw_y = conf_data.raw_y;

% extract percentile infos:
ft_p_gpinfo = ft_perc_data.Group;
ft_p_PMA = ft_perc_data.PMA;
ft_p_5 = ft_perc_data.smoothed_p5;
ft_p_25 = ft_perc_data.smoothed_p25;
ft_p_50 = ft_perc_data.smoothed_p50;
ft_p_75 = ft_perc_data.smoothed_p75;
ft_p_95 = ft_perc_data.smoothed_p95;


pt_p_gpinfo = pt_perc_data.Group;
pt_p_PMA = pt_perc_data.PMA;
pt_p_5 = pt_perc_data.smoothed_p5;
pt_p_25 = pt_perc_data.smoothed_p25;
pt_p_50 = pt_perc_data.smoothed_p50;
pt_p_75 = pt_perc_data.smoothed_p75;
pt_p_95 = pt_perc_data.smoothed_p95;

% Create a figure
% figure;

% Separate data by Group (preterm and full-term)
group1 = Group == 0; %'preterm';  % Assuming 'Group' is categorical
group2 = Group == 1;%'full-term';

PMA_gp1 = PMA(group1);
upper_gp1 = upper(group1);
lower_gp1 = lower(group1);
brain_volume_gp1 = brain_volume(group1);
upper_95_gp1 = upper_95(group1);
lower_95_gp1 = lower_5(group1);

[~,gp1_order]= sort(PMA_gp1);

PMA_gp2 = PMA(group2);
upper_gp2 = upper(group2);
lower_gp2 = lower(group2);
brain_volume_gp2 = brain_volume(group2);
upper_95_gp2 = upper_95(group2);
lower_95_gp2 = lower_5(group2);

[~,gp2_order]= sort(PMA_gp2);

% scatters.
scatter(PMA(group2),raw_y(group2),'filled','markerfacecolor',c_deepblue,'markerfacealpha',0.25); hold on;
scatter(PMA(group1),raw_y(group1),'filled','markerfacecolor',c_deeppink,'markerfacealpha',0.25); hold on;



% Plot the 95% confidence intervals as shaded regions
fill([PMA_gp1(gp1_order); flipud(PMA_gp1(gp1_order))], ...
    [upper_gp1(gp1_order); flipud(lower_gp1(gp1_order))], c_lightpink, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
fill([PMA_gp2(gp2_order); flipud(PMA_gp2(gp2_order))], ...
    [upper_gp2(gp2_order); flipud(lower_gp2(gp2_order))], c_lightblue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

% Plot the 5% and 95% intervals (90% CI)
fill([PMA_gp1(gp1_order); flipud(PMA_gp1(gp1_order))], ...
    [upper_95_gp1(gp1_order); flipud(lower_95_gp1(gp1_order))], c_lightpink, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
fill([PMA_gp2(gp2_order); flipud(PMA_gp2(gp2_order))], ...
    [upper_95_gp2(gp2_order); flipud(lower_95_gp2(gp2_order))], c_lightblue, 'FaceAlpha', 0.8, 'EdgeColor', 'none');

% Plot the predicted brain volume for Preterm group
plot(PMA_gp1(gp1_order), brain_volume_gp1(gp1_order), 'color',c_deeppink*0.75, 'LineWidth', 1.5);
hold on;

% Plot the predicted brain volume for Full-term group
plot(PMA_gp2(gp2_order), brain_volume_gp2(gp2_order), 'color',c_deepblue*0.75, 'LineWidth', 1.5);

% plot all the percentiles.
plot(ft_p_PMA, ft_p_5,'color',c_deepblue,'linestyle','-.','linewidth',1.2); hold on;
%plot(ft_p_PMA, ft_p_25,'color','b','linestyle','--','linewidth',1.5); hold on;
%plot(ft_p_PMA, ft_p_75,'color','b','linestyle','--','linewidth',1.5); hold on;
plot(ft_p_PMA, ft_p_95,'color',c_deepblue,'linestyle','-.','linewidth',1.2); hold on;

plot(pt_p_PMA, pt_p_5,'color',c_deeppink,'linestyle','-.','linewidth',1.2); hold on;
%plot(pt_p_PMA, pt_p_25,'color','r','linestyle','--','linewidth',1.5); hold on;
%plot(pt_p_PMA, pt_p_75,'color','r','linestyle','--','linewidth',1.5); hold on;
plot(pt_p_PMA, pt_p_95,'color',c_deeppink,'linestyle','-.','linewidth',1.2); hold on;

% Add labels and title
%xlabel('PMA');
%ylabel('Brain Volume');

max_value = max(pt_p_95);
% set proper rannge for exponent,
ax = gca; ax.YAxis.Exponent = floor(log10(abs(max_value)));

% Customize the plot
grid on;
hold off;
xlim([30,52])

% add significance info,
cur_p  = p_res{cur_V_ind-1,4};

if cur_p<0.001
    p_txt = 'p<0.001';
elseif cur_p<0.005
    p_txt = 'p<0.005';
elseif cur_p<0.01
    p_txt = 'p<0.01';
elseif cur_p<0.05
    p_txt = 'p<0.05';
else
    p_txt = '';
end

ylow = min(ylim);
yhigh = max(ylim);
set(findall(gca,'-property','fontweight'),'fontweight','bold');
set(findall(gca,'-property','fontsize'),'fontsize',11);
title(cur_region_name,'fontsize',13,'fontweight','bold');
text(32,(yhigh-ylow)*0.8 + ylow,p_txt,'fontsize',12,'fontweight','bold');