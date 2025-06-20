clear; clc; close all;

% COLOR PROFILE,
% color profiling:
c_deepblue = [87,111,160]/256;
c_lightblue = [167,185,215]/256;
c_deeppink = [181,121,121]/256;
c_lightpink = [222,163,162]/256;

tt = table2array(readtable('../raw_data/full_res_for_gamm00_use.csv'));

ptag = tt(:,1);
GA = tt(:,85);
PMA = tt(:,86);

% correct GA.
GA = (GA-floor(GA))/7*10 + floor(GA);

scatter(GA(ptag==0),PMA(ptag==0),56,'filled','markerfacecolor',c_deeppink,'markerfacealpha',0.25); hold on;
scatter(GA(ptag==1),PMA(ptag==1),56,'filled','markerfacecolor',c_deepblue,'markerfacealpha',0.25); 
xlabel('GA at birth (weeks)')
ylabel('PMA at MRI scan (weeks)')

set(findall(gcf,'-property','fontweight'),'fontweight','bold');
set(findall(gcf,'-property','fontsize'),'fontsize',14);

grid off
xlim([24,43])
ylim([30,54])