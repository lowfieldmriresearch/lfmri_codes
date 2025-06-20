clear;clc;close all;
fullres = readmatrix('../mat_for_nutr_analysis_zscore.csv');
full_lab_res = fullres(:,49:86);
mri_res = fullres(:,2:48);
mri_res(:,48) = mean(mri_res(:,14:47),2);

% need to polish the style of boxplot.
% color profiling:
c_deepblue = [87,111,160]/256;
c_lightblue = [167,185,215]/256;
c_deeppink = [181,121,121]/256;
c_lightpink = [222,163,162]/256;

%focus on albumin
close all;
mri_ind = 48;
lab_ind = 28;
lab_res = full_lab_res(:,lab_ind);
lab_nonnan_ind = (~isnan(lab_res).*lab_res>0);
pt_subject_index = find((fullres(:,1)==0).*(lab_nonnan_ind));
ft_subject_index = find((fullres(:,1)==1).*(lab_nonnan_ind));



%
figure,

% layer #1 , scatter plot
scatter(lab_res(ft_subject_index),mri_res(ft_subject_index,mri_ind),60,'filled','markerfacecolor',c_deepblue,'markerfacealpha',0.15,'markeredgecolor','none'); hold on;
scatter(lab_res(pt_subject_index),mri_res(pt_subject_index,mri_ind),60,'filled','markerfacecolor',c_deeppink,'markerfacealpha',0.15,'markeredgecolor','none'); hold on;
% 

% layer #3, boxplot each interval,
vmin = 0.75; vmax = 2.15;
thewidth = (vmax-vmin)/25;
func_multi_bplot(lab_res(pt_subject_index),mri_res(pt_subject_index,mri_ind),linspace(vmin,vmax,10)+thewidth,0.1,c_deeppink*0.75+0.25,thewidth)
func_multi_bplot(lab_res(ft_subject_index),mri_res(ft_subject_index,mri_ind),linspace(vmin,vmax,10),0.1,c_deepblue*0.75+0.25,thewidth)

xlim([0.5,2.2])

% write back corr value for R GAMM.
albumin_res(1,ft_subject_index) = 1;
albumin_res(1,pt_subject_index) = 0;

albumin_res(2,ft_subject_index) = lab_res(ft_subject_index);
albumin_res(2,pt_subject_index) = lab_res(pt_subject_index);

albumin_res(3,ft_subject_index) = mri_res(ft_subject_index,mri_ind);
albumin_res(3,pt_subject_index) = mri_res(pt_subject_index,mri_ind);

exc_ind = find(min(albumin_res(2:3,:))==0);
albumin_res(:,exc_ind) = [];
writematrix(albumin_res','cysc_vs_aveBV_for_GAMM.csv');

%%
% % layer #2, GAMM fitted curve and CI,
gamm_res = readmatrix('cysc_fitted_GAMM_curve_for_ft.csv');
gamm_res(:,4) =gamm_res(:,4);
gamm_res(:,5) =gamm_res(:,5);
plot(gamm_res(:,1),gamm_res(:,2),'color','b','linewidth',2); hold on;
plot(gamm_res(:,1),gamm_res(:,4),'color','b','linewidth',2,'linestyle','--');
plot(gamm_res(:,1),gamm_res(:,5),'color','b','linewidth',2,'linestyle','--');
fill([gamm_res(:,1)', fliplr(gamm_res(:,1)')], [gamm_res(:,4)', fliplr(gamm_res(:,5)')], c_deepblue*0.75, 'FaceAlpha',0.35 );

gamm_res = readmatrix('cysc_fitted_GAMM_curve_for_Pt.csv');
gamm_res(:,4) =gamm_res(:,4);
gamm_res(:,5) =gamm_res(:,5);
plot(gamm_res(:,1),gamm_res(:,2),'color',c_deeppink*0.5,'linewidth',2);
plot(gamm_res(:,1),gamm_res(:,4),'color',c_deeppink*0.5,'linewidth',2,'linestyle','--');
plot(gamm_res(:,1),gamm_res(:,5),'color',c_deeppink*0.5,'linewidth',2,'linestyle','--');
fill([gamm_res(:,1)', fliplr(gamm_res(:,1)')], [gamm_res(:,4)', fliplr(gamm_res(:,5)')], c_deeppink*0.75, 'FaceAlpha',0.35 );


%%
ylabel('Brain volume Z-Score');
xlabel('Cystatin C (mg/L)')
ylim([-4,4]);

set(findall(gcf,'-property','fontweight'),'fontweight','bold');
set(findall(gcf,'-property','fontsize'),'fontsize',14);
title('Cystatin C','fontsize',16)

%xlim([2,46])
ylim([-3.5,5])
grid on




