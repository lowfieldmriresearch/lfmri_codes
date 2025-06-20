clear; clc; close all;

region_names = readtable('V_var_names2.csv');
p_res = readtable('ANOVA_results.csv');
num_k = 3;

% COLOR PROFILE,
% color profiling:
c_deepblue = [87,111,160]/256;
c_lightblue = [167,185,215]/256;
c_deeppink = [181,121,121]/256;
c_lightpink = [222,163,162]/256;
%% anatomical grouping,

regions = { ...
    'Lateral Ventricle', 'Third Ventricle', 'Fourth Ventricle', 'Corpus Callosum', ...
    'Amygdala', 'White Matter', 'Hippocampus', 'Thalamus', 'Caudate', 'Putamen', ...
    'Pallidum', 'Accumbens', 'Brainstem', 'bankssts', 'caudalanteriorcingulate', ...
    'caudalmiddlefrontal', 'cuneus', 'entorhinal', 'fusiform', 'inferiorparietal', ...
    'inferiortemporal', 'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal', ...
    'lingual', 'medialorbitofrontal', 'middletemporal', 'parahippocampal', 'paracentral', ...
    'parsopercularis', 'parsorbitalis', 'parstriangularis', 'pericalcarine', 'postcentral', ...
    'posteriorcingulate', 'precentral', 'precuneus', 'rostralanteriorcingulate', 'rostralmiddlefrontal', ...
    'superiorfrontal', 'superiorparietal', 'superiortemporal', 'supramarginal', 'frontalpole', ...
    'temporalpole', 'transversetemporal', 'insula' ...
};

% Helper function that preserves the input order
get_indices_in_specified_order = @(names) cellfun(@(x) find(strcmp(x, regions)), names);

% Define the categories in your specified order
% Define the categories in order: based on the anatomical groupings,
categories = struct( ...
    'Subcortical', get_indices_in_specified_order({'Amygdala', 'Hippocampus', 'Accumbens','Pallidum', 'Putamen', 'Caudate','Thalamus', ...
         'Brainstem', 'Corpus Callosum'}), ...
    'Ventricles_CSF', get_indices_in_specified_order({'Lateral Ventricle', 'Third Ventricle', 'Fourth Ventricle'}), ...
    'Frontal', get_indices_in_specified_order({'caudalmiddlefrontal', 'rostralmiddlefrontal', 'superiorfrontal', ...
        'parsopercularis', 'parsorbitalis', 'parstriangularis', ...
        'lateralorbitofrontal', 'medialorbitofrontal', 'frontalpole','precentral'}), ...
    'Temporal', get_indices_in_specified_order({'middletemporal', 'inferiortemporal', 'superiortemporal', ...
        'entorhinal', 'fusiform', 'parahippocampal', 'temporalpole', ...
        'transversetemporal', 'bankssts'}), ...
    'Parietal', get_indices_in_specified_order({'inferiorparietal', 'superiorparietal', 'supramarginal', ...
        'postcentral', 'paracentral','precuneus'}), ...
    'Occipital', get_indices_in_specified_order({'cuneus', 'lingual', 'pericalcarine', 'lateraloccipital'}), ...
    'Cingulate', get_indices_in_specified_order({'caudalanteriorcingulate', 'rostralanteriorcingulate', ...
        'posteriorcingulate', 'isthmuscingulate'}), ...
    'Insular_Misc', get_indices_in_specified_order({'insula'}), ...
    'WhiteMatter', get_indices_in_specified_order({'White Matter'}) ...
);

%%
% making subcortical part,
cur_working_index = categories.Subcortical;
figure('position',[48,813,1614,600])
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,ind)
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
    
end
set(gcf, 'Renderer', 'painters');
saveas(gcf,'./fig_in_paper/fig1.fig')

% making CSF part,
cur_working_index = categories.Ventricles_CSF;
figure('position',[48,813,1614,600])
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,ind)
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
end
saveas(gcf,'./fig_in_paper/fig2.fig')

%%
% making Frontel part,
cur_working_index = categories.Frontal;
figure('position',[48,813,1614,600])
panel_ind_list = [1,2,3,4,7,8,9,13,14,15];
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,panel_ind_list(ind))
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
end
set(gcf, 'Renderer', 'painters');
saveas(gcf,'./fig_in_paper/fig3.fig')

%%
% making Temporal part,
cur_working_index = categories.Temporal;
figure('position',[48,813,1614,600])
panel_ind_list = [1,2,3,7,8,9,13,14,15];
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,panel_ind_list(ind))
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
end
set(gcf, 'Renderer', 'painters');
saveas(gcf,'./fig_in_paper/fig4.fig')

% making Parietal part,
cur_working_index = categories.Parietal;
figure('position',[48,813,1614,600])
panel_ind_list = [1:6];
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,panel_ind_list(ind))
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
end
set(gcf, 'Renderer', 'painters');
saveas(gcf,'./fig_in_paper/fig5.fig')

% making Occipital part,
cur_working_index = categories.Occipital;
figure('position',[48,813,1614,600])
panel_ind_list = [1,2,3,4];
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,panel_ind_list(ind))
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
end
set(gcf, 'Renderer', 'painters');
saveas(gcf,'./fig_in_paper/fig6.fig')


% making Cingulate part,
cur_working_index = categories.Cingulate;
figure('position',[48,813,1614,600])
panel_ind_list = [1,2,3,4];
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,panel_ind_list(ind))
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
end
saveas(gcf,'./fig_in_paper/fig7.fig')

% making Insular_Misc part,
cur_working_index = categories.Insular_Misc;
figure('position',[48,813,1614,600])
panel_ind_list = [1,2];
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,panel_ind_list(ind))
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
end
set(gcf,'color','w')
saveas(gcf,'./fig_in_paper/fig8.fig')

% making Insular_Misc part,
cur_working_index = categories.WhiteMatter;
figure('position',[48,813,1614,600])
panel_ind_list = [1];
for ind = 1:length(cur_working_index)
    % figure part shape,
    subplot(3,6,panel_ind_list(ind))
    cur_V_ind = cur_working_index(ind) + 1;
    cur_region_name = region_names{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    conf_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_example.csv']);
    ft_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    pt_perc_data = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    cur_btwgp_pvalue = readtable(['./saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/GAMM_group_pvalue.csv']);
    % run it.
    inter_code_1
end
set(gcf,'color','w')
saveas(gcf,'./fig_in_paper/fig9.fig')
