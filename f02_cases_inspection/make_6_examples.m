clear; clc; close all;

tt = readtable('/Users/szx/Desktop/paper2_preparation/2025_mar_rewrite/02_rep_casereview/GA.20250418.xlsx');

load('/Users/szx/Desktop/paper2_preparation/2025_mar_rewrite/02_rep_casereview/full_volumetric_stats.mat');
load('/Users/szx/Desktop/paper2_preparation/2025_mar_rewrite/02_rep_casereview/folder_used_in_curve.mat');
folder_in_table = tt{:,17};
subset_tt = tt;
folder_has_mri = full_volumetric_stats(:,1);



%%

% arrange the selected brain region and corresponding boxplots panels.
% mri brain region grouping: raw_mri_res -> mri_res:
tt = readtable('/Users/szx/Desktop/gp_subfolder/M-CRIB_atlas-master/M-CRIB_labels_ITK_format.txt');

region_name{1} = 'Lateral Ventricle';
region_index{1} = [3,15];
region_name{2} = 'Third Ventricle';
region_index{2} = 8;
region_name{3} = 'Fourth Ventricle';
region_index{3} = 9;
region_name{4} = 'Corpus Callosum';
region_index{4} = 31;
region_name{5} = 'Amyglada';
region_index{5} = [11,21];
region_name{6} = 'White Matter';
region_index{6} = [2,14];
region_name{7} = 'Hippocampus';
region_index{7} = [10,20];
region_name{8} = 'Thalamus';
region_index{8} = [4,16];
region_name{9} = 'Caudate';
region_index{9} = [5,17];
region_name{10} = 'Putamen';
region_index{10} = [6,18];
region_name{11} = 'Pallidum';
region_index{11} = [7,19];
region_name{12} = 'Accumbens';
region_index{12} = [13,22];
region_name{13} = 'Brainstem';
region_index{13} = [30];

for ind = 1:34
    cur_rn = tt{ind+32,8};
    cur_rn = strrep(cur_rn{1},'ctx-lh-','');
    region_name{ind+13} = cur_rn;
    region_index{ind+13} = [ind+32,ind+67];
end

%% anatomical grouping,

% Define the list of brain regions
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

get_indices_in_specified_order = @(names) cellfun(@(x) find(strcmp(x, regions)), names);

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
grouped_index = [categories.Subcortical,categories.Ventricles_CSF,categories.Frontal,categories.Temporal, categories.Parietal, categories.Occipital, categories.Cingulate, categories.Insular_Misc, categories.WhiteMatter];

%%
% extract one measurement using the index,
% load the GAMM estimated percentiles,

% making subcortical part,
var_name = readtable('../01_growth_curve/V_var_names.csv');
num_k = 3;
% load pt percentile,
for ind = 1:47
    % figure part shape,
    cur_V_ind = ind+1;
    cur_region_name = var_name{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    pt_perc_data = readtable(['../01_growth_curve/saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_pt.csv']);
    conf_cd(:,1) =  pt_perc_data{:,2};
    conf_cd(:,2:6) =  pt_perc_data{:,14:18};
    PT_CONF_CD{ind} = conf_cd;
    clear conf_cd;
end

% load ft percentile,
for ind = 1:47
    % figure part shape,
    cur_V_ind = ind+1;
    cur_region_name = var_name{cur_V_ind-1,1};
    cur_region_name = cur_region_name{1,1};
    % load GAMM fitting resutls,
    pt_perc_data = readtable(['../01_growth_curve/saved_GAMMs_K',num2str(num_k),'/V',num2str(cur_V_ind),'/smoothed_percentiles_ft.csv']);
    conf_cd(:,1) =  pt_perc_data{:,2};
    conf_cd(:,2:6) =  pt_perc_data{:,14:18};
    FT_CONF_CD{ind} = conf_cd;
    clear conf_cd;
end

%%
search_root = '/Volumes/WD_BLACK/lfmri_data_proc/finished_data/';
close all;

sel_ind_list = [28,83,282,48,256,106];
k = 1;
for ind = sel_ind_list
    % this step find the raw data,
    cur_folder_name = subset_tt{ind,17};
    found_index = find(strcmp(cur_folder_name,folder_has_mri));

    %
    % also I want to extract the QC quickpic,
    fn1 = dir([search_root,'*/',cur_folder_name{1},'/T2_AX/qc/qc_poriors.png']);
    fn2 = dir([search_root,'*/*/',cur_folder_name{1},'/T2_AX/qc/qc_poriors.png']);
    fn = cat(1,fn1,fn2);
    if(length(fn)>0)
        mkdir(['./qcs/Index_',num2str(ind,'%.3d'),'_',cur_folder_name{1}]);
        qc_img_path = [fn(1).folder,filesep,fn(1).name];
        copyfile(qc_img_path,['./qcs/Index_',num2str(ind,'%.3d'),'_',cur_folder_name{1},'/qc.png']);
    end
    % maxprob seg,
    fn1 = dir([search_root,'*/',cur_folder_name{1},'/T2_AX/temp2ind/def_maxprob_hardseg.nii.gz']);
    fn2 = dir([search_root,'*/*/',cur_folder_name{1},'/T2_AX/temp2ind/def_maxprob_hardseg.nii.gz']);
    fn = cat(1,fn1,fn2);
    if(length(fn)>0)
        mkdir(['./qcs/Index_',num2str(ind,'%.3d'),'_',cur_folder_name{1}]);
        qc_img_path = [fn(1).folder,filesep,fn(1).name];
        copyfile(qc_img_path,['./qcs/Index_',num2str(ind,'%.3d'),'_',cur_folder_name{1},'/def_maxprob_hardseg.nii.gz']);
    end

    % raw img,
    fn1 = dir([search_root,'*/',cur_folder_name{1},'/T2_AX/RAW/*interp*masked*vol*.nii.gz']);
    fn2 = dir([search_root,'*/*/',cur_folder_name{1},'/T2_AX/RAW/*interp*masked*vol*.nii.gz']);
    fn = cat(1,fn1,fn2);
    if(length(fn)>0)
        mkdir(['./qcs/Index_',num2str(ind,'%.3d'),'_',cur_folder_name{1}]);

        qc_img_path = [fn(1).folder,filesep,fn(1).name];
        copyfile(qc_img_path,['./qcs/Index_',num2str(ind,'%.3d'),'_',cur_folder_name{1},'/interp_10x_masked_volume.nii.gz']);
    end


    %
    cur_rawmri_measures = full_volumetric_stats{found_index,2};
    if(subset_tt{ind,20} == 0)
        disp('This one no GA record??')
        continue;
    end
    GA = subset_tt{ind,20};
    PMA = days(datetime(subset_tt{ind,15}) - datetime(subset_tt{ind,8}))/7 + subset_tt{ind,20};

    GA_list(ind) = GA;
    PMA_list(ind) = PMA;

    if(GA<37)
        CONF_CD = PT_CONF_CD;
        disp(['Cur GA = ',num2str(GA),' Preterm'])
    else
        CONF_CD = FT_CONF_CD;
        disp(['Cur GA = ',num2str(GA),' Term'])
    end
    m47 = func_101to47_regions(cur_rawmri_measures);

    M47{ind} = m47;

    % compute the percentile of the current m47 data,
    for brInd = 1:47
        conf_cd = CONF_CD{brInd};
        p5 = interp1(conf_cd(:,1),conf_cd(:,2),PMA);
        p25 = interp1(conf_cd(:,1),conf_cd(:,3),PMA);
        p50 = interp1(conf_cd(:,1),conf_cd(:,4),PMA);
        p75 = interp1(conf_cd(:,1),conf_cd(:,5),PMA);
        p95 = interp1(conf_cd(:,1),conf_cd(:,6),PMA);

        % Estimate mean and std (assuming normal distribution)
        mu = p50;
        sigma1 = (p95 - p5) / (2 * 1.645);
        sigma2 = (p75 - p25) / 1.34896;
        sigma = (sigma1+sigma2)/2;

        % Compute z-score
        cur_z_score = (m47(brInd) - mu) / sigma;
        z_score_list(ind,brInd) = cur_z_score;
    end

    % make example inputs,
    example_tt = readtable('/Users/szx/Desktop/paper2_preparation/2025_mar_rewrite/02_rep_casereview/brain-age-app2/example_input.csv');
    example_tt{1,1} = PMA;
    example_tt{1,2:48} = m47;

    writetable(example_tt,['/Users/szx/Desktop/paper2_preparation/2025_mar_rewrite/02_rep_casereview/brain-age-app2/example_inputs/example_',num2str(k),'.csv']);

    k=k+1;



end

%% make figure,

% make nature colorbar,
load('/Users/szx/Desktop/paper2_preparation/mcrib_ver2/full_brain_3d_model_LR.mat');
MODELS_LR = MODELS;
load('/Users/szx/Desktop/paper2_preparation/mcrib_ver2/full_brain_3d_model.mat');
MODELS_L = MODELS;
thiscolormap = read_json_colormap('percentiles_bluered.json');

sel_ind_list = [28,83,282,48,256,106];
for ind_inspect = 106
    kk=1;
    close all;

    % show the bar plot style2.
    figure('position',[35 421 1403 454])
    showcase_z_score = z_score_list(ind_inspect,:);
    % inspect one example,
    cur_percentile = normcdf(showcase_z_score*kk)* 100;
    b = bar(cur_percentile(grouped_index)-50);
    for ind = 1:length(cur_percentile)
        barColors(ind,:) = thiscolormap(ceil(cur_percentile(grouped_index(ind))),:);
    end
    % Assign colors to each bar
    b.FaceColor = 'flat';         % important: enable per-bar coloring
    b.CData = barColors;
    ylabel('Volume percentiles')
    xlabel('')
    colormap(thiscolormap)
    clim([0 100])
    set(findall(gcf,'-property','fontweight'),'fontweight','bold');
    set(findall(gcf,'-property','fontsize'),'fontsize',18);
    ylim([-50,50]); xlim([0,48]);
    set(gcf, 'Renderer', 'painters');
    xticks(1:47)
    xticklabels(regions(grouped_index))
    xtickangle(45);  % Sets the x-tick labels to 45 degrees 
    yticks([-50,0,50]);
    yticklabels({'0th','50th','100th'})

    set(findall(gcf,'-property','fontweight'),'fontweight','bold');
    set(findall(gcf,'-property','fontsize'),'fontsize',16);
    saveas(gcf,'test1.png')
    %
    % show 3d brain modelres,pt1
    figure('position',[ 158         814        1163         454])
    ha = tight_subplot(1,3,0,0,0);
    axes(ha(1))
    ccmap = thiscolormap;
    %ccmap = flip(ccmap,1);
    for ind = 14:47
        fv = MODELS_L{ind};
        if(cur_percentile(ind)>10)
        else
            patch(fv,'edgecolor','none','facecolor',ccmap(ceil(cur_percentile(ind)),:),'SpecularStrength', 0.1);
        end
    end
    fv = MODELS_LR{48};
    patch(fv,'edgecolor','none','facecolor',[1,1,1]*0.5,'facealpha',0.13);

    axis equal; axis off;
    view([163,8.8]);
    light

    %
    axes(ha(3))

    for ind = 14:47
        fv = MODELS_L{ind};
        if(cur_percentile(ind)>10)
        else
            patch(fv,'edgecolor','none','facecolor',ccmap(ceil(cur_percentile(ind)),:),'SpecularStrength', 0.1);
        end
    end
    fv = MODELS_LR{48};
    patch(fv,'edgecolor','none','facecolor',[1,1,1]*0.5,'facealpha',0.13);
    axis equal; axis off;

    view([33,10]);
    light


    %
    axes(ha(2))

    for ind = 1:13
        fv = MODELS_LR{ind};
        if(cur_percentile(ind)>10)
        else
            patch(fv,'edgecolor','none','facecolor',ccmap(ceil(cur_percentile(ind)),:),'SpecularStrength', 0.1);
        end
    end
    fv = MODELS_LR{48};
    patch(fv,'edgecolor','none','facecolor',[1,1,1]*0.5,'facealpha',0.13);
    axis equal; axis off;
    view([100.3,7.8]);
    light


    % also write example input:
    example_tt = readtable('/Users/szx/Desktop/paper2_preparation/2025_mar_rewrite/02_rep_casereview/brain-age-app2/example_input.csv');
    example_tt{1,1} = days(datetime(subset_tt{ind_inspect,15}) - datetime(subset_tt{ind_inspect,8}))/7 + subset_tt{ind_inspect,20};






end
%  9,3,9,9,6,4,4,2,1
% ticks at
xticks([1,9,10,12,13,21,22,30,31,36,37,40,41,44,45,46]);

%% select 3 representative region to show the growth curve,
% cur_percentile(grouped_index),
% close all;
% COLOR PROFILE,
c_deepblue = [87,111,160]/256;
c_lightblue = [167,185,215]/256;
c_deeppink = [181,121,121]/256;
c_lightpink = [222,163,162]/256;

show_ind_list = [1,19,40];
figure('position',[48,813,1614,600])
for ind = 1:3
    subplot(3,6,ind);
    cur_show_region = grouped_index(show_ind_list(ind));

    PMA = PMA_list(ind_inspect);
    GA = GA_list(ind_inspect);

    if(GA<37)
        CONF_CD = PT_CONF_CD;
        disp(['Cur GA = ',num2str(GA),' Preterm'])
        c_color = c_deeppink;
    else
        CONF_CD = FT_CONF_CD;
        disp(['Cur GA = ',num2str(GA),' Term'])
        c_color = c_deepblue;
    end
    m47 = M47{ind_inspect};
    bv = m47(cur_show_region);

    cur_curve = CONF_CD{cur_show_region};

    pt_pma = cur_curve(:,1);
    pt_pert = cur_curve(:,2:6);

    fill([pt_pma; flipud(pt_pma)], ...
        [pt_pert(:,1); flipud(pt_pert(:,5))], ...
        c_color, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); hold on;

    % Shaded area for preterm: 25–75%
    fill([pt_pma; flipud(pt_pma)], ...
        [pt_pert(:,2); flipud(pt_pert(:,4))], ...
        c_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none'); hold on;


    for ind = 1:5
        plot(pt_pma,pt_pert(:,ind),'color',c_color,'linewidth',0.8); hold on;
    end
    plot(pt_pma,pt_pert(:,3),'color',c_color,'linewidth',1.5); hold on;

    grid
    max_value = max(pt_pert(:));
    ax = gca; ax.YAxis.Exponent = floor(log10(abs(max_value)));

    xlim([31,50])
    %xlim([36,54])

    scatter(PMA,bv,80,'filled','diamond')
    set(findall(gcf,'-property','fontsize'),'fontsize',13);
    set(findall(gcf,'-property','fontweight'),'fontweight','bold');

    title(regions{cur_show_region})
    
end
