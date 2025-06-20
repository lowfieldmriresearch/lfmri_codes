clear; clc; close all;

%% some of the cases do NOT has full serum biomedical profile data,
tt = readtable('../raw_data/full_res_for_gamm00_use.csv');

% how many cases with full serum results,
lab_res = table2array(tt(:,49:84));
good_cases_index = find(~isnan(sum(lab_res,2)).*(sum(lab_res,2)>0));

disp(['',num2str(length(good_cases_index)),' cases had full serum biomedical profile data out of 1245 full data'])

%%
% normal MRI volmetric measurement to z-score based on the estimated
% percentiles in GAMM previously,
% load MRI res from saved Percentile files,
Vdataroot = '../01_growth_curve/saved_GAMMs_K3/';
original_full_table = tt;

% arrange a new full table similar to the input, but MRI volmetric -> z-score
fullresmat = zeros(1245,38);
for ind = 1:47
    vtable = readtable([Vdataroot,'V',num2str(ind+1),'/GAMM_example.csv']);
    vtable = vtable(1:1245,:);
    pt_table = readtable([Vdataroot,'V',num2str(ind+1),'/smoothed_percentiles_pt.csv']);
    ft_table = readtable([Vdataroot,'V',num2str(ind+1),'/smoothed_percentiles_ft.csv']);

    pt_subject_ind = find(vtable{:,2}==0);
    ft_subject_ind = find(vtable{:,2}==1);

    pt_pma = vtable{vtable{:,2}==0,1};
    pt_bv = vtable{vtable{:,2}==0,8};

    ft_pma = vtable{vtable{:,2}==1,1};
    ft_bv = vtable{vtable{:,2}==1,8};

    % z-score for pt,
    ptc_pma = pt_table{:,2};
    ptc_pcurve = pt_table{:,14:18};
    ptc_match_pcurve = interp1(ptc_pma,ptc_pcurve,pt_pma);
    mu = ptc_match_pcurve(:,3);
    sigma1 = (ptc_match_pcurve(:,5) - ptc_match_pcurve(:,1)) / 3.29;   % Using 5th and 95th percentiles
    sigma2 = (ptc_match_pcurve(:,4) - ptc_match_pcurve(:,2)) / 1.348; % Using 25th and 75th percentiles
    sigma12 = (sigma1 + sigma2) / 2;
    pt_zscore = (pt_bv - mu) ./ sigma12;

    % z-score for ft,
    ftc_pma = ft_table{:,2};
    ftc_pcurve = ft_table{:,14:18};
    ftc_match_pcurve = interp1(ftc_pma,ftc_pcurve,ft_pma);
    mu = ftc_match_pcurve(:,3);
    sigma1 = (ftc_match_pcurve(:,5) - ftc_match_pcurve(:,1)) / 3.29;   % Using 5th and 95th percentiles
    sigma2 = (ftc_match_pcurve(:,4) - ftc_match_pcurve(:,2)) / 1.348; % Using 25th and 75th percentiles
    sigma12 = (sigma1 + sigma2) / 2;
    ft_zscore = (ft_bv - mu) ./ sigma12;

    % results of normalized z-score,
    fullresmat(pt_subject_ind,ind+1) = pt_zscore;
    fullresmat(ft_subject_ind,ind+1) = ft_zscore;

    % results of actual brain volumes,
    fullresmat2(pt_subject_ind,ind+1) = pt_bv;
    fullresmat2(ft_subject_ind,ind+1) = ft_bv;


end


%%
fullresmat(:,[1,49:86]) = original_full_table{:,[1,49:86]};
fullresmat2(:,[1,49:86]) = original_full_table{:,[1,49:86]};


% assemble vairable names,
serum_factor_names = table2array(readtable('lab_names_eng_abbr.csv'));
serum_factor_names = serum_factor_names(1:38,:);
brain_region_names = table2array(readtable('../01_growth_curve/V_var_names.csv'));

table_naming_row = [{'pt=0'};brain_region_names;serum_factor_names];
T1 = array2table(fullresmat, 'VariableNames', table_naming_row'); % Create table
writetable(T1, 'mat_for_nutr_analysis_zscore.csv'); % Save as CSV file

T2 = array2table(fullresmat2, 'VariableNames', table_naming_row'); % Create table
writetable(T2, 'mat_for_nutr_analysis_rawvolume.csv'); % Save as CSV file