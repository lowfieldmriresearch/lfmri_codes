function m47 = func_101to47_regions(m101)

%% brain region arrangement based on mCRIB ver1.
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

for ind = 1:47
    m47(ind) = sum(m101(region_index{ind}));
end

end