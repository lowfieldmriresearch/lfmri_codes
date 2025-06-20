clear; clc; close all;

d_shap = readmatrix('average_shap_values.csv');
ave_shap = mean(abs(d_shap),1);


var_name = readtable('lab_names_eng_abbr.csv');

[B,I] = sort(ave_shap,'Descend');

bar(ave_shap(I),0.65,'facealpha',0.3)
xticks(1:36);
xticklabels(var_name{I,1})
xlim([0,15.5])
