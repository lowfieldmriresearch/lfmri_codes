import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lightgbm as lgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

# Load your data
raw_data_mat = pd.read_csv('mat_for_nutr_analysis_zscore.csv').to_numpy()
lab_var_name = pd.read_csv('lab_names_eng_abbr.csv', encoding='utf-8-sig', usecols=[0])
# 2. Clean feature names
import re
def clean_feature(name):
    name = str(name).strip()  # remove spaces at beginning/end
    name = re.sub(r'[^A-Za-z0-9_]', '_', name)  # replace weird chars with _
    return name
lab_var_name = lab_var_name.applymap(clean_feature)
lab_var_name = lab_var_name[:36]
feature_names = lab_var_name.iloc[:,0].tolist()
# Check
lab_var_name = feature_names

ptft_tag = raw_data_mat[:, 0]
mri_res = raw_data_mat[:, 1:48]
lab_res = raw_data_mat[:, 48:84]

# data scree, remove the NaN and suspecious outliers,
lab_res[lab_res == 0] = np.nan

# Prepare input
lab_input = pd.DataFrame(lab_res)
import shap
all_shap_values = []
for brInd in range(0,47):
    mri_input = pd.DataFrame(mri_res[:, brInd])  # looks like you're picking the 12th brain region?

    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(lab_input, mri_input, test_size=0.2, random_state=42)

    # Train a LightGBM model
    train_data = lgb.Dataset(X_train, label=y_train)
    test_data = lgb.Dataset(X_test, label=y_test, reference=train_data)

    params = {
        'objective': 'regression',
        'metric': 'l2',  # mean squared error
        'num_leaves': 31,
        'learning_rate': 0.05,
        'verbose': -1,
    }

    model = lgb.train(
        params,
        train_data,
        valid_sets=[train_data, test_data],
        num_boost_round=200,
    )

    # Get feature importance
    importance = model.feature_importance(importance_type='gain')  # same as XGBoost 'total_gain'
    feature_names = model.feature_name()

    # SHAP.
    explainer = shap.Explainer(model)
    shap_values = explainer.shap_values(lab_input)
    all_shap_values.append(shap_values)

    # shap.summary_plot(shap_values, lab_input, feature_names=lab_var_name, max_display=36,show=False)
    # plt.savefig(f"shap_summary_{brInd}.png")  # Save each figure
    # plt.close()  # Then close it
    #
    # shap_df = pd.DataFrame(shap_values, columns=lab_var_name)
    # shap_df.to_csv("shap_values_.csv", index=False)
    #
    # filename = f"shap_values_region#{brInd}.csv"  # e.g., shap_values_1.csv, shap_values_2.csv, etc.
    # save unit.
    #shap_df.to_csv('./LGB_shap_res/'+filename, index=False)

# plot average shap res,
all_shap_values = np.array(all_shap_values)  # Shape: (47, n_samples, n_features)
# Then compute the average across models (i.e., across the 47 regions)
average_shap_values = np.mean(all_shap_values, axis=0)
np.savetxt("average_shap_values.csv", average_shap_values, delimiter=",")

# tune colormap
from matplotlib.colors import LinearSegmentedColormap
custom_colors = ['#A29DD1','#565656' ,'#8EBFB9']  # light gray → teal → purple

custom_cmap = LinearSegmentedColormap.from_list("my_soft_cmap", custom_colors)

fig = plt.figure()
ax = fig.add_subplot(111)
shap.summary_plot(average_shap_values, lab_input,feature_names=lab_var_name,max_display=15,show=False, cmap=custom_cmap)
# Set aspect ratio manually
ax.set_aspect(aspect=0.25)
plt.savefig(f"shap_summary_all.svg", bbox_inches='tight')  # Save each figure
plt.close()  # Then close it

