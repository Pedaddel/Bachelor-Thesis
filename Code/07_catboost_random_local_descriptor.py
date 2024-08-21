import pandas as pd
import numpy as np
import catboost as cb
from pathlib import Path
import ast
from sklearn.metrics import mean_absolute_error, mean_squared_error

df = pd.read_csv(r'/tmp/pycharm_project_396/Processed_data/dataset_show_all_combined_with_ecfp.csv')
df = df[['ee_major', 'combined_local_substrate','combined_local_product','combined_global_substrate','combined_global_product','final_FP']]


# Convert string-formatted lists in 'final_FP' to actual lists
df['final_FP'] = df['final_FP'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

# Flatten the 'final_FP' lists into multiple columns
final_fp_flattened = pd.DataFrame(df['final_FP'].tolist(), index=df.index)
final_fp_flattened.columns = [f"final_FP_{i}" for i in range(final_fp_flattened.shape[1])]

# Drop the original 'final_FP' column and join the flattened columns back to the DataFrame
df = df.drop('final_FP', axis=1).join(final_fp_flattened)

'''
####
# Define a safe literal_eval function with error handling
def safe_literal_eval(x):
    if isinstance(x, str):
        try:
            return ast.literal_eval(x)
        except (ValueError, SyntaxError) as e:
            print(f"Error parsing: {x} -> {e}")
            return None  # or handle as needed (e.g., return an empty list)
    return x

# Apply the safe parsing function to 'combined_local_product'
df['combined_global_product'] = df['combined_global_product'].apply(safe_literal_eval)

# Drop rows where 'combined_local_product' is None if necessary
df = df.dropna(subset=['combined_global_product'])

# Flatten the 'combined_local_product' lists into multiple columns
x_flattened = pd.DataFrame(df['combined_global_product'].tolist(), index=df.index)

# Assign new column names
x_flattened.columns = [f"combined_global_product_{i}" for i in range(x_flattened.shape[1])]

# Drop the original 'combined_local_product' column and join the flattened columns back to the DataFrame
df = df.drop('combined_global_product', axis=1).join(x_flattened)

###
'''

df['combined_global_substrate'] = df['combined_global_substrate'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

# Flatten the 'final_FP' lists into multiple columns
x_flattened = pd.DataFrame(df['combined_global_substrate'].tolist(), index=df.index)
x_flattened.columns = [f"combined_global_substrate_{i}" for i in range(x_flattened.shape[1])]

# Drop the original 'final_FP' column and join the flattened columns back to the DataFrame
df = df.drop('combined_global_substrate', axis=1).join(x_flattened)



# List of dataframes
df= df.sample(frac=1,random_state=6).reset_index(drop=True)
df_new=np.array_split(df, 5)
MAE_test = []
MAE_train = []
RMSE_test = []
RMSE_train = []
results_corr = []
results_st_dev = []
results_mean_corr = []

for i, seed in enumerate([1,2,3,4,5]):
    train = pd.concat([df_new[j] for j in range(5) if j != i])
    test = df.drop(train.index)





    cb_r = cb.CatBoostRegressor(iterations=1000,
                                depth=8)

    train = train.reset_index()
    test = test.reset_index()
    feature_cols = [col for col in train.columns if col.startswith('final_FP')]
    #or col.startswith('final_FP')
    cb_r.fit(train[feature_cols],train['ee_major'])

    predictions_train = cb_r.predict(train[feature_cols])
    true_values_train = train['ee_major'].to_list()
    predictions_test = cb_r.predict(test[feature_cols])
    true_values_test = test['ee_major'].to_list()

    # Calculate MAE
    mae_train = mean_absolute_error(true_values_train, predictions_train)
    mae_test = mean_absolute_error(true_values_test, predictions_test)

    # Calculate RMSE
    rmse_train = np.sqrt(mean_squared_error(true_values_train, predictions_train))
    rmse_test = np.sqrt(mean_squared_error(true_values_test, predictions_test))

    corr = np.corrcoef(predictions_test, true_values_test)[0][1] ** 2

    results_corr.append(corr)
    MAE_train.append(mae_train)
    MAE_test.append(mae_test)
    RMSE_train.append(rmse_train)
    RMSE_test.append(rmse_test)

mean_R2 = np.mean(results_corr)

st_dev = np.std(results_corr)

for _ in range(len(results_corr)):
    results_mean_corr.append(mean_R2)
    results_st_dev.append(st_dev)

# Create a DataFrame with all results
df_results = pd.DataFrame({
    'Results_corr': results_corr,
    'Results_stdev': results_st_dev,
    'Results_mean_corr': results_mean_corr,
    'MAE_train': MAE_train,
    'MAE_test': MAE_test,
    'RMSE_train': RMSE_train,
    'RMSE_test': RMSE_test})

path_to_save = Path(r'/tmp/pycharm_project_823/Calculated_results/eleventh_calculation')
filename_to_save = Path('random_show_peter_final_FP.csv')

combined_path = path_to_save / filename_to_save

df_results.to_csv(combined_path, index=False)

print(MAE_train)
print(MAE_test)
print(RMSE_train)
print(RMSE_test)
print(results_corr)
print(st_dev)
print(mean_R2)