import json

def check_integrity(dataframe, log):
    df = dataframe.copy()
    
    for col in df.columns:
        assert not df[col].isnull().all()

    column_data = {}
    for col in df.columns:
        nan_percentage = df[col].isnull().sum() / len(df)
        notnan_percentage = 1 - nan_percentage
        
        non_nan_values = df[col].dropna().value_counts(normalize=True).to_dict()
        
        column_data[col] = {
            'NaN_percentage': nan_percentage,
            'notNaN_percentage': notnan_percentage,
            'unique_values': non_nan_values
        }

    with open("Data/column_data.json", "w") as f:
        json.dump(column_data, f, indent=4)

