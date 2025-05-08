import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt

# Load and preprocess the data
def preprocess_data(df):
    # Remove outliers using IQR method
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    df = df[~((df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR))).any(axis=1)]

    # Feature engineering
    df['Voltage_RMS_to_Mean_Ratio'] = df['Voltage_RMS'] / df['Voltage_Mean']
    df['Duty_RMS_to_Mean_Ratio'] = df['Duty_RMS'] / df['Duty_Mean']
    df['Voltage_Duty_Interaction'] = df['Voltage_Mean'] * df['Duty_Mean']
    # df['Power_Feature'] = df['Voltage_Mean'] * df['Ilf_Mean']

    # Calculate derivatives
    df['Voltage_RMS_derivative'] = np.gradient(df['Voltage_RMS'])
    df['Duty_RMS_derivative'] = np.gradient(df['Duty_RMS'])

    # Normalize data
    scaler = MinMaxScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

    # Create time-based features
    df_scaled['time'] = np.arange(len(df_scaled))
    df_scaled['sin_time'] = np.sin(2 * np.pi * df_scaled['time'] / len(df_scaled))
    df_scaled['cos_time'] = np.cos(2 * np.pi * df_scaled['time'] / len(df_scaled))

    return df_scaled

# Load the data
data = pd.read_csv('DataSetfinal.csv')
preprocessed_data = preprocess_data(data)

# Split features and target
X = preprocessed_data.drop('Ilf_Mean', axis=1)
y = preprocessed_data['Ilf_Mean']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Create and train the Random Forest model
rf_model = RandomForestRegressor(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)

# Make predictions
y_pred = rf_model.predict(X_test)

# Calculate MSE and R-squared
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Mean Squared Error: {mse}")
print(f"R-squared Score: {r2}")

# Feature importance
feature_importance = pd.DataFrame({'feature': X.columns, 'importance': rf_model.feature_importances_})
feature_importance = feature_importance.sort_values('importance', ascending=False)

# Plot feature importance
plt.figure(figsize=(10, 6))
plt.bar(feature_importance['feature'][:10], feature_importance['importance'][:10])
plt.title('Top 10 Feature Importances')
plt.xlabel('Features')
plt.ylabel('Importance')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()

# Scatter plot of predicted vs actual values
plt.figure(figsize=(10, 6))
plt.scatter(y_test, y_pred, alpha=0.5)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.title('Actual vs Predicted Ilf_Mean')
plt.tight_layout()
plt.show()
