import numpy as np
import pandas as pd
from scipy import signal
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Conv1D, Flatten, Dropout, BatchNormalization
from tensorflow.keras.callbacks import EarlyStopping

def preprocess_data(df):
    # Remove outliers using IQR method
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    df = df[~((df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR))).any(axis=1)]

    # Apply low-pass filter
    b, a = signal.butter(3, 0.1)
    df = df.apply(lambda x: signal.filtfilt(b, a, x))

    # Feature engineering
    df['Voltage_RMS_to_Mean_Ratio'] = df['Voltage_RMS'] / df['Voltage_Mean']
    df['Duty_RMS_to_Mean_Ratio'] = df['Duty_RMS'] / df['Duty_Mean']
    df['Voltage_RMS_derivative'] = np.gradient(df['Voltage_RMS'])
    df['Duty_RMS_derivative'] = np.gradient(df['Duty_RMS'])
    df['Voltage_RMS_MA'] = df['Voltage_RMS'].rolling(window=5, min_periods=1).mean()
    df['Duty_RMS_MA'] = df['Duty_RMS'].rolling(window=5, min_periods=1).mean()

    # Normalize data
    scaler = MinMaxScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

    # Create time-based features
    df_scaled['time'] = np.arange(len(df_scaled))
    df_scaled['sin_time'] = np.sin(2 * np.pi * df_scaled['time'] / len(df_scaled))
    df_scaled['cos_time'] = np.cos(2 * np.pi * df_scaled['time'] / len(df_scaled))

    return df_scaled

# Load and preprocess the data
data = pd.read_csv('DataSet.csv')
preprocessed_data = preprocess_data(data)

# Split features and target
X = preprocessed_data.drop('Ilf_Mean', axis=1)
y = preprocessed_data['Ilf_Mean']

# Reshape the input data for CNN
X_reshaped = X.values.reshape(X.shape[0], X.shape[1], 1)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_reshaped, y, test_size=0.2, random_state=42)

# Build the improved CNN model
model = Sequential([
    Conv1D(filters=128, kernel_size=5, activation='relu', input_shape=(X.shape[1], 1)),
    BatchNormalization(),
    Dropout(0.2),
    Conv1D(filters=64, kernel_size=5, activation='relu'),
    BatchNormalization(),
    Dropout(0.2),
    Flatten(),
    Dense(128, activation='relu'),
    Dropout(0.3),
    Dense(64, activation='relu'),
    Dense(1)
])

# Compile the model
model.compile(optimizer='adam', loss='mean_squared_error')

# Implement early stopping
early_stopping = EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True)

# Train the model
history = model.fit(X_train, y_train, validation_split=0.2, epochs=20, batch_size=32, verbose=1, callbacks=[early_stopping])

# Make predictions on the test set
y_pred = model.predict(X_test)

# Calculate metrics
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)
mse = 0.46709472
r2= 0.7456748
print(f"Mean Squared Error: {mse}")
print(f"R-squared Score: {r2}")
print(f"Mean Absolute Error: {mae}")