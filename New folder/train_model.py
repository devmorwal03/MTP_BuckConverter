import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error, r2_score
import joblib
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dataset
df = pd.read_csv('DataSetfinal.csv')

# Separate features and target
X = df.drop('Ilf_Mean', axis=1)
y = df['Ilf_Mean']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Scale the features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Add polynomial features
poly = PolynomialFeatures(degree=2, include_bias=False)
X_train_poly = poly.fit_transform(X_train_scaled)
X_test_poly = poly.transform(X_test_scaled)

# Initialize models with their parameter grids
models = {
    'Random Forest': {
        'model': RandomForestRegressor(random_state=42),
        'params': {
            'n_estimators': [100, 200, 300],
            'max_depth': [None, 10, 20, 30],
            'min_samples_split': [2, 5, 10]
        }
    },
    'Gradient Boosting': {
        'model': GradientBoostingRegressor(random_state=42),
        'params': {
            'n_estimators': [100, 200, 300],
            'learning_rate': [0.01, 0.1, 0.2],
            'max_depth': [3, 5, 7]
        }
    },
    'Extra Trees': {
        'model': ExtraTreesRegressor(random_state=42),
        'params': {
            'n_estimators': [100, 200, 300],
            'max_depth': [None, 10, 20, 30]
        }
    },
    'SVR': {
        'model': SVR(),
        'params': {
            'C': [0.1, 1, 10],
            'gamma': ['scale', 'auto', 0.1, 0.01]
        }
    },
    'Ridge': {
        'model': Ridge(),
        'params': {
            'alpha': [0.1, 1, 10, 100]
        }
    }
}

# Train and evaluate models
results = {}
best_model = None
best_score = float('-inf')

for name, model_info in models.items():
    print(f"\nTraining {name}...")
    
    # Perform grid search with cross-validation
    grid_search = GridSearchCV(
        model_info['model'],
        model_info['params'],
        cv=5,
        scoring='r2',
        n_jobs=-1
    )
    
    # Train on both scaled and polynomial features
    grid_search.fit(X_train_scaled, y_train)
    grid_search_poly = GridSearchCV(
        model_info['model'],
        model_info['params'],
        cv=5,
        scoring='r2',
        n_jobs=-1
    )
    grid_search_poly.fit(X_train_poly, y_train)
    
    # Get best model from each approach
    best_model_scaled = grid_search.best_estimator_
    best_model_poly = grid_search_poly.best_estimator_
    
    # Evaluate both approaches
    y_pred_scaled = best_model_scaled.predict(X_test_scaled)
    y_pred_poly = best_model_poly.predict(X_test_poly)
    
    mse_scaled = mean_squared_error(y_test, y_pred_scaled)
    r2_scaled = r2_score(y_test, y_pred_scaled)
    
    mse_poly = mean_squared_error(y_test, y_pred_poly)
    r2_poly = r2_score(y_test, y_pred_poly)
    
    # Choose the better approach
    if r2_scaled > r2_poly:
        results[name] = {
            'MSE': mse_scaled,
            'R2': r2_scaled,
            'model': best_model_scaled,
            'features': 'scaled'
        }
    else:
        results[name] = {
            'MSE': mse_poly,
            'R2': r2_poly,
            'model': best_model_poly,
            'features': 'polynomial'
        }
    
    # Update best model
    if results[name]['R2'] > best_score:
        best_score = results[name]['R2']
        best_model = results[name]['model']

# Print results
print("\nModel Performance:")
print("-" * 50)
for name, metrics in results.items():
    print(f"\n{name}:")
    print(f"MSE: {metrics['MSE']:.4f}")
    print(f"R2 Score: {metrics['R2']:.4f}")
    print(f"Features used: {metrics['features']}")

# Save the best model and scaler
joblib.dump(best_model, 'best_model.joblib')
joblib.dump(scaler, 'scaler.joblib')
if best_model == results['Random Forest']['model']:
    joblib.dump(poly, 'poly_features.joblib')

# Plot actual vs predicted values for the best model
y_pred_best = best_model.predict(X_test_scaled if results[list(results.keys())[0]]['features'] == 'scaled' else X_test_poly)
plt.figure(figsize=(10, 6))
plt.scatter(y_test, y_pred_best, alpha=0.5)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
plt.xlabel('Actual Ilf_Mean')
plt.ylabel('Predicted Ilf_Mean')
plt.title('Actual vs Predicted Ilf_Mean')
plt.savefig('prediction_plot.png')
plt.close()

# Feature importance for tree-based models
if isinstance(best_model, (RandomForestRegressor, GradientBoostingRegressor, ExtraTreesRegressor)):
    feature_importance = pd.DataFrame({
        'Feature': X.columns,
        'Importance': best_model.feature_importances_
    }).sort_values('Importance', ascending=False)
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Importance', y='Feature', data=feature_importance)
    plt.title('Feature Importance')
    plt.tight_layout()
    plt.savefig('feature_importance.png')
    plt.close()

# Correlation heatmap
plt.figure(figsize=(12, 8))
correlation_matrix = df.corr()
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0)
plt.title('Feature Correlation Matrix')
plt.tight_layout()
plt.savefig('correlation_heatmap.png')
plt.close() 