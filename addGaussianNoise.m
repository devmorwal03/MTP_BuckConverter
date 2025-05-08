% Load the dataset
data = readtable('DataSet.csv');

% List of features to add noise
features_to_modify = {'Voltage_RMS', 'Voltage_Mean', 'Voltage_Peak_Value', ...
                     'Duty_Peak_Value', 'Duty_RMS', 'Duty_Mean', 'Duty_Clearance_factor'};

% Define noise parameters (Gaussian noise: mean = 0, std = 0.1)
noise_mean = 0;
noise_std = 0.1;

% Add Gaussian noise to each feature
for i = 1:length(features_to_modify)
    column_name = features_to_modify{i};
    if ismember(column_name, data.Properties.VariableNames)
        original_data = data.(column_name);
        noise = noise_mean + noise_std * randn(size(original_data));
        data.(column_name) = original_data + noise;  % Add noise to the original data
    else
        disp(['Column not found: ', column_name]);  % Display missing columns
    end
end

% Save the modified dataset with noise added as DataSet10.csv
writetable(data, 'DataSet10.csv');
disp('Gaussian noise added and dataset saved as DataSet10.csv.');
