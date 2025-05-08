
%% new code
%% Load full data from MAT files
data_voltage = load('final_voltage.mat');  % Load everything into the 'data_voltage' structure
final_voltage = data_voltage.final_voltage;  % Access final_voltage from the structure

data_duty = load('final_duty.mat');  % Load everything into the 'data_duty' structure
final_duty = data_duty.final_duty;  % Access final_duty from the structure

% Check the size of the data
numSubarrays_voltage = length(final_voltage);
numSubarrays_duty = length(final_duty);

disp(['Number of voltage subarrays: ', num2str(numSubarrays_voltage)]);
disp(['Number of duty subarrays: ', num2str(numSubarrays_duty)]);

%% Voltage Feature Extraction
Voltage_RMS = zeros(numSubarrays_voltage, 1);
Voltage_Mean = zeros(numSubarrays_voltage, 1);
Voltage_Peak_Value = zeros(numSubarrays_voltage, 1);

for i = 1:numSubarrays_voltage
    if size(final_voltage{i}, 1) < 486000
        disp(['Voltage Subarray ', num2str(i), ' has less than 486000 rows. Skipping.']);
        continue;
    end

    voltage_segment = final_voltage{i}(486000:end);

    % Feature Calculations
    Voltage_RMS(i) = rms(voltage_segment)*10^4;
    Voltage_Mean(i) = mean(voltage_segment)*10^4;
    Voltage_Peak_Value(i) = max(voltage_segment)*10^4;
end

%% Create Voltage Feature Table
valid_indices_voltage = Voltage_RMS ~= 0;  % Ensure no zero entries
voltage_data = table(Voltage_RMS(valid_indices_voltage), Voltage_Mean(valid_indices_voltage), ...
    Voltage_Peak_Value(valid_indices_voltage), ...
    'VariableNames', {'Voltage_RMS', 'Voltage_Mean', 'Voltage_Peak_Value'});

% Save the extracted voltage features
writetable(voltage_data, 'voltage_features.csv');
disp('Voltage feature extraction complete.');

%% Duty Feature Extraction
Duty_RMS = zeros(numSubarrays_duty, 1);
Duty_Mean = zeros(numSubarrays_duty, 1);
Duty_Peak_Value = zeros(numSubarrays_duty, 1);
Duty_Clearance_factor = zeros(numSubarrays_duty, 1);

for i = 1:numSubarrays_duty
    if size(final_duty{i}, 1) < 101222
        disp(['Duty Subarray ', num2str(i), ' has less than 101222 rows. Skipping.']);
        continue;
    end

    duty_segment = final_duty{i}(101222:end);

    % Feature Calculations
    Duty_RMS(i) = rms(duty_segment);
    Duty_Mean(i) = mean(duty_segment);
    Duty_Peak_Value(i) = max(duty_segment);
    Duty_Clearance_factor(i) = Duty_Peak_Value(i) / Duty_RMS(i);
end

%% Create Duty Feature Table
valid_indices_duty = Duty_RMS ~= 0;  % Ensure no zero entries
duty_feature_table = table(Duty_Peak_Value(valid_indices_duty), Duty_RMS(valid_indices_duty), ...
    Duty_Mean(valid_indices_duty), Duty_Clearance_factor(valid_indices_duty), ...
    'VariableNames', {'Duty_Peak_Value', 'Duty_RMS', 'Duty_Mean', 'Duty_Clearance_factor'});

% Save the extracted duty features
writetable(duty_feature_table, 'duty_features.csv');
disp('Duty feature extraction complete.');

%% Merge with IlfArray
load('IlfFinal.mat', 'IlfArray');

IlfTable = array2table(IlfArray, 'VariableNames', {'Ilf_Mean'});

% Ensure same number of rows between IlfTable and duty_feature_table
numRows = min(height(IlfTable), height(duty_feature_table));  % Find the minimum number of rows
IlfTable = IlfTable(1:numRows, :);
duty_feature_table = duty_feature_table(1:numRows, :);

% Merge and save the final dataset
merged_table = [duty_feature_table, IlfTable];
writetable(merged_table, 'DataSet6.csv');
disp('Final dataset saved in DataSet6.csv.');

%% Normalize Voltage Features and Other Data
% Using Min-Max scaling to bring values to the range [0, 1]
% You can also use z-score standardization if preferred.

% Apply normalization to the numeric columns of the table
numeric_columns = {'Voltage_RMS', 'Voltage_Mean', 'Voltage_Peak_Value', ...
                   'Duty_Peak_Value', 'Duty_RMS', 'Duty_Mean', 'Duty_Clearance_factor', ...
                    'Ilf_Mean'};

% Combine voltage data with duty features before normalizing
merged_table = [voltage_data, merged_table];



%% Method 1
% % Normalize using Min-Max scaling
% for col = numeric_columns
%     column_name = col{1};
%     if ~strcmp(column_name, 'Ilf_Mean')  % Don't normalize Ilf_Mean
%         col_data = merged_table.(column_name);
%         min_val = min(col_data);
%         max_val = max(col_data);
%         merged_table.(column_name) = (col_data - min_val) / (max_val - min_val);
%     end
% end
% 
% % Save the normalized dataset
% writetable(merged_table, 'DataSet7.csv');
% disp('Normalized final dataset saved in DataSet7.csv.');



%% Method 2 of normalization
% Ensure no duplicate column names
voltage_data = removevars(voltage_data, intersect(voltage_data.Properties.VariableNames, merged_table.Properties.VariableNames));

% Merge voltage_data and merged_table
merged_table = [voltage_data, merged_table];

%% Normalize using L2 normalization
numeric_columns = {'Voltage_RMS', 'Voltage_Mean', 'Voltage_Peak_Value'};

for col = numeric_columns
    column_name = col{1};
    col_data = merged_table.(column_name);
    l2_norm = sqrt(sum(col_data .^ 2));  % Calculate L2 norm of the column
    if l2_norm ~= 0
        merged_table.(column_name) = col_data / l2_norm;  % Normalize the column
    end
end

% Save the normalized dataset
writetable(merged_table, 'DataSet.csv');
disp('L2-normalized final dataset saved in DataSet.csv.');
