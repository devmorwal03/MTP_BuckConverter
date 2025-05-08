%% 
clc
close all
clear all

%% Converter parameters
% Operating point
Rload = 2: 0.5: 100; % Load resistances
Uout = 140: 5: 180; % Output voltages
% Rload = 2;
% Uout = 140;
tStart = 0.3;
tEnd = 0.4;

%% Code for auto simulation
lenLoadVal = length(Rload);
lenVoltage = length(Uout);
totalNumCases = lenVoltage * lenLoadVal;

% Define table columns
IlfArray = zeros(totalNumCases, 1);
Ilf_preArray = zeros(totalNumCases, 1);
dutyArray = zeros(totalNumCases, 1);
outputVoltage = zeros(totalNumCases, 1);
Xarray = zeros(totalNumCases, 1);
final_voltage = [];
final_duty = [];

count = 1; % Initialize counter for storing data

display(['Total Number of cases: ', num2str(totalNumCases)]);

for grdval = 1:lenVoltage
    for loadR = 1:lenLoadVal
        Uout0 = Uout(grdval);  % Output voltage for this iteration
        RloadVal = Rload(loadR); % Load resistance for this iteration

        Vdc0 = 350; % Input voltage
        D0 = Uout0 / Vdc0; % Duty cycle
        Iload0 = Uout0 / RloadVal; % Load current
        iLf0 = Iload0; % Initial inductor current
        uRef0 = Uout0; % Reference voltage
        Imax = 2 * Iload0; % Maximum current
        
        P = bodeoptions;
        P.FreqUnits = 'Hz';
        P.PhaseWrapping = 'on';
        P.Grid = 'on';
        w = logspace(2, 5, 10000);
        s = tf('s');
        
        % Filter parameters
        fs = 48e3; % Switching frequency
        fsamp = 48e3; % Sampling frequency
        Tsamp = 1 / fsamp;
        Lf = 0.5 * 1 / (fs) * 1 / (0.2 * Iload0) * (1 - 0.5) * Vdc0;
        Cf = 1 / (Lf * (0.05 * 2 * pi * fs)^2);
        
        Cdc = 163e-6;
        rCdc = 0.2;
        
        Lg = 50e-6;
        rLg = 0;
        iLg0 = 0;
        
        f0grid = 1 / (2 * pi * sqrt(Cdc * Lg));
        
        % Total plant delay
        Td = 1.75 / fs;
        [num, den] = pade(Td, 2); % Approximation for dead-time
        Gt = tf(num, den);
        
        % Current controller design
        kuff = 1;
        Gi = (Cf * Gt * s) / (1 - Gt * kuff + Cf * Lf * s^2);
        
        wBWi = 2 * pi * fs * 0.07;
        Tni = Cf * Lf / Td;
        kpi = wBWi * Lf;
        Ri = kpi * (1 + s * Tni) / (s * Tni);
        
        Giol = Ri * Gi;
        Gicl = feedback(Giol, 1);
        
        % Voltage controller design
        Gu = 1 / (s * Cf) * Gicl;
        
        wBWu = 0.1 * wBWi;
        kpu = wBWu * Cf;
        Tnu = 10 * 1 / (wBWu);
        Ru = kpu * (1 + s * Tnu) / (s * Tnu);
        
        Guol = Ru * Gu;

        % Simulation of the buck converter
        disp(['** Start Simulation with Output Voltage: ', num2str(Uout0), ' and Load Resistance: ', num2str(RloadVal)]);
        sim('buck_dev.slx', 0.4); % Ensure the Simulink model name is correct
        % disp(['** Perturbation Simulation Started: ', num2str(Uout0), ' and Load Resistance: ', num2str(RloadVal)]);
        % sim('buck_dev_perturb.slx', 1);
        
        
        % Load the data from the .mat file
        IndCurr = load('buck_converter_Ilf.mat');
        % IndCurr_perturb = load('buck_converter_Ilf_per.mat');

        % Assuming llf is a time series object, access its Time and Data properties
        t = IndCurr.Ilf.Time;              % Extract the time values
        data = IndCurr.Ilf.Data;            % Extract the data matrix
        
        % Perturb
        t_perturb = IndCurr.Ilf.Time;              % Extract the time values
        data_perturb = IndCurr.Ilf.Data;            % Extract the data matrix
        current__perturb = data_perturb(:, 1); 

        % Extract specific columns based on the provided description
        current = data(:, 1);               % Extract the current column
        duty = data(:, 2);                  % Extract the duty column
        voltage = data(:, 3);               % Extract the voltage column

        final_voltage = [final_voltage; {voltage}]; % Append voltage row-wise
        final_duty = [final_duty; {duty}];          % Append duty row-wise
        
        % Save the updated arrays to a .mat file after each iteration
        save('final_data.mat', 'final_voltage', 'final_duty', '-v7.3');

        % Find indices within the specified time range
        timeRangeIndices = find(t >= tStart & t <= tEnd);

        % Calculate mean values over the specified time range
        IlfArray(count, 1) = mean(current(timeRangeIndices));   % Mean inductor current in the time range
        dutyArray(count, 1) = mean(duty(timeRangeIndices));     % Mean duty cycle in the time range
        outputVoltage(count, 1) = mean(voltage(timeRangeIndices));% Mean voltage in the time range
        Ilf_preArray(count, 1) = mean(current__perturb(timeRangeIndices));
        save('IlfFinal.mat', 'IlfArray', '-v7.3');

        % Calculate Xarray value
        Xarray(count, 1) = 5 / (Ilf_preArray(count, 1) - IlfArray(count, 1));

        disp(['Saved result. Completed ', num2str(count), ' runs...']);
        count = count + 1;
    end
end

save('IlfFinal.mat', 'IlfArray', '-v7.3');

% Create a table for training data
trainingData = table(outputVoltage, dutyArray, IlfArray, Ilf_preArray, Xarray);
disp('Data is saved in trainingData file.');

% Plot bode diagrams
% figure();
% bode(Gi, w, P);
% grid on;
% hold on;
% bode(Ri, w, P);
% bode(Giol, w, P);
% legend('G_{iL}', 'R_{i}', 'R_{i}*G_{iL}');
% title('Compensated Open-Loop Transfer Function of Converter Current');
% 
% figure();
% bode(Gu, w, P);
% grid on;
% hold on;
% bode(Ru, w, P);
% bode(Guol, w, P);
% legend('G_{uC}', 'R_{u}', 'R_{u}*G_{uC}');
% title('Compensated Open-Loop Transfer Function of Converter Output Voltage');
