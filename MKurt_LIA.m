
% This program demonstrates the MKurt-LIA method for different measurement tasks.
% Some pre-set parameters, such as the theoretical fault characteristic frequency, are set according to the outer race failure condition of bearings 1 and 2.
% The mkurt spectrum method comes from citation 26 of the article.
% Bearing data are from the article citation 30.
% Meng's mail: meng1.zhang@mail.polimi.it

%% selection of bearing
folder_path = 'Bearing1'; 
files = dir(fullfile(folder_path, '*.csv'));
num_files = length(files);
direc=1; %use horizon sensor

%% single measurement on 1-min data

data = csvread('Bearing1/1.csv', 1, 0); % assigning time segment

signal = data(:, direc);

data_length = length(signal);

time = (0:data_length-1) / 25600;

% assigning reference frequency
frequency = 107.9074; % Hz

% reference signals
sin_wave = sin(2 * pi * frequency * time');
cos_wave = cos(2 * pi * frequency * time');

% PSD
result_signal_sin = signal .* sin_wave;
result_signal_cos = signal .* cos_wave;

% LPF
integration_average_sin = sum(result_signal_sin) / length(result_signal_sin);
integration_average_cos = sum(result_signal_cos) / length(result_signal_cos);

% calculate amplitude on reference frequency
amp = 2*(sqrt(integration_average_sin^2+integration_average_cos^2));


%% check the original signal

csv_files = dir(fullfile(folder_path, '*.csv'));

file_names = {csv_files.name};
file_numbers = cellfun(@(x) sscanf(x, '%d.csv'), file_names);

[~, sorted_indices] = sort(file_numbers);

sorted_csv_files = csv_files(sorted_indices);

combined_data = [];

for i = 1:length(sorted_csv_files)
    file_path = fullfile(folder_path, sorted_csv_files(i).name);
    
    data = csvread(file_path, 1, 0);

    combined_data = [combined_data; data];
end

total_length = size(combined_data, 1);

time = (0:total_length-1) / 60/(32768/60); 
figure;

plot(time, combined_data(:, 1), 'LineWidth', 1.1, 'Color', 'b');set(gca, 'LineWidth', 1.1); 
xlabel('time [min]');
ylabel('amplitude (g)');
xlim([0 num_files]);
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');


%% long-term measurement using theoretical fault characteristic frequency as reference frequency

csv_files = dir(fullfile(folder_path, '*.csv'));
num_files = length(csv_files);

file_list = dir(fullfile(folder_path, '*.csv'));

file_names = {file_list.name};
file_numbers = cellfun(@(x) str2double(x(1:end-4)), file_names);

[~, sorted_idx] = sort(file_numbers);

sorted_file_list = file_list(sorted_idx);

amp_values_unadj = zeros(1, num_files);

for i = 1:length(sorted_file_list)
    file_path = fullfile(folder_path, sorted_file_list(i).name);

    data = csvread(file_path, 1, 0);
    signal = data(:, direc);

    data_length = length(signal);
    time = (0:data_length-1) / 25600; 

    frequency = 107.9074; % using theoretical fault characteristic frequency as reference frequency

    sin_wave = sin(2 * pi * frequency * time');
    cos_wave = cos(2 * pi * frequency * time');

    result_signal_sin = signal .* sin_wave;
    result_signal_cos = signal .* cos_wave;

    integration_average_sin = sum(result_signal_sin) / length(result_signal_sin);
    integration_average_cos = sum(result_signal_cos) / length(result_signal_cos);

    amp = 2 * (sqrt(integration_average_sin^2 + integration_average_cos^2));

    amp_values_unadj(i) = amp;
    disp(file_path); 
end

%% long-term measurement using MKurt spectrum frequency as reference frequency

csv_files = dir(fullfile(folder_path, '*.csv'));
num_files = length(csv_files);

file_list = dir(fullfile(folder_path, '*.csv'));

file_names = {file_list.name};
file_numbers = cellfun(@(x) str2double(x(1:end-4)), file_names);

[~, sorted_idx] = sort(file_numbers);

sorted_file_list = file_list(sorted_idx);

amp_values = zeros(1, num_files);

all_frequencies = zeros(1, num_files);

for i = 1:length(sorted_file_list)

    file_path = fullfile(folder_path, sorted_file_list(i).name);

    data = csvread(file_path, 1, 0);
    signal = data(:, direc);

    data_length = length(signal);
    time = (0:data_length-1) / 25600;
   
    % Mkurt for calculating fault frequency
    window = ones(1,1); 
    range = [25600/112:0.1:25600/104];
    L=1000;
    [T MKurt f y T_best MKurt_best f_best y_best t] = momeda_spectrum(signal,L,window,range,0);
    
    % find the peak
    [max_val, max_idx] = max(MKurt);
    value_at_max_idx = range(max_idx);
    frequency=25600/value_at_max_idx;
    all_frequencies(i) = frequency;

    sin_wave = sin(2 * pi * frequency * time');
    cos_wave = cos(2 * pi * frequency * time');

    result_signal_sin = signal .* sin_wave;
    result_signal_cos = signal .* cos_wave;
    
    integration_average_sin = sum(result_signal_sin) / length(result_signal_sin);
    integration_average_cos = sum(result_signal_cos) / length(result_signal_cos);
    
    amp = 2 * (sqrt(integration_average_sin^2 + integration_average_cos^2));
    
    amp_values(i) = amp;

    disp(file_path);
end

figure;
plot(all_frequencies, 'LineWidth', 1.1, 'color', 'b');set(gca, 'LineWidth', 1.1); 
xlabel('time [min]');
ylabel('frequency [Hz]');
xlim([0 num_files]);
ylim([106 110])
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

figure;
plot(amp_values, 'LineWidth', 1.1, 'color', 'b');set(gca, 'LineWidth', 1.1); 
hold on; 
plot(amp_values_unadj, 'LineWidth', 1.1, 'Color', 'r');
xlabel('time [min]');
ylabel('amplitude (g)');
xlim([0 num_files]);
legend('proposed method', 'original method', 'Location', 'NorthWest');
legend('Color', 'none');
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');



