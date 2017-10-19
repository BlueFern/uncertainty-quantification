%% Script to extract and plot averaged time series for CBF data

clear all

% Load data from file - should be in same folder
load neurovascular_data_for_tim_david.mat

% Different types of stimuli
ISI_vector = [0.6,1,2,3,4,6,8];
stimulus_duration_vector = [2,8,16];
time = cbf_tim_vector;

% Put averaged time series into a vector, to extract:
% Averaged_Stimulations_CBF(:, ISI_index, stimulus_duration_index)

Averaged_Stimulations_CBF = zeros(length(time),7,3);

for stimulus_duration_index = 1:length(stimulus_duration_vector)
    figure(stimulus_duration_index);
    hold on
    
    % Obtain stimulus duration
    stimulus_duration = stimulus_duration_vector(stimulus_duration_index)
    
    for ISI_index = 1:length(ISI_vector)

        % Get ISI
        ISI = ISI_vector(ISI_index);

        % Sum all trials for all animals
        sum_cbf = zeros(size(cbf_tim_vector));
        for animal = 1:11
            for experiment = 1:10
                sum_cbf = sum_cbf+cbf_data(:,ISI_index,stimulus_duration_index,experiment,animal)';
            end
        end
        average_cbf = sum_cbf./110;  % Obtain average CBF
        delta_cbf = average_cbf - 1; % To obtain change in CBF

        %Add to plot
        plot(time, delta_cbf);
        
        %Add to vector
        Averaged_Stimulations_CBF(:,ISI_index,stimulus_duration_index) = delta_cbf;
        
    end
    hold off

    % Label plot
    xlabel('Time [s]');
    ylabel('\Delta CBF')
    title(['CBF with initial duration ' num2str(stimulus_duration) ' sec']);
    legend('0.6','1','2','3','4','6','8')

end