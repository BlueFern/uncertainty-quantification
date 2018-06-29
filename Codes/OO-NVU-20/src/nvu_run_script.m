%% Demonstration script for new NVU model

%% Version control
% This is the latest version of the NVU model: Version 2.0
% K+, NO, Astrocytic Ca2+, TRPV4, ECS, Neuron

%% Construct NVU
% The NVU consists of a number of submodules, implemented as MATLAB
% classes, presently an astrocyte, a lumped SMC/EC model, and a model of
% the wall mechanics. 
%
% The parameters of each of these submodules are
% specified when the modules are constructed, here in the call to NVU, see
% for example the |SMCEC| part of the |NVU| call below:
%
% Options for the ODE solver (currently |ode15s|) are provided by
% specifying the |odeopts| parameter. The code works fine with default
% tolerances.

%clear all

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM1 = 90; XLIM2 = 250;
FIG_NUM = 1;

CURRENT_STRENGTH    = 0.042;  % Max strength of current input in mA/cm2
CURRENT_TYPE        = 3;      % Types of current input. 1: normal, 2: two stimulations, 3: obtained from input data
ISI = 7;    % INDEX for time period between stimulations [0.6,1,2,3,4,6,8]
stim = 3;   % INDEX for length of initial stimulation [2,8,16]
NEURONAL_START      = 100;      % Start of neuronal stimulation

% Used if CURRENT_STRENGTH = 1 or 2
NEURONAL_END        = 102;      % End of neuronal stimulation 

% Not currently used
ECS_START       = 100000000;      % Start of ECS K+ input
ECS_END         = 1000000000;      % End of ECS K+ input

J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations

GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
O2SWITCH        = 1;        % 0: ATP is plentiful, 1: ATP is limited (oxygen-limited regime, default)

% % Output start time and estimated finish time based on how many sec of
% % stimulation there will be and speed of computer (need to check accuracy
% % when loading data)
% %(home: 0.38, work: 0.88)
% estimatedMinutes = (NEURONAL_END - NEURONAL_START)*0.88;
% timeStart = datetime('now');
% estimatedTimeEnd = timeStart + minutes(estimatedMinutes); 
% fprintf('Start time is %s\nEstimated end time is %s\n', char(timeStart), char(estimatedTimeEnd));

% Load initial NVU
nv = NVU(Neuron('V_maxNOS', 1*25e-3, 'SC_coup', 11.5, 'CurrentType', CURRENT_TYPE, 'O2switch', O2SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, 't0_ECS', ECS_START, 'ECS_input', 9), ...
    Astrocyte('R_decay', 0.15, 'trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END, 'Rk_switch', 0), ...
    WallMechanics('wallMech', 1.7), ...
    SMCEC('J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);

% Adjust time vector
nv.neuron.params.dt = 0.001; dt = nv.neuron.params.dt;
nv.T = 0:dt:XLIM2;
numTimeSteps = length(nv.T);

% Load input data from file, save to input_data and put into Neuron
load neurovascular_data_for_tim_david.mat
actual_ISI = info.isi_duration(ISI);
actual_stim = info.condition_stim_duration(stim);
sum_neural = zeros(size(neural_tim_vector));
for animal = 1:11
    for experiment = 1:10
        sum_neural = sum_neural+neural_data(:,ISI,stim,experiment,animal)'; % Sum all data 
    end
end
mean_neural = sum_neural./110;  % Average the neural data over all animals and experiments, animals*experiments=110
neural_tim_vector_shifted = neural_tim_vector + NEURONAL_START;    % Shift so stimulation begins at NEURONAL_START
interp_neural = interp1(neural_tim_vector_shifted, mean_neural, nv.T); % Interpolate so there is data for all timesteps for NVU
interp_neural(isnan(interp_neural))=0.02;   % Remove NaNs     
nv.neuron.input_data = interp_neural;   % Replace dummy in Neuron with input data

%% Run the simulation
nv.simulate() 

% % Run this to take the ICs as the end of the last simulation run i.e. steady state ICs
% ICs = (nv.U(end, :))';
% nv.u0 = ICs;
% nv.simulateManualICs() 

% Find a point in time 20 sec before neuronal stimulation has begun (preNeuronalStimTime1)
% and the point in time when stimulation begins (preNeuronalStimTime2) and get the values for 
% CBF, CBV, HBR, CMRO2 at that time, then find the midpoint between the min and max and normalise with that value. 
% Done in this way so that it also works when the variables are oscillatory (i.e. J_PLC = 0.3)
np = nv.neuron.params; % Shortcut for neuron parameters
preNeuronalStimTime1 = floor((NEURONAL_START-20)*numTimeSteps/XLIM2);
preNeuronalStimTime2 = floor((NEURONAL_START)*numTimeSteps/XLIM2);
CBF = nv.out('CBF'); CBV = (nv.out('CBV'))'; HBR = (nv.out('HBR'))'; CMRO2 = nv.out('CMRO2');
CBF_0 = 0.5*( max(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBF(preNeuronalStimTime1:preNeuronalStimTime2)) );
CBV_0 = 0.5*( max(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CBV(preNeuronalStimTime1:preNeuronalStimTime2)) );
HBR_0 = 0.5*( max(HBR(preNeuronalStimTime1:preNeuronalStimTime2)) + min(HBR(preNeuronalStimTime1:preNeuronalStimTime2)) );
CMRO2_0 = 0.5*( max(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) + min(CMRO2(preNeuronalStimTime1:preNeuronalStimTime2)) );

% Normalised variables
CBF_N = CBF./CBF_0;
CBV_N = CBV./CBV_0;
HBR_N = HBR./HBR_0;
CMRO2_N = CMRO2./CMRO2_0;

HBT_N = CBF_N .* HBR_N ./ CMRO2_N;                                          % Total hemoglobin (normalised)
HBO_N = (HBT_N - 1) - (HBR_N - 1) + 1;                                      % Oxyhemoglobin (normalised)
BOLD_N = 100 * np.V_0 * ( np.a_1 * (1 - HBR_N) - np.a_2 * (1 - CBV_N) );    % BOLD (percentage increase from 0)

%% Plot experimental and model CBF from data file
sum_cbf = zeros(size(cbf_tim_vector));
for animal = 1:11
    for experiment = 1:10
        sum_cbf = sum_cbf+cbf_data(:,ISI,stim,experiment,animal)';
    end
end
mean_cbf = (sum_cbf./110) - 1;
cbf_tim_vector_shifted = cbf_tim_vector + NEURONAL_START;    % Shift so stimulation begins at NEURONAL_START
figure;
plot(cbf_tim_vector_shifted, mean_cbf, ':', nv.T, (nv.out('CBF')-CBF_0)./CBF_0, 'LineWidth', 1);
ylabel('\Delta CBF')
xlabel('Time [s]')
xlim([90 150])
ylim([-0.05 0.3])
title(['CBF with initial duration ' num2str(actual_stim) ', ISI ' num2str(actual_ISI)] );
p1=patch([100 100+actual_stim 100+actual_stim 100],[-0.05 -0.05 0.3 0.3],'k');
set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[-0.05 -0.05 0.3 0.3],'k');
set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
% 

%% Plot CBF
% figure(FIG_NUM);
% hold all;
% plot(nv.T, nv.out('CBF')./CBF_0, 'k', 'LineWidth', 1);
% ylabel('CBF [-]');
% xlim([XLIM1 XLIM2])
% ylim([0.97 1.4])
% title('Normalised CBF');
% xlabel('time (s)');
% p1=patch([100 116 116 100],[0.97 0.97 1.4 1.4],'k');
% set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
% p2=patch([124 125 125 124],[0.97 0.97 1.4 1.4],'k');
% set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
% 
% %% Plot current
% figure(200);
% hold all;
% plot(nv.T, nv.out('current'));
% ylabel('current');

%% Plot hemoglobin
% 
% figure(FIG_NUM);
% plot(nv.T, HBO_N-1, 'r', nv.T, HBR_N-1, 'b', nv.T, HBT_N-1, 'g', 'LineWidth', 2);
% p=patch([100 104 104 100],[-0.15 -0.15 0.25 0.25],'k');
% set(p,'FaceAlpha',0.1,'EdgeColor', 'none');
% xlim([95 130]);
% ylabel('\Delta Concentration');
% xlabel('Time [s]');
% legend('HbO','HbR','HbT');
% 
% figure(FIG_NUM+1);
% plot(nv.T, BOLD_N, 'k', 'LineWidth', 2);
% p=patch([100 104 104 100],[-0.6 -0.6 1.7 1.7],'k');
% set(p,'FaceAlpha',0.1,'EdgeColor', 'none');
% xlim([95 130]);
% ylim([-0.6 1.7]);
% ylabel('\Delta BOLD');
% xlabel('Time [s]');

% figure(FIG_NUM+2);
% plot( nv.T, nv.out('J_KDR_sa'),  nv.T, nv.out('J_KA_sa'), nv.T, nv.out('J_Kleak_sa'), nv.T, nv.out('J_Kpump_sa'), nv.T, nv.out('J_KDR_d'), nv.T, nv.out('J_KA_d'), nv.T, nv.out('J_Kleak_d'), nv.T, nv.out('J_Kpump_d'), nv.T, nv.out('J_NMDA_K_d'));
% xlim([XLIM1 XLIM2]);
% legend('KDR_s','KA_s','Kleak_s','Kpump_s','KDR_d','KA_d','Kleak_d','Kpump_d','NMDA_d')

figure(FIG_NUM+1);
subplot(2,2,1);
    hold all;
    plot(nv.T, nv.out('nNOS_act_n'), 'LineWidth', 1);
    ylabel('nNOS');
    xlim([XLIM1 XLIM2])
    xlabel('Time [s]')
    p1=patch([100 100+actual_stim 100+actual_stim 100],[0.3 0.3 0.9 0.9],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[0.3 0.3 0.9 0.9],'k');
    set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
subplot(2,2,2);
    hold all;
    plot(nv.T, nv.out('NO_n'), 'LineWidth', 1);
    ylabel('NO_n');
    xlim([XLIM1 XLIM2])
    xlabel('Time [s]')
    p1=patch([100 100+actual_stim 100+actual_stim 100],[0.15 0.15 0.45 0.45],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[0.15 0.15 0.45 0.45],'k');
    set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
subplot(2,2,3);
    hold all;
    plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
    ylabel('Radius');
    xlim([XLIM1 XLIM2])
    xlabel('Time [s]')
    p1=patch([100 100+actual_stim 100+actual_stim 100],[22.5 22.5 25.5 25.5],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[22.5 22.5 25.5 25.5],'k');
    set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
subplot(2,2,4);
    hold all;
    plot(nv.T, CBF_N, 'LineWidth', 1);
    ylabel('CBF');
    xlim([XLIM1 XLIM2])
    xlabel('Time [s]')
    p1=patch([100 100+actual_stim 100+actual_stim 100],[1 1 1.4 1.4],'k');
    set(p1,'FaceAlpha',0.1,'EdgeColor', 'none');
    p2=patch([100+actual_stim+actual_ISI 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI+1 100+actual_stim+actual_ISI],[1 1 1.4 1.4],'k');
    set(p2,'FaceAlpha',0.1,'EdgeColor', 'none');
    
% Plot variables
% 
figure(FIG_NUM+1000);
subplot(3,3,1);
    hold all;
    plot(nv.T, nv.out('v_sa'), 'LineWidth', 1);
    ylabel('v_{sa} [mV]');
    xlim([XLIM1 XLIM2])
subplot(3,3,2);
    hold all;
    plot(nv.T, nv.out('current'), 'LineWidth', 1);
    ylabel('current');
    xlim([XLIM1 XLIM2])
subplot(3,3,3);
    hold all;
    plot(nv.T, nv.out('O2')*1e3, 'LineWidth', 1);
    ylabel('O2 [\muM]');
    xlim([XLIM1 XLIM2])
subplot(3,3,4);
    hold all;
    plot(nv.T, (nv.out('CBF')-CBF_0)./CBF_0, 'LineWidth', 1);
    ylabel('\Delta CBF / CBF_0');
    xlim([XLIM1 XLIM2])
subplot(3,3,5);
    hold all;
    plot(nv.T, nv.out('R')*1e6, 'LineWidth', 1);
    ylabel('Radius [\mum]');
    xlim([XLIM1 XLIM2])
subplot(3,3,6);
    hold all;
    plot(nv.T, nv.out('K_s')/1e3, 'LineWidth', 1);
    ylabel('K_s');
    xlim([XLIM1 XLIM2])
subplot(3,3,7);
    hold all;
    plot(nv.T, CBV_N, 'LineWidth', 1);
    ylabel('CBV');
    xlim([XLIM1 XLIM2])
subplot(3,3,8);
    hold all;
    plot(nv.T, HBR_N, 'LineWidth', 1);
    ylabel('HBR');
    xlim([XLIM1 XLIM2])
subplot(3,3,9);
    hold all;
    % BOLD using normalised variables
    plot(nv.T, BOLD_N, 'LineWidth', 1);  
    ylabel('\Delta BOLD (%)');
    xlim([XLIM1 XLIM2])
% %  % Plot figures - whatever you want

%% Plot for paper
% figure(10);
% subplot(4,2,1);
%     hold all;
%     plot(nv.T, nv.out('v_sa'), 'k', 'LineWidth', 1);
%     ylabel('v_{sa} [mV]');
%     xlim([XLIM1 XLIM2])
%     ylim([-100 20])
%     title('A. Soma membrane potential');
%     xlabel('time (s)');
%     rectangle('Position',[300, -100, 16, 5],'FaceColor',[0 0 0])
%     rectangle('Position',[320, -100, 2, 5],'FaceColor',[0 0 0])
% subplot(4,2,2);
%     hold all;
%     plot(nv.T, nv.out('K_e'), 'k', 'LineWidth', 1);
%     ylabel('K_e [mM]');
%     xlim([XLIM1 XLIM2])
%     title('B. Extracellular potassium');
%     xlabel('time (s)');
%     rectangle('Position',[300, 2, 16, 0.2],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 2, 2, 0.2],'FaceColor',[0 0 0])
% subplot(4,2,3);
%     hold all;
%     plot(nv.T, nv.out('R')*1e6, 'k', 'LineWidth', 1);
%     ylabel('Radius [\mum]');
%     ylim([22.4 25])
%     xlim([XLIM1 XLIM2])
%     title('C. Radius');
%     xlabel('time (s)');
%     rectangle('Position',[300, 22.4, 16, 0.1],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 22.4, 2, 0.1],'FaceColor',[0 0 0])
% subplot(4,2,4);
%     hold all;
%     plot(nv.T, (nv.out('CBF')-CBF_0)./CBF_0, 'k', 'LineWidth', 1);
%     ylabel('\Delta CBF / CBF_0 [-]');
%     xlim([XLIM1 XLIM2])
%     ylim([-0.1 0.4])
%     title('D. Normalised \DeltaCBF');
%     xlabel('time (s)');
%     rectangle('Position',[300, -0.1, 16, 0.02],'FaceColor',[0 0 0])
%     rectangle('Position',[320, -0.1, 2, 0.02],'FaceColor',[0 0 0])
% subplot(4,2,5);
%     hold all;
%     plot(nv.T, nv.out('O2'), 'k', 'LineWidth', 1);
%     ylabel('O_2 [mM]');
%     xlim([XLIM1 XLIM2])
%     ylim([0.0265 0.03])
%     title('E. Oxygen concentration');
%     xlabel('time (s)');
%     rectangle('Position',[300, 0.0265, 16, 0.00015],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 0.0265, 2, 0.00015],'FaceColor',[0 0 0])
% subplot(4,2,6);
%     hold all;
%     plot(nv.T, nv.out('CBV')./CBV_0, 'k', 'LineWidth', 1);
%     ylabel('CBV [-]');
%     xlim([XLIM1 XLIM2])
%     ylim([0.97 1.1])
%     title('F. Normalised CBV');
%     xlabel('time (s)');
%     rectangle('Position',[300, 0.97, 16, 0.005],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 0.97, 2, 0.005],'FaceColor',[0 0 0])
% subplot(4,2,7);
%     hold all;
%     plot(nv.T, nv.out('HBR')./HBR_0, 'k', 'LineWidth', 1);
%     ylabel('HBR [-]');
%     xlim([XLIM1 XLIM2])
%     ylim([0.86 1.05])
%     title('G. Normalised deoxyhemoglobin');
%     xlabel('time (s)');
%     rectangle('Position',[300, 0.86, 16, 0.008],'FaceColor',[0 0 0])
%     rectangle('Position',[320, 0.86, 2, 0.008],'FaceColor',[0 0 0])
% subplot(4,2,8);
%     hold all;
%     plot(nv.T, 100 * np.V_0 * ( np.a_1 * (1 - nv.out('HBR')/HBR_0) - np.a_2 * (1 - nv.out('CBV')/CBV_0) ), 'k', 'LineWidth', 1);
%     ylabel('\Delta BOLD (%)');
%     xlim([XLIM1 XLIM2])
%     ylim([-0.9 1.5])
%     title('H. BOLD response');
%     xlabel('time (s)');
%     rectangle('Position',[300, -0.9, 16, 0.09],'FaceColor',[0 0 0])
%     rectangle('Position',[320, -0.9, 2, 0.09],'FaceColor',[0 0 0])

timeEnd = datetime('now');
fprintf('End time is %s\n', char(timeEnd));
