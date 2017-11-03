function [ QoI,time ] = f( v )
%% This function is an adaptation of the script "nvu_run_script.m"
%% The input v is a vector which allows us to vary the parameters in the model.
%% The length of v corresponds to the number of parameters we vary.
%% If v(i) = 1 then parameter i is fixed to its nominal value.
%% If we want to vary parameters i with p% uncertainty about its nominal value
%% then we sample v(i) from a uniform distribution on the interval [1-p,1+p].

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

odeopts = odeset('RelTol', 1e-04, 'AbsTol', 1e-04, 'MaxStep', 0.5, 'Vectorized', 1);

XLIM2 = 150; % End of simulation

% For current type 1 use max current strength 0.022
% For current type 3 use max current strength 0.042
% For current type 4 use max current strength 0.036

CURRENT_STRENGTH    = 0.036;    % Max strength of current input in mA/cm2
NEURONAL_START      = 100;      % Start of neuronal stimulation
CURRENT_TYPE        = 4;        % Types of current input. 1: normal, 2: two stimulations (second stimulation is 8 sec after and 1 sec long), 3: obtained from experimental input data, 4: whisker pad (from experiment) + locus coeruleus (pain pathway)

% Used if CURRENT_STRENGTH = 1 or 2
NEURONAL_END        = 102;      % End of neuronal stimulation 

% Used if CURRENT_STRENGTH = 3 or 4
ISI = 7;                        % INDEX for time period between stimulations [0.6,1,2,3,4,6,8]
stim = 3;                       % INDEX for length of initial stimulation [2,8,16]

% Used if CURRENT_STRENGTH = 4: scaling for the two stimulation components,
% alpha for whisker pad and beta for locus coeruleus/pain. Default for both
% is 1. Can make parameters in vector v if desired.
% I_total = alpha * I_Wh + beta * I_LC
alpha = 1;
beta = 1;

% Not currently used
ECS_START       = 100000000;      % Start of ECS K+ input
ECS_END         = 1000000000;      % End of ECS K+ input

J_PLC           = 0.11;      % Jplc value in EC: 0.11 for steady state, 0.3 for oscillations
GLU_SWITCH      = 1;        % Turn on glutamate input (for NO and Ca2+ pathways)
NO_PROD_SWITCH  = 1;        % Turn on Nitric Oxide production 
TRPV_SWITCH     = 1;        % Turn on TRPV4 Ca2+ channel from AC to PVS
O2SWITCH        = 1;        % 0: ATP is plentiful, 1: ATP is limited (oxygen-limited regime, default)

% Partition v into its part in each component
% v(1) through v(99) in Neuron
% v(100) through v(104) in Astrocyte
% v(105) throught v(132) in SMCEC

% The Neuron channels parameters are in v(1) throught v(70)
% The buffer parameters are in v(71) throught v(74)
% The NO pathway parameters are in v(75) through v(132)

% I give a subset of v as an input argument to each component class,
% Neuron, Astrocyte, and SMCEC. By doing this they begin indexing their
% own parameters at 1 within each class. This allows us to change
% things in one class without having to touch the others.

% Load initial NVU
nv = NVU(Neuron(v(1:99),'V_maxNOS', 1*25e-3, 'SC_coup', 11.5, 'CurrentType', CURRENT_TYPE, 'O2switch', O2SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 'Istrength', CURRENT_STRENGTH, 'GluSwitch', GLU_SWITCH, 'NOswitch', NO_PROD_SWITCH, 't0_ECS', ECS_START, 'ECS_input', 9), ...
    Astrocyte(v(100:104),'R_decay', 0.15, 'trpv_switch', TRPV_SWITCH, 'startpulse', NEURONAL_START, 'lengthpulse', NEURONAL_END - NEURONAL_START, 't0_ECS', ECS_START, 'tend_ECS', ECS_END, 'Rk_switch', 0), ...
    WallMechanics('wallMech', 1.7), ...
    SMCEC(v(105:132),'J_PLC', J_PLC, 'NOswitch', NO_PROD_SWITCH), 'odeopts', odeopts);

% Adjust time vector
nv.neuron.params.dt = 0.1; dt = nv.neuron.params.dt;
nv.T = 0:dt:XLIM2;
% numTimeSteps = length(nv.T);

%% Load whisker stimulation input data from file, save as I_Wh and put into Neuron.m
if nv.neuron.params.CurrentType == 3 || nv.neuron.params.CurrentType == 4
    
    load neurovascular_data_for_tim_david.mat
    actual_ISI = info.isi_duration(ISI);
    actual_stim = info.condition_stim_duration(stim);
    sum_neural_wh = zeros(size(neural_tim_vector));
    for animal = 1:11
        for experiment = 1:10
            sum_neural_wh = sum_neural_wh+neural_data(:,ISI,stim,experiment,animal)'; % Sum all data 
        end
    end
    mean_neural_wh = sum_neural_wh./110;  % Average the neural data over all animals and experiments, animals*experiments=110
    neural_tim_vector_shifted = neural_tim_vector + NEURONAL_START;    % Shift so stimulation begins at NEURONAL_START
    interp_neural_wh = interp1(neural_tim_vector_shifted, mean_neural_wh, nv.T); % Interpolate so there is data for all timesteps for NVU
    interp_neural_wh(isnan(interp_neural_wh))=0.02;   % Remove NaNs     
    
    nv.neuron.input_data = interp_neural_wh;   % Replace dummy in Neuron with input data, leave it at that for CURRENT_TYPE = 3
    I_Wh = interp_neural_wh;                   % Save as whisker pad current I_Wh
end

%% Construct additional pain pathway stimulation I_LC and input total current I_total to Neuron.m
if nv.neuron.params.CurrentType == 4
    
    % Construct first stimulation, set to work for any NEURONAL_START or any
    % initial stimulus duration (2,8,16 s)
    time_1 = linspace(0, actual_stim, 10000); 
    I_1 = 0.0006 * time_1.^2 + 0.1;             % Chosen for the shape! Can be modified if needed
    time_1 = time_1 + NEURONAL_START;           % Shift so starts at NEURONAL_START
    I_1 = interp1(time_1, I_1, nv.T);           % Make same size as nv.T
    I_1(isnan(I_1)) = 0;                        % Replace NaNs with zeros
    
    % Construct second stimulation with duration 2 sec
    time_2 = linspace(0, 1, 10000);
    I_2 = 0.1*ones(size(time_2));
    time_2 = time_2 + NEURONAL_START + actual_stim + actual_ISI;
    I_2 = interp1(time_2, I_2, nv.T);
    I_2(isnan(I_2)) = 0;
    
    % Add together
    I_LC = I_1 + I_2;
    
    % Total current (whisker pad plus LC)
    I_total = alpha*I_Wh + beta*I_LC;
    
    % Input to Neuron.m
    nv.neuron.input_data = I_total;
end

%%
% figure;
% plot(nv.T, alpha*I_Wh, nv.T, beta*I_LC, nv.T, I_total);
% legend('I_Wh','I_LC','I_{total}')
% xlim([95 140])

%% Run the simulation
nv.simulate() 

timeEnd = datetime('now');
fprintf('End time is %s\n', char(timeEnd));

time = nv.T;
QoI = nv.U;


end

