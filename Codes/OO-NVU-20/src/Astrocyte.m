classdef Astrocyte < handle
    % The 'Astrocyte' code contains the following sections of the model:
    % The Neuron, Synaptic cleft, the Astrocyte and the Perivascular Space
    % Currently there is no content under the Neuron sub-section
    % Please refer to the relevient sections in the documentation for
    % full information on the equations and variable names.
    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = Astrocyte(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices(self);
            self.u0 = initial_conditions(self.index,self);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices(self);
        end

        function [du, varargout] = rhs(self, t, u, J_KIR_i, R, J_VOCC_i, NO_n, NO_i, J_K_NEtoSC, Glu)
            % Initalise inputs and parameters
            t = t(:).';
            p = self.params;
            idx = self.index;

            R_k = u(idx.R_k, :);
            K_p = u(idx.K_p, :);
            Ca_p = u(idx.Ca_p, :);
            N_Na_k = u(idx.N_Na_k, :);
            N_K_k = u(idx.N_K_k, :);
            N_Cl_k = u(idx.N_Cl_k, :);
            N_HCO3_k = u(idx.N_HCO3_k, :);
            N_Na_s = u(idx.N_Na_s, :);
            N_K_s = u(idx.N_K_s, :);
            N_HCO3_s = u(idx.N_HCO3_s, :);
            w_k = u(idx.w_k, :);
            I_k = u(idx.I_k, :);
            
            Ca_k = u(idx.Ca_k, :);
            h_k = u(idx.h_k, :);
            s_k = u(idx.s_k, :);
            m_k = u(idx.m_k, :); 
            eet_k = u(idx.eet_k, :);
            NO_k = u(idx.NO_k, :);
            du = zeros(size(u));
            
            %% Synaptic Cleft:

            % Electroneutrality condition
            N_Cl_s = N_Na_s + N_K_s - N_HCO3_s;
            
            % Volume-surface ratio 
            %TODO: set as constant and remove to simplify equations
            R_s = p.R_tot - R_k;
            
            % Scale concentrations to get actual concentrations in uM!!
            K_s = N_K_s ./ R_s;
            Na_s = N_Na_s ./ R_s;
            Cl_s = N_Cl_s ./ R_s;
            HCO3_s = N_HCO3_s ./ R_s;
            Na_k = N_Na_k ./ R_k;
            K_k = N_K_k ./ R_k;
            Cl_k = N_Cl_k ./ R_k;
            HCO3_k = N_HCO3_k ./ R_k;
            
            % Input of K+ to the SC (assuming that the SC is a small part of the ECS and everything that happens to the ECS also happens to the SC)
            J_K_NEtoSC_k = J_K_NEtoSC * 1000 .* R_s; % Convert from mM/s to uMm/s
 
            
            %% Astrocyte
            % Volume-surface Ratio scaling ODE
            du(idx.R_k, :) = p.Rk_switch * p.L_p * (Na_k + K_k + Cl_k + HCO3_k - Na_s - Cl_s - K_s - HCO3_s + p.X_k ./ R_k);
            
            % Nernst potentials ( in V)
            E_K_k = p.R_g * p.T / (p.z_K * p.F) * log(K_s ./ K_k);
            E_Na_k = p.R_g * p.T / (p.z_Na * p.F) * log(Na_s ./ Na_k);
            E_Cl_k = p.R_g * p.T / (p.z_Cl * p.F) * log(Cl_s ./ Cl_k);
            E_NBC_k = p.R_g * p.T / (p.z_NBC * p.F) * log((Na_s .* HCO3_s.^2) ./ (Na_k .* HCO3_k.^2));
            E_BK_k = p.reverseBK + p.switchBK *(p.R_g * p.T / (p.z_K * p.F) * log(K_p ./ K_k)); % nerst potential BK, either constant or as a function of K_k and K_p. [V]
            E_TRPV_k = p.R_g * p.T / (p.z_Ca * p.F) * log(Ca_p./Ca_k); % Nernst potential TRPV
            
            % Flux through the Sodium Potassium pump
            J_NaK_k = p.J_NaK_max * Na_k.^1.5 ./ (Na_k.^1.5 + p.K_Na_k^1.5) .* K_s ./ (K_s + p.K_K_s);
            
            % Membrane voltage
            g_BK_k = p.G_BK_k*1e-12 / p.A_ef_k;
            g_TRPV_k = (p.G_TRPV_k* 1e-12)/(p.A_ef_k);%mho m^-2
            v_k = (p.g_Na_k * E_Na_k + p.g_K_k * E_K_k + g_TRPV_k * m_k .* E_TRPV_k + p.g_Cl_k * E_Cl_k + p.g_NBC_k * E_NBC_k + g_BK_k * w_k .* E_BK_k - J_NaK_k * p.F / p.C_correction) ./ ...
                (p.g_Na_k + p.g_K_k + p.g_Cl_k + p.g_NBC_k + g_TRPV_k * m_k + g_BK_k * w_k);
            
            % Fluxes
            J_N_BK_k = g_BK_k / p.F * w_k .* (v_k - E_BK_k) * p.C_correction; % scaled BK flux (uM m /s)
            J_BK_p = J_N_BK_k ./ (R_k * p.VR_pa); % K+ influx into the PVS (uM/s)
            J_BK_k = J_N_BK_k ./ R_k; % K+ efflux from the AC (uM/s)
            J_K_k = p.g_K_k / p.F * (v_k - E_K_k) * p.C_correction;
            J_Na_k = p.g_Na_k / p.F * (v_k - E_Na_k) * p.C_correction;
            J_NBC_k = p.g_NBC_k / p.F * (v_k - E_NBC_k) * p.C_correction;
            J_KCC1_k = p.g_KCC1_k / p.F * p.R_g * p.T / p.F .* log((K_s .* Cl_s) ./ (K_k .* Cl_k)) * p.C_correction;
            J_NKCC1_k = p.g_NKCC1_k / p.F * p.R_g * p.T / p.F .* log((Na_s .* K_s .* Cl_s.^2) ./ (Na_k .* K_k .* Cl_k.^2)) * p.C_correction;
            
            %% Calcium Equations
            % Flux
            J_IP3 = p.J_max * ( I_k ./ (I_k + p.K_I) .*  Ca_k ./ (Ca_k + p.K_act) .* h_k).^3 .* (1 - Ca_k ./ s_k);
            J_ER_leak = p.P_L * (1 - Ca_k ./ s_k);
            J_pump = p.V_max * Ca_k.^2 ./ (Ca_k.^2 + p.k_pump^2);
            I_TRPV_k = p.G_TRPV_k * m_k .* (v_k-E_TRPV_k) * p.C_correction; % current TRPV
            J_TRPV_k = -0.5 * I_TRPV_k / (p.C_astr_k * p.gamma_k); 

            rho = p.rho_min + (p.rho_max - p.rho_min)/p.Glu_max * Glu;
            
            % Other equations
            B_cyt = 1 ./ (1 + p.BK_end + p.K_ex * p.B_ex ./ (p.K_ex + Ca_k).^2);
            G = (rho + p.delta) ./ (p.K_G + rho + p.delta);
            v_3 = p.v_6 - p.v_5 / 2 * tanh((Ca_k - p.Ca_3) / p.Ca_4);
            
            %% Parent Calcium equations
            
            w_inf = 0.5 * (1 + tanh((v_k + p.eet_shift * eet_k - v_3) / p.v_4));
            phi_w = p.psi_w * cosh((v_k - v_3) / (2 * p.v_4));
		   
            %% TRPV Channel open probability equations
            H_Ca_k = Ca_k ./ p.gam_cai_k + Ca_p ./ p.gam_cae_k;
            eta = (R - p.R_0_passive_k) ./ (p.R_0_passive_k);
            minf_k = (1 ./ (1 + exp(-(eta - p.epshalf_k) ./ p.kappa_k))) .* ((1 ./ (1 + H_Ca_k)) .* (H_Ca_k + tanh((v_k - p.v1_TRPV_k) ./ p.v2_TRPV_k))); 
            J_VOCC_k = J_VOCC_i;
            
            % NO pathway
            tau_nk = p.x_nk ^ 2 ./  (2 * p.D_cNO);
            tau_ki = p.x_ki ^ 2 ./  (2 * p.D_cNO);
            p_NO_k = 0;
            c_NO_k = p.k_O2_k * NO_k.^2 * p.O2_k; % [uM/s]
            d_NO_k = (NO_n - NO_k) ./ tau_nk + (NO_i - NO_k) ./ tau_ki;

            %% Conservation Equations
            % Differential Equations in the Astrocyte
            du(idx.N_K_k, :)    = -J_K_k + 2*J_NaK_k + J_NKCC1_k + J_KCC1_k - J_N_BK_k;
            du(idx.N_Na_k, :)   = -J_Na_k - 3*J_NaK_k + J_NKCC1_k + J_NBC_k;
            du(idx.N_HCO3_k, :) = 2*J_NBC_k;
            du(idx.N_Cl_k, :)   = du(idx.N_Na_k, :) + du(idx.N_K_k, :) - du(idx.N_HCO3_k, :);
            
            % Differential Calcium Equations in Astrocyte
            du(idx.Ca_k, :)     = B_cyt .* (J_IP3 - J_pump + J_ER_leak + J_TRPV_k/p.r_buff);           
            du(idx.s_k, :)      = -(B_cyt .* (J_IP3 - J_pump + J_ER_leak)) ./ (p.VR_ER_cyt);
            du(idx.h_k, :)      = p.k_on * (p.K_inh - (Ca_k + p.K_inh) .* h_k);
            du(idx.I_k, :)      = p.r_h * G - p.k_deg * I_k;
            du(idx.m_k, :)      = p.trpv_switch .* ((minf_k - m_k) ./ p.t_TRPV_k) ; % TRPV open probability, the trpv_switch can switch the turn the channel on or off.
            du(idx.eet_k, :)    = p.V_eet * max(Ca_k - p.Ca_k_min, 0) - p.k_eet * eet_k;
            du(idx.w_k, :)      = phi_w .* (w_inf - w_k);
            
            % Differential Equations in the Perivascular space
            du(idx.K_p, :)      = J_N_BK_k ./ (R_k * p.VR_pa) + J_KIR_i ./ p.VR_ps - p.R_decay * (K_p - p.K_p_min) ;
            du(idx.Ca_p, :)     = (-J_TRPV_k ./ p.VR_pa) + (J_VOCC_k ./ p.VR_ps) - p.Ca_decay_k .* (Ca_p - p.Capmin_k); % calcium concentration in PVS
           
            % Differential Equations in the Synaptic Cleft
            du(idx.N_K_s, :)    = J_K_k - 2 * J_NaK_k - J_NKCC1_k - J_KCC1_k + J_K_NEtoSC_k;
            du(idx.N_Na_s, :)   = -du(idx.N_Na_k, :) - J_K_NEtoSC_k;
            du(idx.N_HCO3_s, :) = -du(idx.N_HCO3_k, :);
 
            % NO pathway:
            du(idx.NO_k, :) = p_NO_k - c_NO_k + d_NO_k;
             
            du = bsxfun(@times, self.enabled, du);
            if nargout == 2

                Uout = zeros(self.n_out, size(u, 2));
                Uout(self.idx_out.J_K_NEtoSC, :) = J_K_NEtoSC;
                Uout(self.idx_out.v_k, :) = v_k;
                Uout(self.idx_out.K_s, :) = K_s;
                Uout(self.idx_out.K_p, :) = K_p;
                Uout(self.idx_out.J_N_BK_k, :) = J_N_BK_k;
                Uout(self.idx_out.rho, :) = rho;
                Uout(self.idx_out.B_cyt, :) = B_cyt;
                Uout(self.idx_out.G, :) = G;
                Uout(self.idx_out.v_3, :) = v_3;
                Uout(self.idx_out.J_IP3, :) = J_IP3;
                Uout(self.idx_out.J_pump, :) = J_pump;
                Uout(self.idx_out.J_ER_leak, :) = J_ER_leak;
                Uout(self.idx_out.J_TRPV_k, :) = J_TRPV_k;
                Uout(self.idx_out.E_BK_k, :) = E_BK_k;
                Uout(self.idx_out.J_NBC_k, :) = J_NBC_k;
                Uout(self.idx_out.E_Na_k, :) = E_Na_k;
                Uout(self.idx_out.E_K_k, :) = E_K_k;
                Uout(self.idx_out.E_NBC_k, :) = E_NBC_k;
                Uout(self.idx_out.E_Cl_k, :) = E_Cl_k;
                Uout(self.idx_out.I_TRPV_k, :) = I_TRPV_k;
                Uout(self.idx_out.J_VOCC_k, :) = J_VOCC_k;
                Uout(self.idx_out.J_K_k, :) = J_K_k;
                Uout(self.idx_out.J_Na_k, :) = J_Na_k;
                Uout(self.idx_out.K_k, :) = K_k;
                Uout(self.idx_out.w_inf, :) = w_inf;
                Uout(self.idx_out.phi_w, :) = phi_w;
                Uout(self.idx_out.J_BK_p, :) = J_BK_p;
                Uout(self.idx_out.J_BK_k, :) = J_BK_k;
                Uout(self.idx_out.E_TRPV_k, :) = E_TRPV_k;
                Uout(self.idx_out.Na_k, :) = Na_k;
                Uout(self.idx_out.J_KCC1_k, :) = J_KCC1_k;
                Uout(self.idx_out.J_NKCC1_k, :) = J_NKCC1_k;
                Uout(self.idx_out.J_NaK_k, :) = J_NaK_k;
                
                varargout = {Uout};
            end
        end
        function [K_p, NO_k] = shared(self, ~, u)
            p = self.params;
            idx = self.index;
            K_p = u(self.index.K_p, :);
            NO_k = u(idx.NO_k, :);
        end
        
       function Glu = input_Glu(self, t) 
            p = self.params;
            Glu = p.GluSwitch * (p.Glu_max - p.Glu_min) * ( ...
            0.5 * tanh((t - p.t_0_Glu) / p.theta_L_Glu) - ...
            0.5 * tanh((t - p.t_2_Glu) / p.theta_R_Glu)) + p.Glu_min;
       end   
        
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end   
end

function idx = indices(self)
% Index of state variables
    idx.R_k = 1;
    idx.K_p = 2;
    idx.N_Na_k = 3;
    idx.N_K_k = 4;
    idx.N_Cl_k = 5;
    idx.N_HCO3_k = 6;
    idx.N_Na_s = 7;
    idx.N_K_s = 8;
    idx.N_HCO3_s = 9;
    idx.w_k = 10;
    idx.I_k = 11;
    idx.Ca_k = 12;
    idx.h_k = 13;
    idx.s_k = 14;
    idx.eet_k = 15;
    idx.m_k = 16;
    idx.Ca_p = 17;
    idx.NO_k = 18;
end

function [idx, n] = output_indices(self)
    % Index of all other output parameters
    idx.K_k  = 1;
    idx.v_k  = 2;
    idx.K_s  = 3;
    idx.K_p  = 4;
    idx.J_N_BK_k  = 5;
    idx.rho  = 6;
    idx.B_cyt  = 7;
    idx.G  = 8;
    idx.v_3  = 9;
    idx.w_inf  = 10;
    idx.phi_w  = 11;
    idx.J_IP3  = 12;
    idx.J_pump  = 13;
    idx.J_ER_leak  = 14;
    idx.J_TRPV_k  = 15;
    idx.E_BK_k  = 16;
    idx.J_NBC_k  = 17;
    idx.E_Na_k  = 18; 
    idx.E_K_k  = 19; 
    idx.E_NBC_k  = 20; 
    idx.E_Cl_k  = 21; 
    idx.I_TRPV_k  = 22; 
    idx.J_VOCC_k  = 23;
    idx.J_K_k  = 24; 
    idx.J_Na_k  = 25;
    idx.J_BK_p = 26;
    idx.J_BK_k = 27;
    idx.E_TRPV_k = 28;
    idx.Cak = 29;
    idx.Na_k = 30;
    idx.J_K_NEtoSC = 31;
    idx.J_KCC1_k = 32;
    idx.J_NKCC1_k = 33;
    idx.J_NaK_k = 34;
    n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
    parser = inputParser();

    % ECS K+ input parameters
    parser.addParameter('ECS_input', 9); 
    parser.addParameter('t0_ECS', 10000);
    parser.addParameter('tend_ECS', 20000);
    
    % Switches to turn on and off some things
    parser.addParameter('rhoSwitch', 1); 
    parser.addParameter('blockSwitch', 1); 
    parser.addParameter('GluSwitch', 1); 
    parser.addParameter('trpv_switch', 1); 
    parser.addParameter('Rk_switch', 1);
    
    % SC decay
    parser.addParameter('Ks_decay', 0.05);      % [s^-1] (M.E.)
    parser.addParameter('Ks_min', 3130);         % [uM] (M.E.)
 

    % Buffer in SC parameters
    parser.addParameter('Mu', 8e-4);            % [m/s]
    parser.addParameter('B0', 500);             % Effective total buffer concentration [mM]
    
    parser.addParameter('rho_min', 0.1);       % microM (one vesicle, Santucci2008)
    parser.addParameter('rho_max', 0.7);  
    
    parser.addParameter('Glu_max', 1846);       % microM (one vesicle, Santucci2008)
    parser.addParameter('Glu_min', 0);          % microM
    parser.addParameter('theta_L_Glu', 1);      % slope of Glu input 
    parser.addParameter('theta_R_Glu', 1);      % slope of Glu input 
    
    parser.addParameter('G_BK_k', 225); % pS (later converted to mho m^-2)
    parser.addParameter('v_4', 8e-3); %V
    parser.addParameter('v_5', 15e-3); %V
    parser.addParameter('v_6', -55e-3); %V
    
    %%% Old parameters %%%
%     parser.addParameter('G_BK_k', 4300); % pS (later converted to mho m^-2)
%     parser.addParameter('v_4', 14.5e-3); %V
%     parser.addParameter('v_5', 8e-3); %V
%     parser.addParameter('v_6', 22e-3); %V
    
    parser.addParameter('reverseBK', 0);
    parser.addParameter('switchBK', 1);
    parser.addParameter('Ca_4', 0.35); % uM 

    % ECS constants
    parser.addParameter('VR_se', 1);
    parser.addParameter('VR_pe', 0.001);
    parser.addParameter('tau', 0.7);
    parser.addParameter('tau2', 2.8);

    % Scaling Constants
    parser.addParameter('L_p', 2.1e-9); % m uM^-1 s^-1
    parser.addParameter('X_k', 12.41e-3); % uM m
    parser.addParameter('R_tot', 8.79e-8); % m

    % Input K+ and glutamate signal
    parser.addParameter('startpulse', 200); % s
    parser.addParameter('lengthpulse', 200); % s
    parser.addParameter('lengtht1', 10); % s
    parser.addParameter('alpha', 2);% [-]
    parser.addParameter('beta', 5);% [-]
    
    % Calcium in the Astrocyte Equations Constants
    parser.addParameter('Amp', 0.7);
    parser.addParameter('base', 0.1);
    parser.addParameter('theta_L', 1);
    parser.addParameter('theta_R', 1);
    parser.addParameter('delta', 1.235e-2); 
    parser.addParameter('BK_switch', -2.5e-10);
    parser.addParameter('VR_ER_cyt', 0.185)
    parser.addParameter('k_on', 2); %uM s^-1
    parser.addParameter('K_inh', 0.1); %uM
    parser.addParameter('r_h', 4.8); % uM
    parser.addParameter('k_deg', 1.25); % s^-1
    parser.addParameter('V_eet', 72); % uM
    parser.addParameter('k_eet', 7.2); % uM
    parser.addParameter('Ca_k_min', 0.1); % uM
    parser.addParameter('Ca_3', 0.4);
    parser.addParameter('eet_shift', 2e-3);
    parser.addParameter('K_I', 0.03); % uM

    %TRPV4
    parser.addParameter('Capmin_k', 2000); %uM
    parser.addParameter('C_astr_k', 40);%pF
    parser.addParameter('gamma_k', 834.3);%mV/uM
    parser.addParameter('gam_cae_k', 200); %uM
    parser.addParameter('gam_cai_k', 0.01); %uM
     parser.addParameter('epshalf_k', 0.1); % 
    parser.addParameter('kappa_k', 0.1);
    parser.addParameter('v1_TRPV_k', 0.120); %V
    parser.addParameter('v2_TRPV_k', 0.013); %V
    parser.addParameter('t_TRPV_k', 0.9); 
    parser.addParameter('R_0_passive_k', 20e-6); 
    parser.addParameter('Ca_decay_k', 0.5);
    parser.addParameter('G_TRPV_k', 50); %pS
    parser.addParameter('r_buff', 0.05); % Rate at which Ca2+ from the TRPV4 channel at the endfoot is buffered compared to rest of channels on the astrocyte body [-]
    
    % Perivascular space
    parser.addParameter('VR_pa', 0.001);% [-]
    parser.addParameter('VR_ps', 0.001);% [-]
    parser.addParameter('R_decay', 0.15);% s^-1
    parser.addParameter('K_p_min', 3e3);% uM

    % Fluxes Constants
    parser.addParameter('F', 9.65e4); %C mol^-1; Faraday's constant
    parser.addParameter('R_g', 8.315); %J mol^-1 K^-1; Gas constant
    parser.addParameter('T', 300); % K; Temperature
    parser.addParameter('g_K_k', 40); %mho m^-2
    parser.addParameter('g_Na_k', 1.314); % mho m^-2
    parser.addParameter('g_NBC_k', 7.57e-1); % mho m^-2
    parser.addParameter('g_KCC1_k', 1e-2); % mho m^-2
    parser.addParameter('g_NKCC1_k', 5.54e-2); % mho m^-2
    parser.addParameter('J_NaK_max', 1.42e-3); % uM m s^-1
    parser.addParameter('K_Na_k', 10000); % uM
    parser.addParameter('K_K_s', 1500); % uM

    parser.addParameter('A_ef_k', 3.7e-9); % m2
    parser.addParameter('C_correction', 1e3); % [-]
    parser.addParameter('J_max', 2880); %uM s^-1
    parser.addParameter('K_act', 0.17); %uM
    parser.addParameter('P_L', 0.0804); %uM
    parser.addParameter('V_max', 20); %uM s^-1
    parser.addParameter('k_pump', 0.24); %uM

    % Additional Equations; Astrocyte Constants
    parser.addParameter('g_Cl_k', 8.797e-1); % mho m^-2
    parser.addParameter('z_K', 1);% [-]
    parser.addParameter('z_Na', 1);% [-]
    parser.addParameter('z_Cl', -1);% [-]
    parser.addParameter('z_NBC', -1);% [-]
    parser.addParameter('z_Ca', 2);% [-]
    parser.addParameter('BK_end', 40);% [-]
    parser.addParameter('K_ex', 0.26); %uM
    parser.addParameter('B_ex', 11.35); %uM
    parser.addParameter('K_G', 8.82); %uM
    parser.addParameter('psi_w', 2.664); %s^-1

    % NO Pathway
    parser.addParameter('D_cNO', 3300);         % [um^2 s^-1] ; Diffusion coefficient NO (Malinski1993)
    parser.addParameter('x_nk', 25);            % [um] ;  (M.E.)
    parser.addParameter('x_ki', 25);            % [um] ;  (M.E.)
    parser.addParameter('k_O2_k', 9.6e-6);      % [uM^-2 s^-1] ;  (Kavdia2002)
    parser.addParameter('O2_k', 200);           % [uM] ;  (M.E.)

    
    parser.parse(varargin{:})
    params = parser.Results;
    params.t_0 = params.startpulse;
    params.t_1 = params.t_0 + params.lengtht1;
    params.t_2 = params.t_0 + params.lengthpulse;
    params.t_3 = params.t_1 + params.lengthpulse;
    params.gab = factorial(params.alpha + params.beta - 1);
    params.ga = factorial(params.alpha - 1);
    params.gb = factorial(params.beta - 1);
    params.t_0_Glu = params.t_0;
    params.t_2_Glu = params.t_2;
    params.startpulse_Glu = params.startpulse;
    params.lengthpulse_Glu = params.lengthpulse;

end

function u0 = initial_conditions(idx,self)
    % Inital estimations of parameters from experimental data
    p = self.params;
    u0 = zeros(length(fieldnames(idx)), 1);
    u0(idx.R_k) = 0.06e-6; %0.0621e-6; 
    u0(idx.N_Na_k) = 0.0010961;
    u0(idx.N_K_k) = 0.0055247;
    u0(idx.N_HCO3_k) = 0.00054791;
    u0(idx.N_Cl_k) = 0.00046402;
    u0(idx.N_Na_s) = 0.00420714;
    u0(idx.N_K_s) = 7.9445e-5;
    u0(idx.N_HCO3_s) = 4.72678e-4;
    u0(idx.K_p) = 3045.1;
    u0(idx.w_k) = 1.703e-4;
    u0(idx.Ca_k) = 0.1612; 
    u0(idx.s_k) = 480.8;
    u0(idx.h_k) = 0.3828;
    u0(idx.I_k) = 0.048299;
    u0(idx.eet_k) = 0.6123;
    u0(idx.m_k) = 0.5710;
    u0(idx.Ca_p) = 1746.4;
    u0(idx.NO_k) = 0.1106;
end
