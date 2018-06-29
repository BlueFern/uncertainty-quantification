classdef SMCEC < handle
    % The 'SMCEC' code contains the following sections of the model:
    %   The Smooth Muscle cell and the Endothelial Cell
    %   Please refer to the relevient sections in the documentation for 
    %   full information on the equations and variable names.
    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = SMCEC(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, R, h, K_p, NO_k, O2)
            % Initalise inputs and parameters
            p = self.params;
            idx = self.index;
            
            Ca_i = u(idx.Ca_i, :);
            s_i = u(idx.s_i, :);
            v_i = u(idx.v_i, :);
            w_i = u(idx.w_i, :);
            I_i = u(idx.I_i, :);
            NO_i = u(idx.NO_i, :);
            E_b = u(idx.E_b, :);
            E_6c = u(idx.E_6c, :);
            cGMP_i = u(idx.cGMP_i, :);
            
            Ca_j = u(idx.Ca_j, :);
            s_j = u(idx.s_j, :);
            v_j = u(idx.v_j, :);
            I_j = u(idx.I_j, :);
            eNOS_act_j = u(idx.eNOS_act_j, :);
            NO_j = u(idx.NO_j, :);
            
            
            %% SMC fluxes
            J_IP3_i = p.F_i * I_i.^2 ./ (p.K_r_i^2 + I_i.^2);           % IP3 channel
            J_SR_uptake_i = p.B_i * Ca_i.^2 ./ (p.c_b_i^2 + Ca_i.^2);   % SERCA pump
            J_CICR_i = p.C_i * s_i.^2 ./ (p.s_c_i^2 + s_i.^2) .* Ca_i.^4 ./ (p.c_c_i^4 + Ca_i.^4);
            J_extrusion_i = p.D_i * Ca_i .* (1 + (v_i - p.v_d) / p.R_d_i);
            J_SR_leak_i = p.L_i * s_i;
            J_VOCC_i =p.G_Ca_i .* (v_i - p.v_Ca1_i) ./ (1 + exp(-(v_i - p.v_Ca2_i) ./ p.R_Ca_i));
            J_NaCa_i = p.G_NaCa_i * Ca_i ./ (Ca_i + p.c_NaCa_i) .* (v_i - p.v_NaCa_i);
            J_stretch_i = p.G_stretch ./ (1 + exp(-p.alpha_stretch*(p.trans_p_mmHg*R./h - p.sigma_0))) .* (v_i - p.E_SAC);
            J_Cl_i = p.G_Cl_i * (v_i - p.v_Cl_i);
            J_NaK_i = p.F_NaK_i;
            J_K_i   = p.G_K_i * w_i .* (v_i - p.v_K_i);

            [J_KIR_i] = self.shared(t, u, K_p);
            
            J_degrad_i = p.k_d_i * I_i;
            
            %% EC fluxes
            J_IP3_j = p.F_j * I_j.^2 ./ (p.K_r_j^2 + I_j.^2);
            J_ER_uptake_j = p.B_j * Ca_j.^2 ./ (p.c_b_j^2 + Ca_j.^2);  
            J_CICR_j = p.C_j * s_j.^2 ./ (p.s_c_j^2 + s_j.^2) .* Ca_j.^4 ./ (p.c_c_j^4 + Ca_j.^4);
            J_extrusion_j = p.D_j * Ca_j;
            J_stretch_j = p.G_stretch ./ (1 + exp(-p.alpha_stretch*(p.trans_p_mmHg*R./h - p.sigma_0))) .* (v_j - p.E_SAC);
            J_ER_leak_j = p.L_j * s_j;
            J_cation_j = p.G_cat_j * (p.E_Ca_j - v_j) * 0.5 .* (1 + tanh((log10(Ca_j) - p.m_3_cat_j) / p.m_4_cat_j));
            J_BK_Ca_j = 0.2 * (1 + tanh( ((log10(Ca_j) - p.c) .* (v_j - p.bb_j) - p.a_1_j) ./ (p.m_3b_j * (v_j + p.a_2_j*(log10(Ca_j) - p.c) - p.bb_j).^2 + p.m_4b_j)));
            J_SK_Ca_j = 0.3 * (1 + tanh((log10(Ca_j) - p.m_3s_j) / p.m_4s_j));
            J_K_j = p.G_tot_j * (v_j - p.v_K_j) .* (J_BK_Ca_j + J_SK_Ca_j);
            J_R_j = p.G_R_j * (v_j - p.v_rest_j);
            J_degrad_j = p.k_d_j * I_j;
            
            %% Coupling
            V_coup_i = -p.G_coup * (v_i - v_j);
            J_IP3_coup_i = -p.P_IP3 * (I_i - I_j);
            J_Ca_coup_i = -p.P_Ca * (Ca_i - Ca_j);

            c_w_i = 1/2 * (1 + tanh((cGMP_i - 10.75)/0.668) );
            
            K_act_i = (Ca_i + c_w_i).^2 ./ ((Ca_i + c_w_i).^2 + p.beta_i * exp(-(v_i - p.v_Ca3_i) / p.R_K_i)); 
 
            tau_wss = R/2 * p.delta_p_L; % from Dormanns 2016
            
            % NO pathway 
            tau_ki = p.x_ki ^ 2 ./  (2 * p.D_cNO);
            tau_ij = p.x_ij ^ 2 ./  (2 * p.D_cNO);
            p_NO_i = 0;
            c_NO_i = p.k_dno * NO_i;
            d_NO_i = (NO_k - NO_i) ./ tau_ki + (NO_j - NO_i) ./ tau_ij;
            
            k4 = p.C_4 * cGMP_i.^2;
            E_5c = 1 - E_b - E_6c;

            V_max_pde = p.k_pde * cGMP_i;

%             p_NO_j = p.NOswitch * ( p.V_NOj_max * eNOS_act_j * p.O2_j / (p.K_mO2_j + p.O2_j) * p.LArg_j / (p.K_mArg_j + p.LArg_j) );
%             c_NO_j = p.k_O2 * NO_j.^2 * p.O2_j;

            O2_j = O2*1e3;  % Oxygen in EC taken as O2 from lumen (diffusion very fast so plausible!) instead of constant
            p_NO_j = p.NOswitch * ( p.V_NOj_max .* eNOS_act_j .* O2_j / (p.K_mO2_j + O2_j) .* p.LArg_j / (p.K_mArg_j + p.LArg_j) );
            c_NO_j = p.k_O2 .* NO_j.^2 .* O2_j;
            
            J_lumen = - NO_j * 4 * p.D_cNO ./ (25.^2); 
            d_NO_j = (NO_i - NO_j) ./ tau_ij + J_lumen; 

            W_wss = p.W_0 * (tau_wss + sqrt(16 * p.delta_wss^2 + tau_wss.^2) - 4 * p.delta_wss).^2 / (tau_wss + sqrt(16 * p.delta_wss^2 + tau_wss.^2)); 
            F_wss = 1 / (1 + p.alp * exp(-W_wss)) - 1 / (1 + p.alp); 

            Act_eNOS_Ca = p.K_dis * Ca_j / (p.K_eNOS + Ca_j); 
            Act_eNOS_wss = p.g_max * F_wss;  
            
            
            %% Differential Equations
            % Smooth muscle cell
            du(idx.Ca_i, :) = J_IP3_i - J_SR_uptake_i - J_extrusion_i + J_SR_leak_i - J_VOCC_i + J_CICR_i + J_NaCa_i - 0.1*J_stretch_i + J_Ca_coup_i;
            
            du(idx.s_i, :) = J_SR_uptake_i - J_CICR_i - J_SR_leak_i;
            du(idx.v_i, :) = p.gamma_i * ( -J_NaK_i - J_Cl_i - 2*J_VOCC_i - J_NaCa_i - J_K_i -J_stretch_i - J_KIR_i) + V_coup_i;
            du(idx.w_i, :) = p.lambda_i * (K_act_i - w_i);
            du(idx.I_i, :) = J_IP3_coup_i - J_degrad_i;
            du(idx.K_i, :) = J_NaK_i - J_KIR_i - J_K_i;
            
            % Endothelial Cell
            du(idx.Ca_j, :) = J_IP3_j - J_ER_uptake_j + J_CICR_j - J_extrusion_j + J_ER_leak_j + J_cation_j + p.J_0_j - J_stretch_j - J_Ca_coup_i;
            du(idx.s_j, :) = J_ER_uptake_j - J_CICR_j - J_ER_leak_j;
            du(idx.v_j, :) = -1/p.C_m_j * (J_K_j + J_R_j) - V_coup_i;
            du(idx.I_j, :) = p.J_PLC - J_degrad_j - J_IP3_coup_i;           % p.J_PLC or self.input_plc(t)
            
            % NO pathway
            du(idx.NO_i, :) = p_NO_i - c_NO_i + d_NO_i;
            du(idx.E_b, :)  = -p.k1 * E_b .* NO_i + p.k_1 * E_6c + k4 .* E_5c;     
            du(idx.E_6c, :) = p.k1 * E_b .* NO_i - (p.k_1 + p.k2) * E_6c - p.k3 * E_6c .* NO_i;
            du(idx.cGMP_i, :) = p.V_max_sGC * E_5c - V_max_pde .* cGMP_i ./ (p.K_m_pde + cGMP_i); 
            du(idx.eNOS_act_j, :) = (p.gam_eNOS * Act_eNOS_Ca  + (1 - p.gam_eNOS) * Act_eNOS_wss - p.mu2_j * eNOS_act_j); 
            du(idx.NO_j, :) = p_NO_j - c_NO_j + d_NO_j;
            
            du = bsxfun(@times, self.enabled, du);
            
            if nargout == 2
                Uout = zeros(self.n_out, size(u, 2));
                Uout(self.idx_out.V_coup_i, :) = V_coup_i;
                Uout(self.idx_out.J_Ca_coup_i, :) = J_Ca_coup_i;
                Uout(self.idx_out.J_IP3_coup_i, :) = J_IP3_coup_i;
                Uout(self.idx_out.J_stretch_i, :) = J_stretch_i;
                Uout(self.idx_out.J_IP3_i, :) = J_IP3_i;
                Uout(self.idx_out.J_SR_uptake_i, :) = J_SR_uptake_i;
                Uout(self.idx_out.J_CICR_i, :) = J_CICR_i;
                Uout(self.idx_out.J_extrusion_i, :) = J_extrusion_i;
                Uout(self.idx_out.J_SR_leak_i, :) = J_SR_leak_i;
                Uout(self.idx_out.J_VOCC_i, :) = J_VOCC_i;
                Uout(self.idx_out.J_NaCa_i, :) = J_NaCa_i;
                Uout(self.idx_out.J_NaK_i, :) = J_NaK_i;
                Uout(self.idx_out.J_Cl_i, :) = J_Cl_i;
                Uout(self.idx_out.J_K_i, :) = J_K_i;
                Uout(self.idx_out.J_KIR_i, :) = J_KIR_i;
                Uout(self.idx_out.K_act_i, :) = K_act_i;
                Uout(self.idx_out.J_degrad_i, :) = J_degrad_i;
                
                Uout(self.idx_out.V_coup_j, :) = -V_coup_i;
                Uout(self.idx_out.J_Ca_coup_j, :) = -J_Ca_coup_i;
                Uout(self.idx_out.J_IP3_coup_j, :) = -J_IP3_coup_i;
                Uout(self.idx_out.J_stretch_j, :) = J_stretch_j;
                Uout(self.idx_out.J_0_j, :) = p.J_0_j;
                Uout(self.idx_out.J_IP3_j, :) = J_IP3_j;
                Uout(self.idx_out.J_ER_uptake_j, :) = J_ER_uptake_j;
                Uout(self.idx_out.J_CICR_j, :) = J_CICR_j;
                Uout(self.idx_out.J_extrusion_j, :) = J_extrusion_j;
                Uout(self.idx_out.J_ER_leak_j, :) = J_ER_leak_j;
                Uout(self.idx_out.J_cation_j, :) = J_cation_j;
                Uout(self.idx_out.J_BK_Ca_j, :) = J_BK_Ca_j;
                Uout(self.idx_out.J_SK_Ca_j, :) = J_SK_Ca_j;
                Uout(self.idx_out.J_K_j, :) = J_K_j;
                Uout(self.idx_out.J_R_j, :) = J_R_j;
                Uout(self.idx_out.J_degrad_j, :) = J_degrad_j;
                Uout(self.idx_out.c_w_i, :) = c_w_i;
                Uout(self.idx_out.tau_wss, :) = tau_wss;
                Uout(self.idx_out.E_5c, :) = E_5c;
                Uout(self.idx_out.d_NO_j, :) = d_NO_j;
                Uout(self.idx_out.J_lumen, :) = J_lumen;
                Uout(self.idx_out.p_NO_j, :) = p_NO_j;
                Uout(self.idx_out.c_NO_j, :) = -c_NO_j;
                varargout{1} = Uout; 
            end
        end
        
        function [J_KIR_i, Ca_i, J_VOCC_i, NO_i, R_cGMP2] = shared(self, t, u, K_p)
            p = self.params;
            idx = self.index;
            v_i = u(idx.v_i, :);
            Ca_i = u(idx.Ca_i, :);
            NO_i = u(idx.NO_i, :);
            cGMP_i = u(idx.cGMP_i, :);

            R_cGMP2 = cGMP_i.^2 ./ (cGMP_i.^2 + p.K_m_mlcp^2);
            
            v_KIR_i = p.z_1 * K_p - p.z_2;
            g_KIR_i = exp(p.z_5 * v_i + p.z_3 * K_p - p.z_4);
            J_KIR_i = p.F_KIR_i * g_KIR_i / p.gamma_i .* (v_i - v_KIR_i);
            
            J_VOCC_i = p.G_Ca_i .* (v_i - p.v_Ca1_i) ./ (1 + exp(-(v_i - p.v_Ca2_i) ./ p.R_Ca_i));
        end
        
        function jplc = input_plc(self, t)
            % The PLC input (stable to oscillatory)
            PLC_min = 0.18; PLC_max = 0.4; 
            t_up = 300; 
            t_down = 800;
            jplc = PLC_min + (PLC_max - PLC_min) * (0.5 * tanh((t - t_up) / 0.05) - 0.5 * tanh((t - t_down) / 0.05));
            jplc = jplc(:).';
        end
        
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
        
        
    end
end

function idx = indices()
    % Index of parameters needing inital conditions 
    idx.Ca_i = 1;
    idx.s_i = 2;
    idx.v_i = 3;
    idx.w_i = 4;
    idx.I_i = 5;
    idx.K_i = 6;
    idx.NO_i = 7;
    idx.E_b = 8;
    idx.E_6c = 9;
    idx.cGMP_i = 10;
    idx.Ca_j = 11;
    idx.s_j = 12;
    idx.v_j = 13;
    idx.I_j = 14;
    idx.eNOS_act_j = 15;
    idx.NO_j = 16;
    
end

function [idx, n] = output_indices()
    % Index of all other output parameters
    idx.V_coup_i = 1;
    idx.J_Ca_coup_i = 2;
    idx.J_IP3_coup_i = 3;
    idx.J_stretch_i = 4;
    idx.J_IP3_i = 5;
    idx.J_SR_uptake_i = 6;
    idx.J_CICR_i = 7;
    idx.J_extrusion_i = 8;
    idx.J_SR_leak_i = 9;
    idx.J_VOCC_i = 10;
    idx.J_NaCa_i = 11;
    idx.J_NaK_i = 12;
    idx.J_Cl_i = 13;
    idx.J_K_i = 14;
    idx.J_KIR_i = 15;
    idx.K_act_i = 16;
    idx.J_degrad_i = 17;
    idx.V_coup_j = 18;
    idx.J_Ca_coup_j = 19;
    idx.J_IP3_coup_j = 20;
    idx.J_stretch_j = 21;
    idx.J_0_j = 22;
    idx.J_IP3_j = 23;
    idx.J_ER_uptake_j = 24;
    idx.J_CICR_j = 25;
    idx.J_extrusion_j = 26;
    idx.J_ER_leak_j = 27;
    idx.J_cation_j = 28;
    idx.J_BK_Ca_j = 29;
    idx.J_SK_Ca_j = 30;
    idx.J_K_j = 31;
    idx.J_R_j = 32;
    idx.J_degrad_j = 33;
    idx.c_w_i = 34;
    idx.tau_wss = 35;
    idx.E_5c = 36;
    idx.d_NO_j = 37;
    idx.J_lumen = 38;
    idx.p_NO_j = 39;
    idx.c_NO_j = 40;
    
    n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
    parser = inputParser();
    
    parser.addParameter('NOswitch', 1); 
    
    % Smooth Muscle Cell ODE Constants
    parser.addParameter('gamma_i', 1970); %mV uM^-1
    parser.addParameter('lambda_i', 45); % s^-1

    % Endothelial Cell ODE Constants
    parser.addParameter('C_m_j', 25.8); %pF
    parser.addParameter('J_PLC', 0.18); % uMs^-1
    parser.addParameter('J_0_j', 0.029); %constant Ca influx (EC)uMs^-1

    % Smooth Muscle Cell Flux Constants
    parser.addParameter('F_i', 0.23); %uM s^-1      % IP3/RYR channel strength
    parser.addParameter('K_r_i', 1); %uM

    parser.addParameter('B_i', 2.025); %uM s^-1     % SERCA pump strength     
    parser.addParameter('c_b_i', 1.0); %uM

    parser.addParameter('C_i', 55); % uM s^-1
    parser.addParameter('s_c_i', 2.0); %uM
    parser.addParameter('c_c_i', 0.9); %uM

    parser.addParameter('D_i', 0.24); %s^-1
    parser.addParameter('v_d', -100); %mV
    parser.addParameter('R_d_i', 250); %mV

    parser.addParameter('L_i', 0.025); %s^-1

    parser.addParameter('G_Ca_i', 1.29e-3); %uM mV^-1 s^-1
    parser.addParameter('v_Ca1_i', 100); %mV
    parser.addParameter('v_Ca2_i', -24); %mV
    parser.addParameter('R_Ca_i', 8.5); %mV

    parser.addParameter('G_NaCa_i', 3.16e-3); %uM mV^-1 s^-1
    parser.addParameter('c_NaCa_i', 0.5); %uM
    parser.addParameter('v_NaCa_i', -30); %mV

    parser.addParameter('G_stretch', 6.1e-3); % uM mV^-1 s^-1   (Also EC parameter)
    parser.addParameter('alpha_stretch', 7.4e-3); % mmHg^-1     (Also EC parameter)
    parser.addParameter('trans_p_mmHg', 30); % mmHg                  (Also EC parameter) transmural pressure. 30 mmHg = 4000 Pa
    parser.addParameter('sigma_0', 500); % mmHg                 (Also EC parameter)
    parser.addParameter('E_SAC', -18); % mV                     (Also EC parameter)

    parser.addParameter('F_NaK_i', 4.32e-2); %uM s^-1

    parser.addParameter('G_Cl_i', 1.34e-3); %uM mV^-1 s^-1
    parser.addParameter('v_Cl_i', -25); %mV

    parser.addParameter('G_K_i', 4.46e-3); %uM mV^-1 s^-1
    parser.addParameter('v_K_i', -94); %mV

    parser.addParameter('F_KIR_i', 7.5e2); % [-]
    parser.addParameter('k_d_i', 0.1); % s^-1

    % Endothelial Cell Flux Constants
    parser.addParameter('F_j', 0.23); %uM s^-1
    parser.addParameter('K_r_j', 1); %uM
    parser.addParameter('B_j', 0.5); %uM s^-1
    parser.addParameter('c_b_j', 1); %uM
    parser.addParameter('C_j', 5); %uM s^-1
    parser.addParameter('s_c_j', 2); %uM
    parser.addParameter('c_c_j', 0.9); %uM
    parser.addParameter('D_j', 0.24);% s^-1

    % (G_stretch, alpha_stretch, trans_p_mmHg, sigma0, E_SAC are included above in 
    %  SMC flux Constants)

    parser.addParameter('L_j', 0.025); %s^-1 

    parser.addParameter('G_cat_j', 6.6e-4); %uM mV^-1 s^-1
    parser.addParameter('E_Ca_j', 50); %mV
    parser.addParameter('m_3_cat_j', -0.18); %uM
    parser.addParameter('m_4_cat_j', 0.37); %uM

    parser.addParameter('G_tot_j', 6927); %p mho
    parser.addParameter('v_K_j', -80); %m mho

    parser.addParameter('c', -0.4); %uM
    parser.addParameter('bb_j', -80.8); %mV
    parser.addParameter('a_1_j', 53.3); %uM mV
    parser.addParameter('a_2_j', 53.3); % mV uM^-1
    parser.addParameter('m_3b_j', 1.32e-3); %uM mV^-1
    parser.addParameter('m_4b_j', 0.3); %uM mV
    parser.addParameter('m_3s_j', -0.28); %uM
    parser.addParameter('m_4s_j', 0.389); %uM

    parser.addParameter('G_R_j', 955); %p omh
    parser.addParameter('v_rest_j', -31.1); %mV
    parser.addParameter('k_d_j', 0.1); %s^-1

    parser.addParameter('P_Ca', 0.05); %s^-1
    parser.addParameter('P_IP3', 0.05); %s^-1
    parser.addParameter('G_coup', 0.5); %s^-1

    % Additional Equations Constants
    parser.addParameter('beta_i', 0.13); %uM^2
    parser.addParameter('v_Ca3_i', -27); %mV
    parser.addParameter('R_K_i', 12); %mV
    parser.addParameter('z_1', 4.5e-3); %mV
    parser.addParameter('z_2', 112); %mV
    parser.addParameter('z_3', 4.2e-4); %uM mV^-1 s^-1
    parser.addParameter('z_4', 12.6); %uM mV^-1 s^-1
    parser.addParameter('z_5', -7.4e-2); %uM mV^-1 s^-1
    
    % NO pathway
    parser.addParameter('D_cNO', 3300); % [um^2 s^-1] ; Diffusion coefficient NO (Malinski1993)
    parser.addParameter('K_mArg_j', 1.5); % [] ;
    parser.addParameter('K_mO2_j', 7.7); % [] ; Chen2006
    parser.addParameter('k_dno', 0.01); % [s^-1] ;
    parser.addParameter('K_m_mlcp', 5.5); % [uM] ;
    parser.addParameter('V_NOj_max', 1.22); % [s^-1] ; maximum catalytic rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18
    parser.addParameter('O2_j', 200); % [uM] ; O2 concentration in the EC (ME)
    parser.addParameter('LArg_j', 100); % [uM] ;
    parser.addParameter('k_O2', 9.6e-6); % [uM^-2 s^-1] ;
    parser.addParameter('W_0', 1.4); % [Pa^-1] ; shear gating constant (Comerford2008)
    parser.addParameter('delta_wss', 2.86); % [Pa] ; the membrane shear modulus (Comerford2008)
    parser.addParameter('k_1', 100); % [s^{-1}] ;
    parser.addParameter('k1', 2e3); % [uM^-1 s^-1] ;
    parser.addParameter('k2', 0.1); % [s^-1] ;
    parser.addParameter('k3', 3); % [uM^-1 s^-1] ;
    parser.addParameter('V_max_sGC', 0.8520); % [] ;
    parser.addParameter('k_pde', 0.0195); % [s^-1] ;
    parser.addParameter('C_4', 0.011); % [s^{-1} microM^{-2}] ;
    parser.addParameter('K_m_pde', 2); % [uM] ;
    parser.addParameter('gam_eNOS', 0.1); % [-] ;
    parser.addParameter('mu2_j', 0.0167); % [s^-1] ;
    parser.addParameter('K_dis', 9e-2); % [uM s^-1] ;
    parser.addParameter('K_eNOS', 4.5e-1); % [uM] ;    
    parser.addParameter('g_max', 0.06); % [uM s^-1] ;    
    parser.addParameter('alp', 2); % [-] ; zero shear open channel constant (Comerford2008); in Wiesner1997: alp = 3
    parser.addParameter('delta_p_L', 9.1e4); % 9.1e4: ME
    parser.addParameter('x_ki', 25); % [um]  (M.E.)
    parser.addParameter('x_ij', 3.75); % [um]  (Kavdia2002)

    
    parser.parse(varargin{:})
    params = parser.Results;
end

function u0 = initial_conditions(idx)
    u0 = zeros(length(fieldnames(idx)), 1);

    u0(idx.Ca_i) = 0.2637; % 0.1;
    u0(idx.s_i) = 1.1686; % 0.1;
    u0(idx.v_i) = -34.7; % -60;
    u0(idx.w_i) = 0.2206; % 0.1;
    u0(idx.I_i) = 0.275; % 0.1;
    u0(idx.K_i) = 99994.8; % 100e3;
    u0(idx.NO_i) = 0.0541; % 0.05;
    u0(idx.E_b) = 0.4077; % 1/3;
    u0(idx.E_6c) = 0.4396; % 1/3;
    u0(idx.cGMP_i) = 8.2826; % 8;
    u0(idx.Ca_j) = 0.8331; % 0.1;
    u0(idx.s_j) = 0.6266; % 0.1;
    u0(idx.v_j) = -68.27; % -75;
    u0(idx.I_j) = 0.825; % 0.1;
    u0(idx.eNOS_act_j) = 0.4479; % 0.7;
    u0(idx.NO_j) = 0.0528; % 0.05;

end

