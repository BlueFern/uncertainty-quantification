classdef WallMechanics < handle
    % The 'WallMechanics' code contains the following sections of the model:
    %   The Contraction and Mechanical Models.
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
        function self = WallMechanics(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, Ca_i, R_cGMP2)
            % Initalise inputs and parameters
            idx = self.index;
            p = self.params;
            % Make the output function compute these so as to avoid
            % duplicating equations
            [R, h] = self.shared(t, u);

            Mp = u(idx.Mp, :);
            AMp = u(idx.AMp, :);
            AM = u(idx.AM, :);
            
            %% Contraction Equations
            K_1 = p.gamma_cross * Ca_i.^p.n_cross;
            K_6 = K_1;
            K_2 = 58.1395 * p.k_mlcp_b + 58.1395 * p.k_mlcp_c * R_cGMP2;
            K_5 = K_2;
            
            M = 1 - AM - AMp - Mp;
            du(idx.Mp, :) = p.wallMech * ( p.K_4 * AMp + K_1 .* M - (K_2 + p.K_3) .* Mp );
            du(idx.AMp, :) = p.wallMech  * ( p.K_3 * Mp + K_6 .* AM - (p.K_4 + K_5) .* AMp );
            du(idx.AM, :) = p.wallMech  * ( K_5 .* AMp - (p.K_7 + K_6) .* AM );
            
            % Mechanical Equations
            F_r = AMp + AM;
            E = p.E_passive + F_r * (p.E_active - p.E_passive);
            R_0 = p.R_0_passive + F_r * (p.alpha - 1) * p.R_0_passive;

            du(idx.R, :) =  p.R_0_passive / p.eta * ( R * p.trans_p ./ h - E .* (R - R_0) ./ R_0);

            du = bsxfun(@times, self.enabled, du);
            
            if nargout == 2
               Uout = zeros(self.n_out, size(u, 2));
               Uout(self.idx_out.M, :) = M;
               Uout(self.idx_out.F_r, :) = F_r;
               Uout(self.idx_out.R, :) = R;
               Uout(self.idx_out.K_2, :) = K_2;
               varargout{1} = Uout;
            end
        end 

        function [R, h] = shared(self, ~,u)
           
            R = u(self.index.R, :);
            h = 0.1 * R; 
        end  

        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end 
end

function idx = indices()
    % Index of parameters needing inital conditions 
    idx.Mp = 1;
    idx.AMp = 2;
    idx.AM = 3;
    idx.R = 4;
end

function [idx, n] = output_indices()
    % Index of all other output parameters
    idx.M = 1;
    idx.F_r = 2;
    idx.R = 3;
    idx.K_2 = 4;
    n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
    parser = inputParser();
    
    % Parameter for changing the wall mechanics rate constants
    parser.addParameter('wallMech', 1);
    
    % Contraction Equation Constants
    parser.addParameter('K_3', 0.4); % s^-1
    parser.addParameter('K_4', 0.1); % s^-1
    parser.addParameter('K_7', 0.1); % s^-1
    parser.addParameter('gamma_cross', 17); %uM^-3 s^-1
    parser.addParameter('n_cross', 3); % fraction constant of the phosphorylation crossbridge
    
    % Mechanical Equation Constants
    parser.addParameter('eta', 1e4); %Pa s
    parser.addParameter('R_0_passive', 20e-6); % m
    parser.addParameter('trans_p', 4000); % Pa  transmural pressure
    parser.addParameter('E_passive', 66e3); % Pa
    parser.addParameter('E_active', 233e3); % Pa
    parser.addParameter('alpha', 0.6); % [-]
    parser.addParameter('k_mlcp_b', 0.0086); % [s^-1]
    parser.addParameter('k_mlcp_c', 0.0327); % [s^-1]

    parser.parse(varargin{:});
    params = parser.Results;

end
function u0 = initial_conditions(idx)
    u0 = zeros(length(fieldnames(idx)), 1);
    % Inital estimations of parameters from experimental data
    u0(idx.Mp) = 0.0842;
    u0(idx.AMp) = 0.0622;
    u0(idx.AM) = 0.2746;
    u0(idx.R) = 22.97e-6;
end
