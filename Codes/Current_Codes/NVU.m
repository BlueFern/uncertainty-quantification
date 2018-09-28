classdef NVU < handle
    % The 'NVU' Code wraps the three subcodes (Astrocyte, SMCEC and Wall
    %   Mechanics) together to be solved by the ode15s solver for stiff
    %   problems. 
    properties
        neuron
        astrocyte
        wall
        smcec
        i_neuron
        i_astrocyte
        i_wall
        i_smcec
        offsets
        outputs
        n
        u0
        T
        U
        odeopts
    end
    methods 
        function self = NVU(neuron, astrocyte, wall, smcec, varargin)
            % Deal with input parameters
            params = parse_inputs(varargin{:});
            names = fieldnames(params);
            for i = 1:numel(names)
                self.(names{i}) = params.(names{i});
            end
            self.neuron = neuron;
            self.astrocyte = astrocyte;
            self.wall = wall;
            self.smcec = smcec;
            
            % Construct mapping to full state vector
            nn = length(fieldnames(self.neuron.index));
            na = length(fieldnames(self.astrocyte.index));
            nw = length(fieldnames(self.wall.index));
            ns = length(fieldnames(self.smcec.index));
            self.offsets = [0, nn, nn+na, nn+na+ns];
            self.outputs = {[],[],[],[]};
            self.i_neuron = 1:nn;
            self.i_astrocyte = nn + (1:na);
            self.i_smcec = nn + na + (1:ns);
            self.i_wall = nn + na + ns + (1:nw);
            self.n = nn + na + ns + nw;      
            
            self.init_conds()
        end
        function du = rhs(self, t, u)
            % Separate out the model components
            un = u(self.i_neuron, :);
            ua = u(self.i_astrocyte, :);
            us = u(self.i_smcec, :);
            uw = u(self.i_wall, :);
            
            % Evaluate the coupling quantities to be passed between
            % submodels as coupling

            [K_p, NO_k] = self.astrocyte.shared(t, ua);
            [Glu, J_K_NEtoSC, NO_n, O2] = self.neuron.shared(t, un);
            [J_KIR_i, Ca_i, J_VOCC_i, NO_i, R_cGMP2] = self.smcec.shared(t, us, K_p); 
            [R, h] = self.wall.shared(t, uw);

            du = zeros(size(u));
            du(self.i_neuron, :) = self.neuron.rhs(t, un, NO_k, R);
            du(self.i_astrocyte, :) = self.astrocyte.rhs(t, ua, J_KIR_i, R, J_VOCC_i, NO_n, NO_i, J_K_NEtoSC, Glu);
            du(self.i_wall, :) = self.wall.rhs(t, uw, Ca_i, R_cGMP2);
            du(self.i_smcec, :) = self.smcec.rhs(t, us, R, h, K_p, NO_k, O2);
        end
        function init_conds(self)
            self.u0 = zeros(self.n, 1);
            self.u0(self.i_neuron) = self.neuron.u0;
            self.u0(self.i_astrocyte) = self.astrocyte.u0;
            self.u0(self.i_smcec) = self.smcec.u0;
            self.u0(self.i_wall) = self.wall.u0;
        end
        function simulate(self)
            self.init_conds()
            f = @(t, u) self.rhs(t, u);
            tStart = tic;
            [self.T, self.U] = ode15s(f, self.T, self.u0, self.odeopts);
            % Now evaluate all of the additional parameters
            un = self.U(:, self.i_neuron).';
            ua = self.U(:, self.i_astrocyte).';
            us = self.U(:, self.i_smcec).';
            uw = self.U(:, self.i_wall).';

            [K_p, NO_k] = self.astrocyte.shared(self.T, ua);
            [Glu, J_K_NEtoSC, NO_n, O2] = self.neuron.shared(self.T, un);
            [J_KIR_i, Ca_i, J_VOCC_i, NO_i, R_cGMP2] = self.smcec.shared(self.T, us, K_p);
            [R, h] = self.wall.shared(self.T, uw);
                     
            [~, self.outputs{1}] = self.neuron.rhs(self.T, un, NO_k, R);
            [~, self.outputs{2}] = self.astrocyte.rhs(self.T, ua, J_KIR_i, R, J_VOCC_i, NO_n, NO_i, J_K_NEtoSC, Glu);
            [~, self.outputs{3}] = self.smcec.rhs(self.T, us, R, h, K_p, NO_k, O2);
            [~, self.outputs{4}] = self.wall.rhs(self.T, uw, Ca_i, R_cGMP2);

            tEnd = toc(tStart);
            fprintf('Elapsed time is %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
        end
        
        % Use this function to run with user specified initial conditions -
        % currently uses the ICs from the previous run to avoid transient
        % behaviour, see nvu_run_script for details
        function simulateManualICs(self)
            f = @(t, u) self.rhs(t, u);
            tStart = tic;
            [self.T, self.U] = ode15s(f, self.T, self.u0, self.odeopts);
            % Now evaluate all of the additional parameters
            un = self.U(:, self.i_neuron).';
            ua = self.U(:, self.i_astrocyte).';
            us = self.U(:, self.i_smcec).';
            uw = self.U(:, self.i_wall).';

            [K_p, NO_k] = self.astrocyte.shared(self.T, ua);
            [Glu, J_K_NEtoSC, NO_n, O2] = self.neuron.shared(self.T, un);
            [J_KIR_i, Ca_i, J_VOCC_i, NO_i, R_cGMP2] = self.smcec.shared(self.T, us, K_p);
            [R, h] = self.wall.shared(self.T, uw);
                     
            [~, self.outputs{1}] = self.neuron.rhs(self.T, un, NO_k, R);
            [~, self.outputs{2}] = self.astrocyte.rhs(self.T, ua, J_KIR_i, R, J_VOCC_i, NO_n, NO_i, J_K_NEtoSC, Glu);
            [~, self.outputs{3}] = self.smcec.rhs(self.T, us, R, h, K_p, NO_k, O2);
            [~, self.outputs{4}] = self.wall.rhs(self.T, uw, Ca_i, R_cGMP2);

            tEnd = toc(tStart);
            fprintf('Elapsed time is %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
        end
        function u = out(self, input_str)
            success = false;
            modules = {self.neuron, self.astrocyte, self.smcec, self.wall};
            for i = 1:4
                module = modules{i};
                if ismember(input_str, fieldnames(module.index))
                   u = self.U(:, self.offsets(i) + (module.index.(input_str)));
                   success = true;
                   break
                end
                if ismember(input_str, fieldnames(module.idx_out))
                    u = self.outputs{i}(module.idx_out.(input_str), :);
                    success = true;
                    break
                end
            end
            if ~success
                error('NVU:InvalidFieldName', 'No matching field: %s', input_str)
            end
        end
        
    end
end

function params = parse_inputs(varargin)
parser = inputParser();
parser.addParameter('odeopts', odeset());
parser.addParameter('T', linspace(0, 1200, 1000));
parser.parse(varargin{:});
params = parser.Results;
end
