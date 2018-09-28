% Plot all variables & fluxes in categories with labels. Automated except
% for subplot dimensions - need to be changed manually. 

%XLIM1 = 1000; XLIM2 = 1010; 

%% Plot state variables

figure(1);
hold on
set(gcf,'Name', 'Neuron State Variables')
neuron_vars = fieldnames(nv.neuron.index);
i_neuron = size(neuron_vars, 1);
for i = 1:1:i_neuron
    subplot(6,6,i)
    plot(nv.T, nv.out(char(neuron_vars(i))));
    xlabel('Time [s]'); ylabel(neuron_vars(i));
    xlim([XLIM1 XLIM2])
end
hold off

figure(2);
hold on
set(gcf,'Name', 'Astrocyte State Variables')
astrocyte_vars = fieldnames(nv.astrocyte.index);
i_astrocyte = size(astrocyte_vars, 1);
for i = 1:1:i_astrocyte
    subplot(4,5,i)
    plot(nv.T, nv.out(char(astrocyte_vars(i))));
    xlabel('Time [s]'); ylabel(astrocyte_vars(i));
    xlim([XLIM1 XLIM2])
end
hold off

figure(3);
hold on
set(gcf,'Name', 'SMC and EC State Variables')
smcec_vars = fieldnames(nv.smcec.index);
i_smcec = size(smcec_vars, 1);
for i = 1:1:i_smcec
    subplot(4,4,i)
    plot(nv.T, nv.out(char(smcec_vars(i))));
    xlabel('Time [s]'); ylabel(smcec_vars(i));
    xlim([XLIM1 XLIM2])
end
hold off

figure(4);
hold on
set(gcf,'Name', 'Wall Mechanics State Variables')
wall_vars = fieldnames(nv.wall.index);
i_wall = size(wall_vars, 1);
for i = 1:1:i_wall
    subplot(2,2,i)
    plot(nv.T, nv.out(char(wall_vars(i))));
    xlabel('Time [s]'); ylabel(wall_vars(i));
    xlim([XLIM1 XLIM2])
end
hold off



%% Fluxes etc.

figure(5);
hold on
set(gcf,'Name', 'Neuron Fluxes/Algebraic Variables')
neuron_flux_vars = fieldnames(nv.neuron.idx_out);
i_neuron_flux = size(neuron_flux_vars, 1);
for i = 1:1:i_neuron_flux
    subplot(6,5,i)
    plot(nv.T, nv.out(char(neuron_flux_vars(i))));
    xlabel('Time [s]'); ylabel(neuron_flux_vars(i));
    xlim([XLIM1 XLIM2])
end
hold off

figure(6);
hold on
set(gcf,'Name', 'Astrocyte Fluxes/Algebraic Variables')
astrocyte_flux_vars = fieldnames(nv.astrocyte.idx_out);
i_astrocyte_flux = size(astrocyte_flux_vars, 1);
for i = 1:1:i_astrocyte_flux
    subplot(6,6,i)
    plot(nv.T, nv.out(char(astrocyte_flux_vars(i))));
    xlabel('Time [s]'); ylabel(astrocyte_flux_vars(i));
    xlim([XLIM1 XLIM2])
end
hold off

figure(7);
hold on
set(gcf,'Name', 'SMC and EC Fluxes/Algebraic Variables')
smcec_flux_vars = fieldnames(nv.smcec.idx_out);
i_smcec_flux = size(smcec_flux_vars, 1);
for i = 1:1:i_smcec_flux
    subplot(6,7,i)
    plot(nv.T, nv.out(char(smcec_flux_vars(i))));
    xlabel('Time [s]'); ylabel(smcec_flux_vars(i));
    xlim([XLIM1 XLIM2])
end
hold off

figure(8);
hold on
set(gcf,'Name', 'Wall Mechanics Fluxes/Algebraic Variables')
wall_flux_vars = fieldnames(nv.wall.idx_out);
i_wall_flux = size(wall_flux_vars, 1);
for i = 1:1:i_wall_flux
    subplot(2,2,i)
    plot(nv.T, nv.out(char(wall_flux_vars(i))));
    xlabel('Time [s]'); ylabel(wall_flux_vars(i));
    xlim([XLIM1 XLIM2])
end
hold off

