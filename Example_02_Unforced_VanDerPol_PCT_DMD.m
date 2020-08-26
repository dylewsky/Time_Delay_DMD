clear variables;

%% Generate input data: Unforced Van der Pol Oscillator

if ~exist(fullfile('Simulation_Data','vdp_data.mat'),'file')
    run(fullfile('Simulation_Data','vdp_sim.m'));
end
load(fullfile('Simulation_Data','vdp_data.mat'));

n = size(y,2);
dt = t(2)-t(1);

%% PCT DMD

nDelay = 128;
delaySteps = 1;
dmd_rank = 10;
dmd_type = 'exact'; % other DMD formulations require 3rd party packages (see documentation)

[Phi,omega,b,U,S,V] = time_delay_dmd(y.',t,nDelay,delaySteps,dmd_rank,...
    'dmd_type',dmd_type,'plot_power_spec',true);

%% DMD Reconstruction
res_delay_step = 1; % which time delay register to read off when collapsing 
                    % from delay coordinates back down to state-space
                    % coordinates
                    
V_recon = zeros(dmd_rank,length(t));
for j = 1:dmd_rank
    V_recon = V_recon + b(j) * Phi(:,j) * exp(omega(j)*t)';
end
H_recon = U(:,1:dmd_rank) * S(1:dmd_rank,1:dmd_rank) * V_recon;
y_recon = real(H_recon((res_delay_step-1)*n + 1:n,:));


figure('Units','Normalized','OuterPosition',[0 0 1 1],'Name',...
    'DMD Reconstruction')
plot_xlim = [0 100];

subplot(2,1,1)
plot(t,y,'LineWidth',1.5);
title('Original Signal')
xlim(plot_xlim)

subplot(2,1,2)
plot(t,y_recon,'LineWidth',1.5);
title(['Rank-' num2str(dmd_rank) ' DMD Reconstruction'])
xlim(plot_xlim)

%% Mode-Separated DMD Reconstruction
figure('Units','Normalized','OuterPosition',[0 0 1 1],'Name',...
    'Mode-Separated DMD Reconstruction')
plot_xlim = [0 100];
plot_ylim = [-4 4];

subplot(7,1,1)
plot(t,y,'LineWidth',1.5);
ylabel('Original Signal')
xlim(plot_xlim)
ylim(plot_ylim)

for j = 1:5
    subplot(7,1,j+1)
    this_V_recon = zeros(dmd_rank,length(t));
    for k = 2*j-[1 0]
        this_V_recon = this_V_recon + b(k) * Phi(:,k) * exp(omega(k)*t)';
    end
    this_H_recon = U(:,1:dmd_rank) * S(1:dmd_rank,1:dmd_rank) * this_V_recon;
    this_y_recon = real(this_H_recon((res_delay_step-1)*n + 1:n,:));
    plot(t,this_y_recon,'LineWidth',1.5)
    xlim(plot_xlim)
    ylim(plot_ylim)
    ylabel({'DMD:',['Modes ' num2str(2*j-1) ' & ' num2str(2*j)]})
end

subplot(7,1,7)
plot(t,y_recon,'LineWidth',1.5)
xlim(plot_xlim)
ylim(plot_ylim)
ylabel({'Sum: All','DMD Modes'})

