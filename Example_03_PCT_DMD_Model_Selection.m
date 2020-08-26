clear variables;

%% Generate input data: Unforced Van der Pol Oscillator

if ~exist(fullfile('Simulation_Data','vdp_data.mat'),'file')
    run(fullfile('Simulation_Data','vdp_sim.m'));
end
load(fullfile('Simulation_Data','vdp_data.mat'));

n = size(y,2);
dt = t(2)-t(1);

%% PCT DMD with varied delay embedding dimension
nDelay = 2.^(4:7);
delaySteps = 1;
dmd_rank = 10;
dmd_type = 'exact'; % other DMD formulations require 3rd party packages (see documentation)

[all_Phi,all_omega,all_b,all_U,all_S,all_V] = time_delay_dmd(y.',t,nDelay,delaySteps,dmd_rank,...
    'dmd_type',dmd_type);

figure('Units','Normalized','OuterPosition',[0 0 1 1],'Name',...
    'DMD Reconstructions: Varied nDelay')
plot_xlim = [0 100];

subplot(length(nDelay)+1,1,1)
plot(t,y,'LineWidth',1.5);
title('Original Signal')
xlim(plot_xlim)

res_delay_step = 1; % which time delay register to read off when collapsing 
                    % from delay coordinates back down to state-space
                    % coordinates

for j = 1:length(nDelay)
    Phi = all_Phi{j}; omega = all_omega{j}; b = all_b{j};
    U = all_U{j}; S = diag(all_S{j}); V = all_V{j};
    V_recon = zeros(dmd_rank,length(t));
    for k = 1:dmd_rank
        V_recon = V_recon + b(k) * Phi(:,k) * exp(omega(k)*t)';
    end
    H_recon = U(:,1:dmd_rank) * S(1:dmd_rank,1:dmd_rank) * V_recon;
    y_recon = real(H_recon((res_delay_step-1)*n + 1:n,:));

    subplot(length(nDelay)+1,1,j+1)
    plot(t,y_recon,'LineWidth',1.5);
    title(['DMD Reconstruction: nDelay = ' num2str(nDelay(j))])
    xlim(plot_xlim)
end

%% PCT DMD with varied rank
nDelay = 128;
delaySteps = 1;
dmd_rank_list = 2.^(1:4);
dmd_type = 'exact'; % other DMD formulations require 3rd party packages (see documentation)

[all_Phi,all_omega,all_b,all_U,all_S,all_V] = time_delay_dmd(y.',t,nDelay,delaySteps,dmd_rank_list,...
    'dmd_type',dmd_type);

figure('Units','Normalized','OuterPosition',[0 0 1 1],'Name',...
    'DMD Reconstructions: Varied nDelay')
plot_xlim = [0 100];

subplot(length(dmd_rank_list)+1,1,1)
plot(t,y,'LineWidth',1.5);
title('Original Signal')
xlim(plot_xlim)

res_delay_step = 1; % which time delay register to read off when collapsing 
                    % from delay coordinates back down to state-space
                    % coordinates

for j = 1:length(dmd_rank_list)
    dmd_rank = dmd_rank_list(j);
    Phi = all_Phi{j}; omega = all_omega{j}; b = all_b{j};
    U = all_U{j}; S = diag(all_S{j}); V = all_V{j};
    V_recon = zeros(dmd_rank,length(t));
    for k = 1:dmd_rank
        V_recon = V_recon + b(k) * Phi(:,k) * exp(omega(k)*t)';
    end
    H_recon = U(:,1:dmd_rank) * S(1:dmd_rank,1:dmd_rank) * V_recon;
    y_recon = real(H_recon((res_delay_step-1)*n + 1:n,:));

    subplot(length(dmd_rank_list)+1,1,j+1)
    plot(t,y_recon,'LineWidth',1.5);
    title(['DMD Reconstruction: Rank = ' num2str(dmd_rank_list(j))])
    xlim(plot_xlim)
end

%% PCT DMD with varied delay embedding dimension and step size
% Here nDelay and delaySteps are varied such that the embedding period
% stays constant (of duration 256*dt), but the number of data points
% spanning that duration varies
%
% From the plot it will be apparent that reconstructions are not strongly
% affected by this variation as long as nDelay is large enough to
% accommodate the desired DMD rank (nDelay >= rank)

nDelay_list = 2.^(4:8);
delaySteps_list = 256./nDelay_list;
dmd_rank_list = 10*ones(size(nDelay_list));
dmd_type = 'exact'; % other DMD formulations require 3rd party packages (see documentation)

[all_Phi,all_omega,all_b,all_U,all_S,all_V] = time_delay_dmd(y.',t,nDelay,delaySteps,dmd_rank_list,...
    'dmd_type',dmd_type);

figure('Units','Normalized','OuterPosition',[0 0 1 1],'Name',...
    'DMD Reconstructions: Varied nDelay')
plot_xlim = [0 100];

subplot(length(dmd_rank_list)+1,1,1)
plot(t,y,'LineWidth',1.5);
title('Original Signal')
xlim(plot_xlim)

res_delay_step = 1; % which time delay register to read off when collapsing 
                    % from delay coordinates back down to state-space
                    % coordinates

for j = 1:length(dmd_rank_list)
    dmd_rank = dmd_rank_list(j);
    nDelay = nDelay_list(j);
    delaySteps = delaySteps_list(j);
    Phi = all_Phi{j}; omega = all_omega{j}; b = all_b{j};
    U = all_U{j}; S = diag(all_S{j}); V = all_V{j};
    V_recon = zeros(dmd_rank,length(t));
    for k = 1:dmd_rank
        V_recon = V_recon + b(k) * Phi(:,k) * exp(omega(k)*t)';
    end
    H_recon = U(:,1:dmd_rank) * S(1:dmd_rank,1:dmd_rank) * V_recon;
    y_recon = real(H_recon((res_delay_step-1)*n + 1:n,:));

    subplot(length(dmd_rank_list)+1,1,j+1)
    plot(t,y_recon,'LineWidth',1.5);
    title(['DMD Reconstruction: nDelay = ' num2str(nDelay_list(j)) ', delaySteps = ' num2str(delaySteps_list(j))])
    xlim(plot_xlim)
end


