clear variables;

%% Generate input data: Van der Pol Oscillator

if ~exist(fullfile('Simulation_Data','vdp_data.mat'),'file')
    run(fullfile('Simulation_Data','vdp_sim.m'));
end
load(fullfile('Simulation_Data','vdp_data.mat'));

n = size(y,2);
dt = t(2)-t(1);

%% Time Delay Embedding

nDelay = 256; % number of delay embeddings
delaySteps = 1; % number of time steps by which each embedding is shifted
                % from the previous one

[H,t_td] = time_delay_embed(y,t,nDelay,delaySteps);

[U,S,V] = svd(H,'econ');
all_PCTs = reshape(U,n,nDelay,n*nDelay);

%% Plot Principal Component Trajectories

figure('Units','Normalized','OuterPosition',[0 0 1 1],'Name',...
    'Principal Component Trajectories: Van der Pol Oscillator')
subplot(4,4,1:2)
plot(dt*(0:(nDelay-1)),y(1:nDelay,:),'LineWidth',1.5)
xlim([0 dt*(nDelay-1)])
title('Input Signal')
subplot(4,4,3:4)
S_cumuluative = cumsum(diag(S))/sum(diag(S));
plot(0:12,[0; S_cumuluative(1:12)],'ko-')
hold on
plot([0 12],[1 1],'k--')
ylim([0 1.1])
title('Singular Value Spectrum (Cumulative)')

for j = 1:12
    subplot(4,4,j+4)
    plot(dt*(0:(nDelay-1)),squeeze(all_PCTs(:,:,j)),'LineWidth',1.5)
    title(['u_{' num2str(j) '}'])
    xlim([0 dt*(nDelay-1)])
end


