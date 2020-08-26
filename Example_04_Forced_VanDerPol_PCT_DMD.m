clear variables;

%% Generate input data: Forced Van der Pol Oscillator

if ~exist(fullfile('Simulation_Data','vdp_forced_lorenz_data.mat'),'file')
    run(fullfile('Simulation_Data','vdp_forced_lorenz_sim.m'));
end
load(fullfile('Simulation_Data','vdp_forced_lorenz_data.mat'));

n = size(y,2);
dt = t(2)-t(1);

nDelay = 128;
delaySteps = 1;
dmd_rank = 10;
dmd_type = 'exact'; % other DMD formulations require 3rd party packages (see documentation)

[Phi,omega,b,U,S,V] = time_delay_dmd(y.',t,nDelay,delaySteps,dmd_rank,...
    'dmd_type',dmd_type,'plot_power_spec',true);

V = V.';

A = (Phi * diag(omega))/(Phi);

A_res.A = A; A_res.W = Phi; A_res.lambda = omega; A_res.b = b;

%% Unforced DMD Reconstruction
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
title(['Rank-' num2str(dmd_rank) ' Unforced DMD Reconstruction'])
xlim(plot_xlim)

%% Extract Forcing Signal

% As will be evident from the plot, the discovered forcing in this case
% captures the chaotic switching behavior of the true forcing, but also
% contains additional spurious high frequency oscillations. This often
% occurs as a result of the ambiguous definition of "forcing" as all
% dynamical content which the DMD model fails to reproduce. Better forcing
% results can often be obtained using Optimized DMD instead of Exact DMD.

% Choose which forward time stepping routine to use
% linear_step = @euler_step;
% linear_step_mat = @euler_stepwise_update;

% linear_step = @exp_step;
% linear_step_mat = @exp_stepwise_update;

linear_step = @ode45_step;
linear_step_mat = @ode45_stepwise_update;

V_trunc = V(1:dmd_rank,:);
V_trunc_update = linear_step_mat(V_trunc(:,1:end-1),A_res,dt);
V_trunc_update_augment = [V_trunc_update; zeros(size(V)-[dmd_rank 1])]; % pad rows above r with zeros

V_residual = V(:,2:end) - V_trunc_update_augment;
V_residual = real(V_residual); % cut off any machine-error imaginary components
H_residual = U * S * V_residual;
X_residual = H_residual((res_delay_step-1)*n+(1:n),:);

figure('Name','True vs. Discovered Forcing')
subplot(2,3,1:2)
plot(t,u(:,2))
xlim([0 700])
title('True Forcing')
subplot(2,3,4:5)
plot(t(1:size(X_residual,2)),X_residual(2,:))
title('Discovered Forcing')
xlim([0 700])

% Compute Autocorrelation
u_acf = xcorr(u(:,2));
u_acf = u_acf(floor(length(u_acf)/2):end);
xr_acf = xcorr(X_residual(2,:));
xr_acf = xr_acf(floor(length(xr_acf)/2):end);

subplot(2,3,3)
acf_plot_length = 300;
plot(dt*(0:acf_plot_length),u_acf(1:acf_plot_length+1))
title('Autocorrelation')

subplot(2,3,6)
plot(dt*(0:acf_plot_length),xr_acf(1:acf_plot_length+1))
title('Autocorrelation')

%% Simulate Forced Linear System
V_recon_forced = zeros(size(V));
V_recon_forced(:,1) = V(:,1);
for j = 2:size(V,2)
    this_V = V_recon_forced(:,j-1);
    this_V_trunc = this_V(1:dmd_rank);
    V_recon_forced(1:dmd_rank,j) = linear_step(this_V_trunc,A_res,dt);
    V_recon_forced(:,j) = V_recon_forced(:,j) + V_residual(:,j-1);
end

H_recon_forced = U * S * V_recon_forced;
y_recon_forced = H_recon_forced((res_delay_step-1)*n+(1:n),:);

figure('Name','Forced DMD Reconstruction')
subplot(2,1,1)
plot(t,y,'LineWidth',1.5)
title('Original Signal')
xlim(plot_xlim)
subplot(2,1,2)
plot(t(1:size(y_recon_forced,2)),y_recon_forced,'LineWidth',1.5)
title('Forced DMD Reconstruction')
xlim(plot_xlim)

%% Different forward time stepping schemes
% Euler time step
function v = euler_step(V,A_res,dt)
    A = A_res.A;
    v = (eye(size(A))+dt*A) * V(:,end);
end
function V_update = euler_stepwise_update(V,A_res,dt)
    A = A_res.A;
    dV = A * V;
    V_update = V + dt*dV;
end

% Exponential time step
function v = exp_step(V,A_res,dt)
    lambda = A_res.lambda;
    W = A_res.W;

    B = W\V;
    
    b_abs = vecnorm(B,2,2);
    b_ang = mean(angle(B),2);
    b_local = b_abs .* exp(sqrt(-1)*b_ang);
    

    v = W*(exp(lambda*dt).*b_local);
    
end
function V_update = exp_stepwise_update(V,A_res,dt)
    W = A_res.W;
    lambda = A_res.lambda;
    all_b = W\V;
    
    V_update = W*(repmat(exp(lambda*dt),1,size(V,2)).*all_b);
end

% ode45 time step
function v = ode45_step(V,A_res,dt)
    A = A_res.A;
    v0 = V(:,end);
    [~,vp] = ode45(@(t,v) A*v, [0 dt], v0);
    vp = vp.';
    v = vp(:,end);
end
function V_update = ode45_stepwise_update(V,A_res,dt)
    A = A_res.A;
    sim_duration = size(V,2);
    V_update = zeros(size(V));
    for j = 1:sim_duration
        V_update(:,j) = ode45_step(V(:,j),A_res,dt);
    end
end
