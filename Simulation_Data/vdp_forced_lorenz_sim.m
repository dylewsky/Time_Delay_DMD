dt = 0.05;
tspan = 0:dt:1000;


x0 = [2; 0]; % VDP IC
z0 = [-1.0455; -0.5064; 19.5038]; % Lorenz IC
Mu = 1;
A = 0.01; % Forcing amplitude coefficient

y0 = [x0;z0];
[t,y] = ode45(@(t,y) vdp_forced_lorenz_rhs(y,t,A,Mu), tspan, y0);

u = A*[zeros(size(t)) y(:,3)]; % isolated forcing signal

y = y(:,1:2); % cut off Lorenz variables

save('vdp_forced_lorenz_data.mat','t','y','u');




