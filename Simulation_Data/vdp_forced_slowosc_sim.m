% Use Matlab built-in to simulate VDP oscillator with mu=10 out to t=200
dt = 0.05;
tspan = 0:dt:1000;
y0 = [2; 0];
Mu = 1;

A = 0.1; % forcing amplitude
omega = 0.08; % forcing frequency

[t,y] = ode45(@(t,y) vdp_forced_slowosc_rhs(y,t,A,omega,Mu), tspan, y0);


u = A*[zeros(size(t)) sin(omega*t)]; % isolated forcing signal


save('vdp_forced_slowosc_data.mat','t','y','u');




