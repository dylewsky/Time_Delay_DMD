% Simulate Van der Pol oscillator

dt = 0.05;
tspan = 0:dt:1000;
y0 = [2; 0];
Mu = 1;

[t,y] = ode45(@(t,y) vdp_rhs(y,Mu), tspan, y0);

save('vdp_data.mat','t','y');





