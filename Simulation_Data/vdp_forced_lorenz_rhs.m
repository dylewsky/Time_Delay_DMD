function dy = vdp_forced_lorenz_rhs(y,t,A,varargin)
    if isempty(varargin)
        Mu = 1;
    elseif length(varargin) == 1
        Mu = varargin{1};
    else
        disp('invalid input argument')
        return;
    end
    
    % Rossler
    sigma = 10;
    rho = 28;
    beta = 8/3;
    
    lx = y(3);
    ly = y(4);
    lz = y(5);
    
    lorenz_tau = 10;

    dlx = (1/lorenz_tau)*(sigma*(ly - lx));
    dly = (1/lorenz_tau)*(lx*(rho - lz) - ly);
    dlz = (1/lorenz_tau)*(lx*ly - beta*lz);
    
    
   % VDP
   
   dvy = [y(2); Mu*(1-y(1)^2)*y(2)-y(1)]; % Unforced VDP
   dvy(2) = dvy(2) + A*lx; % Lorenz (x) forcing
   
   dy = [dvy; dlx; dly; dlz];
    
end