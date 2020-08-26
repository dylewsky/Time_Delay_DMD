function dy = vdp_forced_slowosc_rhs(y,t,A,omega,varargin)
    if isempty(varargin)
        Mu = 1;
    elseif length(varargin) == 1
        Mu = varargin{1};
    else
        disp('invalid input argument')
        return;
    end
    
    dy = [y(2); Mu*(1-y(1)^2)*y(2)-y(1)+A*sin(omega*t)];
end