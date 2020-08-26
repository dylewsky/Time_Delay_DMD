function [H,t_td] = time_delay_embed(X,t,nDelay,delaySteps)
% [H,t_td] = time_delay_embed(X,t,nDelay,delaySteps)
% Given time series data X and time vector t, embeds data nDelay times
% offset by delaySteps steps. Returns delay embedding Hankel matrix H and
% corresponding truncated time vector t_td
    if size(X,1) == length(t) && size(X,2) ~= length(t)
        X = X.'; %transpose X if necessary
    elseif size(X,1) ~= length(t) && size(X,2) ~= length(t)
        disp('X and t do not have matching dimensions')
        return;
    end
    nVars = size(X,1);
    
    t_td = t(1:length(t)-(nDelay-1)*delaySteps);
    stackmax = nVars*nDelay; % Number of shift-stacked rows
    
    H = zeros(stackmax,size(X,2)-(nDelay-1)*delaySteps);
    for k=1:nDelay
        delayInds = ((k-1)*nVars + 1) : k*nVars;
        H(delayInds,:) = X(:,(k-1)*delaySteps+1:(size(X,2)-(nDelay)*delaySteps + k*delaySteps));
    end
end

