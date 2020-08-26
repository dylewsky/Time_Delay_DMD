function [pks,varargout] = moving_median_findpeaks(y,varargin)
%MOVING_MEDIAN_PEAKS Function used to find Fourier peaks to use as
%initialization points for Optimized DMD eigenspectrum

    if isempty(varargin)
        window = round(length(y)/10);
    else
        window = varargin{1};
    end
    
    y_med = movmedian(y,window);
%     y_med = movmean(y,window);
    
    [pks,locs,w,p] = findpeaks(y-y_med);
    
    
    if (nargout > 1)
        varargout{1} = locs;
    end
    if (nargout > 2)
        varargout{2} = w;
    end
    if (nargout > 3)
        varargout{3} = p;
    end
    if (nargout > 3)
        varargout{4} = y_med;
    end
    
    
end

