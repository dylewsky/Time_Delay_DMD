function [Phi,omega,b,varargout] = time_delay_dmd(y,t,nDelay_in,delaySteps_in,rank_in,varargin)
%TIME_DELAY_DMD Dynamic mode decomposition on a time-delay embedding of the
%data in y. Data is delay embedded and then projected onto Principal
%Component Trajectory coordinates. The DMD model is regressed in this space
%
%   REQUIRED INPUT ARGUMENTS:
%       y: Input data matrix (state dimension x time dimension)
%       t: Vector of sampling times (evenly spaced)
%       nDelay_in: Number of embeddings used to construct Hankel matrix (or
%           a vector of such numbers to iterate over)
%       delaySteps_in: Number of time steps per delay embedding (or
%           a vector of such numbers to iterate over)
%       rank_in: Truncation rank for the SVD on the Hankel matrix (or
%           a vector of ranks to iterate over)
%
%   OPTIONAL NAME/VALUE PAIR ARGUMENTS:
%       'dmd_type' | 'exact' (default), 'exact_discrete', 'opt', or 'window'
%            Optimized DMD requires the optDMD package available here:
%               https://github.com/duqbo/optdmd
%            Windowed DMD requires the oDMD package available here:
%               https://github.com/haozhg/odmd
%       'dmd_domain' | 'global' (default) or 'windowed'. Global computes a
%           single DMD over the full time domain. Windowed computes iterative
%           DMDs on sliding-window subsets of the data and returns a cell
%           matrix of all results. Tunable windowing parameters are 'window_size'
%           and 'window_step'
%       'window_size' | Width of DMD domain window in time steps (Only applies
%           if using Windowed DMD (or windowed domain on another DMD method))
%       'window_step' | Number of time steps to slide the DMD domain window
%           forward from one iteration to the next(Only applies if using
%           dmd_domain == 'windowed' or dmd_type == 'window')
%       'optdmd_[param] | Parameters to be passed to optimized DMD (see
%           optdmd documentation for details)
%       'eig_init' | [optdmd only] Manually supplied initial guess for
%           optimization eigenvalues. Must be a vector of length less than or
%           equal to the DMD rank. If it is less, only the first r values
%           will be used.
%       'fourier_initialize' | [optdmd only] Use an initial guess for
%           optimization derived from a peak-finding routine on the Fourier
%           power spectrum (Default: True)
%       'fourier_peak_method' | 'sum' (default) or 'separate'. 'sum' runs a
%           single peak-finding routine on the summed spectra of all nVars
%           signal elements. 'separate' runs on each separately, and
%           compiles the top peaks based on prominence (runs the risk of
%           redundancy when the same frequency is present on multiple
%           variables)
%       'fourier_peak_subtract_median' | When doing peak-finding on the
%          Fourier spectrum, subtract a moving median from the power first
%          (Default: True)
%       'optdmd_constrain' | [optdmd only] Option for constraining the
%           search space for DMD eigenvalues. 0 (default) leaves it unconstrained, 1
%           forces values to lie in the left half-plane, and 2 forces
%           values to lie along the imaginary axis (real part = 0)
%       'windowdmd_weighting' | [windowed DMD only] Weighting factor for
%           state sample history. Must be in (0,1]. If 1, all samples are
%           weighted equally. Otherwise, samples from p steps ago are
%           weighted by (weighting factor)^p. (Default: 0.9)
%       'plot_power_spec' | Toggle plot Fourier power spectrum with DMD
%           eigenfrequencies overlaid (Default: false)
%       'plot_A_mat' | Toggle plot visualization of A matrix obtained by DMD (i.e.
%           the linear operator which advances V forward in time) (Default:
%           false)
%       'save_results' | Pass a filename string to save DMD results to
%           disk (if left empty, no file will be saved)
%
%   OUTPUT ARGUMENTS
%       Phi, omega, b: DMD eigenvectors, eigenvalues, and weight
%           coefficients
%       Optional additional arguments: U, S, V, the results of SVD on the
%          embedding matrix H
%       Note that if the function is supplied with multiple sets of
%           embedding parameters over which to iterate, all of these outputs
%           will be cell arrays containing a separate result for each iteration

    p = inputParser;
    p.PartialMatching = 0;
    
    
    addRequired(p,'y',@isnumeric)
    addRequired(p,'t',@isnumeric)
    addRequired(p,'nDelay_in',@isnumeric)
    addRequired(p,'delaySteps_in',@isnumeric)
    addRequired(p,'rank_in',@isnumeric)
    
    dmd_type_default = 'opt';
    dmd_type_valid = {'exact','exact_discrete','opt','window'};
    dmd_type_check = @(x) any(validatestring(x,dmd_type_valid));
    addParameter(p,'dmd_type',dmd_type_default,dmd_type_check);
    
    dmd_domain_default = 'global';
    dmd_domain_valid = {'global','windowed'};
    dmd_domain_check = @(x) any(validatestring(x,dmd_domain_valid));
    addParameter(p,'dmd_domain',dmd_domain_default,dmd_domain_check);
    
    window_size_default = round(length(t)/10);
    addParameter(p,'window_size',window_size_default,@isnumeric);
    
    window_step_default = round(length(t)/40);
    addParameter(p,'window_step',window_step_default,@isnumeric);
    
    optdmd_maxiter_default = 50;
    addParameter(p,'optdmd_maxiter',optdmd_maxiter_default,@isnumeric);

    optdmd_tol_default = 1e-4;
    addParameter(p,'optdmd_tol',optdmd_tol_default,@isnumeric);
    
    optdmd_eps_stall_default = 1e-7;
    addParameter(p,'optdmd_eps_stall',optdmd_eps_stall_default,@isnumeric);
    
    optdmd_imode_default = 2;
    addParameter(p,'optdmd_imode',optdmd_imode_default,@isnumeric);
    
    optdmd_ifprint_default = 0;
    addParameter(p,'optdmd_ifprint',optdmd_ifprint_default,@isnumeric);

    eig_init_default = [];
    addParameter(p,'eig_init',eig_init_default,@isnumeric);
    
    fourier_initialize_default = true;
    addParameter(p,'fourier_initialize',fourier_initialize_default,@islogical);
    
    fourier_peak_method_default = 'sum';
    addParameter(p,'fourier_peak_method',fourier_peak_method_default,@ischar);

    optdmd_constrain_default = 0;
    addParameter(p,'optdmd_constrain',optdmd_constrain_default,@isnumeric);
    
    windowdmd_weighting_default = 0.9;
    addParameter(p,'windowdmd_weighting',windowdmd_weighting_default,@isnumeric);
    
    fourier_peak_subtract_median_default = true;
    addParameter(p,'fourier_peak_subtract_median',fourier_peak_subtract_median_default,@islogical);
 
    
    addParameter(p,'save_results','',@ischar);
    addParameter(p,'plot_power_spec',false,@islogical);
    addParameter(p,'plot_A_mat',false,@islogical);
    
    parse(p,y,t,nDelay_in,delaySteps_in,rank_in,varargin{:});
    
    % Assign contents of p.Results to variables
    for j = 1:length(p.Parameters)
        eval(sprintf('%s = p.Results.%s;',p.Parameters{j},p.Parameters{j}));
    end
    
    % Check time spacing
    if range(diff(t)) > 10^-5
        error('Expected evenly-spaced time vector')
    end
    dt = t(2)-t(1);
    
    % Resolve windowed DMD argument ambiguities
    if strcmp('dmd_type','window') && strcmp('dmd_domain','global')
        warning('dmd_type == ''window'' is not valid for dmd_domain == ''global''. Setting dmd_domain = ''windowed''')
        dmd_domain = 'windowed';
    end
    if strcmp('dmd_type','window') && window_step ~= window_step_default
        warning('Cannot set window step interval for dmd_type == ''window''. Defaults to 1.')
    end
    
    %% Parse time delay parameters
    delay_param_lengths = [length(nDelay_in), length(delaySteps_in), length(rank_in)];
    num_iter = max(delay_param_lengths);
    if max(delay_param_lengths) == min(delay_param_lengths)
        % If all are single values or if all are equal-length lists
        nDelay_list = nDelay_in;
        delayStep_list = delaySteps_in;
        rank_list = rank_in;
        if max(delay_param_lengths == 1)
            varParam = 'None';
        else
            varParam = 'All';
        end
    elseif nnz(delay_param_lengths == 1) == 2
        % If all but one have length 1
        if length(nDelay_in) > 1
            nDelay_list = nDelay_in;
            varParam = 'nDelay';
        else
            nDelay_list = nDelay_in*ones(1,max(delay_param_lengths));
        end
        if length(delaySteps_in) > 1
            delayStep_list = delaySteps_in;
            varParam = 'delaySteps';
        else
            delayStep_list = delaySteps_in*ones(1,max(delay_param_lengths));
        end
        if length(rank_in) > 1
            varParam = 'rank';
            rank_list = rank_in;
        else
            rank_list = rank_in*ones(1,max(delay_param_lengths));
        end
    else
        error('Error: input arguments for nDelay, delaySteps, and rank must be either 1) all of equal length, or 2) all length 1 except for one')
    end
    
    
    %% Compute power spectrum of data
    N = length(t);
    Fs = 1/dt; %sampling freq
    freq = 0:Fs/N:Fs/2;
    
    y_power = zeros(length(freq),size(y,1));
    for j = 1:size(y,1)
        ydft = fft(y(j,1:N).');
        ydft = ydft(1:floor(N/2)+1);
        psdy = (1/(Fs*N)) * abs(ydft).^2;
        psdy(2:end-1) = 2*psdy(2:end-1);
        y_power(:,j) = psdy;
    end
    freq_ang = 2*pi*freq;
    if strcmp(fourier_peak_method,'sum')
        y_power = sum(y_power,2);
    end
    
    all_peak_locs = [];
    all_peak_prominence = [];
    all_peak_origin = [];
    for fj = 1:size(y_power,2)
        if fourier_peak_subtract_median
            [~,peak_locs,~,peak_prominence,y_power_med] = moving_median_findpeaks(10*log10(y_power(:,fj)));
        else
            [~,peak_locs,~,peak_prominence] = findpeaks(10*log10(y_power(:,fj)));
        end
        all_peak_locs = [all_peak_locs; peak_locs];
        all_peak_prominence = [all_peak_prominence; peak_prominence];
        all_peak_origin = [all_peak_origin; fj*ones(size(peak_locs))];
    end
    
    [~,p_sort_ind] = sort(all_peak_prominence,'descend');
    peak_locs = all_peak_locs(p_sort_ind);
    peak_origin = all_peak_origin(p_sort_ind);
    [peak_locs,iup,~] = unique(peak_locs,'stable');
    peak_origin = peak_origin(iup);
    clear('all_peak_locs','all_peak_prominence','peak_prominence','iup');
    

    all_U = cell(num_iter,1);
    all_S = cell(num_iter,1);
    all_V = cell(num_iter,1);
    all_S_sum = zeros(1,num_iter);
    all_Phi = cell(num_iter,1);
    all_omega = cell(num_iter,1);
    all_b = cell(num_iter,1);
    all_A = cell(num_iter,1);
    all_window_bounds = cell(num_iter,1);

    for nd = 1:num_iter
        if num_iter > 1
            disp(['Running time delay DMD iteration ' num2str(nd)])
        end
        nDelay = nDelay_list(nd);
        delaySteps = delayStep_list(nd);
        r = rank_list(nd);
        [H, t_td] = time_delay_embed(y,t,nDelay,delaySteps);
        [U,S,V] = svd(H,'econ');
        all_U{nd} = U;
        all_V{nd} = V;
        all_S{nd} = diag(S);
        all_S_sum(nd) = sum(diag(S));
        
        V_trunc = V(:,1:r);
        
        if strcmp(dmd_type,'exact')
            if strcmp(dmd_domain,'windowed')
                % Note that in this case Phi, omega, b are cell arrays
                % containing all windowed results
                [Phi,omega,b,A] = exactdmd_windowed(V_trunc.',t_td,window_size,window_step);
            elseif strcmp(dmd_domain,'global')
                % COMPUTE DERIVATIVES (4TH ORDER CENTRAL DIFFERENCE)
                dV = zeros(length(V_trunc)-5,r);
                for i=3:length(V_trunc)-3
                    for k=1:r
                        dV(i-2,k) = (1/(12 * dt)) * (-V_trunc(i+2,k)+8 * V_trunc(i+1,k)-8 * V_trunc(i-1,k)+V_trunc(i-2,k));
                    end
                end
                % trim first and last two that are lost in derivative
                Vtrim = V_trunc(3:end-3,1:r);

                % BUILD HAVOK REGRESSION MODEL ON TIME DELAY COORDINATES
                A = Vtrim\dV;

                [Phi,omega] = eig(A);
                omega = diag(omega);
                b = Phi\(Vtrim(1,:).');
            end
        elseif strcmp(dmd_type,'exact_discrete')
            if strcmp(dmd_domain,'windowed')
                warning('windowed domain not implemented for dmd_type == ''exact_discrete''')
            elseif strcmp(dmd_domain,'global')
                V1 = V_trunc(1:end-1,:).';
                V2 = V_trunc(2:end,:).';
                
                A = V2/V1;
                
                [Phi,omega] = eig(A);
                omega = diag(omega);
                b = Phi\V_trunc(1,:).';
            end
        elseif strcmp(dmd_type,'opt')
            % Set optdmd options
            lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
            ubc = [Inf*ones(r,1); Inf*ones(r,1)];
            if optdmd_constrain == 1
                lbc = [-Inf*ones(r,1); -Inf*ones(r,1)];
                ubc = [zeros(r,1); Inf*ones(r,1)];
            elseif optdmd_constrain == 2
                lbc = [-1e-6 * ones(r,1); -Inf*ones(r,1)];
                ubc = [1e-6 * ones(r,1); Inf*ones(r,1)];
            end

            copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
            
            if ~isempty(eig_init) % if eig_init was supplied manually
                if size(eig_init,1) ~= 1 % transpose to row vector, if needed
                    eig_init = eig_init.';
                end
                if length(eig_init) < r
                    warning('Length of eig_init is less than DMD rank. DMD optimization will not be manually initialized.')
                    e_init = [];
                else
                    e_init = eig_init(1:r).';
                    if length(eig_init) > r
                        warning('Length of eig_init is greater than DMD rank. Only first r values will be used')
                    end
                end
                if fourier_initialize
                    warning('optdmd will not be initialized with Fourier spectrum peaks because an alternate initialization has been manually provided');
                end
            else
                if fourier_initialize
                    try
                        e_init = sqrt(-1)*freq_ang(peak_locs(1:floor(r/2)));
                        e_init = repelem(e_init,2);
                        e_init(2:2:end) = -e_init(2:2:end);

                        if length(e_init) ~= r
                            e_init = [0, e_init]; % if rank is odd, expect one pure-real eigenvalue
                        end
                        e_init = e_init.';
                    catch ME
                        warning('Not enough Fourier peaks found to initialize optDMD')
                        e_init = [];
                    end
                else
                    e_init = [];
                end
            end
            opts = varpro_opts('maxiter',optdmd_maxiter,'tol',optdmd_tol,'eps_stall',optdmd_eps_stall,'ifprint',optdmd_ifprint);
            imode = optdmd_imode;
            if strcmp(dmd_domain,'windowed')
                % Note that in this case Phi, omega, b are cell arrays
                % containing all windowed results
                [Phi,omega,b] = optdmd_windowed(V_trunc.',t_td,r,imode,window_size,window_step,opts,e_init,[],copts);
                A = cell(size(Phi));
                for kk = 1:length(Phi)
                    A{kk} = (Phi{kk} * diag(omega{kk}))/(Phi{kk});
                end
                window_step_bounds = zeros(length(Phi),2);
                window_step_bounds(:,1) = 1:window_step:(window_step*length(Phi));
                window_step_bounds(:,2) = window_step_bounds(:,1) + window_size;
                all_window_bounds{nd} = window_step_bounds;
            elseif strcmp(dmd_domain,'global')
                [Phi,omega,b] = optdmd(V_trunc.',t_td,r,imode,opts,e_init,[],copts);
                A = (Phi * diag(omega))/(Phi);
                A = real(A); %remove negligible imaginary components
            end
            
        elseif strcmp(dmd_type,'window')
            % Note that here window_step is used as a downsampling factor,
            % so the window covers a time duration of
            % window_step*window_size*dt
            
            window_size_steps = floor(window_size/window_step);
            
            V_trunc_downsample = V_trunc(1:window_step:end,:);
            V_trunc_x = V_trunc_downsample(1:(end-window_size_steps),:).';
            V_trunc_y = V_trunc_downsample((window_size_steps+1):end,:).';
            
            n_window_steps = size(V_trunc_x,2)-window_size_steps;

            A = cell(n_window_steps,1);
            Phi = cell(n_window_steps,1);
            omega = cell(n_window_steps,1);
            b = cell(n_window_steps,1);
            window_step_bounds = zeros(n_window_steps,2);
            window_step_bounds(1,1) = 1;
            
            wdmd = WindowDMD(r,window_size_steps,windowdmd_weighting);
            wdmd.initialize(V_trunc_x(:,1:window_size_steps), V_trunc_y(:,1:window_size_steps));
            % window DMD
            A{1} = wdmd.A;
            [this_Phi,this_lambda] = eig(wdmd.A);
            Phi{1} = this_Phi;
            omega{1} = log(diag(this_lambda))/(dt*window_step); % continuous-time e'vals
            b{1} = this_Phi\V_trunc_x(:,window_size_steps);

            for k = 2:n_window_steps
                window_step_bounds(k,1) = k*window_step;
                wdmd.update(V_trunc_x(:,window_size_steps+k-1), V_trunc_y(:,window_size_steps+k-1));
                A{k} = wdmd.A;
                [this_Phi,this_lambda] = eig(wdmd.A);
                Phi{k} = this_Phi;
                omega{k} = log(diag(this_lambda))/(dt*window_step); % continuous-time e'vals
                b{k} = this_Phi\V_trunc_x(:,window_size_steps+k-1);
            end
            window_step_bounds(:,2) = window_step_bounds(:,1) + window_size;
            all_window_bounds{nd} = window_step_bounds;
        end
        
        all_Phi{nd} = Phi;
        all_omega{nd} = omega;
        all_b{nd} = b;
        all_A{nd} = A;
        
        if plot_power_spec && strcmp(dmd_domain,'global')
            fig_powerspec = figure();
            p_power = plot(freq_ang,10*log10(y_power),'LineWidth',1.5,'DisplayName','Fourier Power Spec.');
            hold on
            
            pos_omega = omega(imag(omega) > 0);
            if ~isempty(pos_omega)
                % Skip frequency plotting if DMD has failed to identify any
                % conjugate frequency pairs
                p_eig = cell(length(pos_omega),1);
                for j = 1:length(pos_omega)
                    p_eig{j} = plot(imag(pos_omega(j))*[1 1],ylim,'b-','DisplayName','DMD Eigenfreqs');
                    hold on
                end
                try
                    legend([p_power, p_eig{1}])
                catch ME
                    legend([p_power(1), p_eig{1}])
                end
                if strcmp(dmd_type,'opt') && ~isempty(e_init)
                    p_init = cell(floor(r/2),1);
                    if length(eig_init) >= r % if eig_init was supplied manually
                        pos_eig_init = eig_init(imag(eig_init)>=0);
                        pos_eig_init = imag(pos_eig_init(1:floor(r/2)));
                        for j = 1:floor(r/2)
                            [~,init_j_loc] = min(abs(freq_ang - pos_eig_init(j)));
                            p_init{j} = plot(freq_ang(init_j_loc),10*log10(y_power(init_j_loc,peak_origin(j))),'ro','DisplayName','Initial Guess Freqs');
                            hold on
                        end
                        clear('init_j_loc');
                    elseif fourier_initialize
                        for j = 1:floor(r/2)
                            p_init{j} = plot(freq_ang(peak_locs(j)),10*log10(y_power(peak_locs(j),peak_origin(j))),'ro','DisplayName','Initial Guess Freqs');
                            hold on
                        end
                    end
                    
                    try
                        legend([p_power,p_eig{1},p_init{1}])
                    catch ME
                        legend([p_power(1),p_eig{1},p_init{1}])
                    end
                end
            end
            if fourier_peak_subtract_median
                plot(freq_ang,y_power_med,'g','LineWidth',2,'DisplayName','Subtracted Fourier median');
            end
            xlabel('\omega')
            ylabel('Power/Frequency')
            title(['Fourier Power Spectrum'])
        end
        
        if plot_A_mat
            fig_A_mat = figure();
            imagesc(A)
            title(['A Matrix (nDelay = ' num2str(nDelay) ', delaySteps = ' num2str(delaySteps) ', r = ' num2str(r)])
            colorbar
            axis equal
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
        end    
        
    end
    
    %% Prepare output
    if num_iter ~= 1
        U = all_U;
        S = all_S;
        V = all_V;
        Phi = all_Phi;
        omega = all_omega;
        b = all_b;
        A = all_A;
    end
    window_bounds = all_window_bounds;
    if (nargout > 3)
        varargout{1} = U;
    end
    if (nargout > 4)
        varargout{2} = S;
    end
    if (nargout > 5)
        varargout{3} = V;
    end
    
    rank = rank_list;
    nDelay = nDelay_list;
    delaySteps = delayStep_list;

    if ~isempty(save_results)
        if strcmp('dmd_domain','global')
            save(save_results,'U','S','Phi','omega','b','A','nDelay','delaySteps','rank','dmd_type','dmd_domain','-v7.3');
        else
            save(save_results,'U','S','Phi','omega','b','A','nDelay','delaySteps','rank','dmd_type','dmd_domain','window_size','window_step','window_bounds','-v7.3');
        end
    end
end

