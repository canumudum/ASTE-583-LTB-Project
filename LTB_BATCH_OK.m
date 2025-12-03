%  Batch filter for LTB Preliminary OD solution

% ---------- 1) Initialization ----------
clear all; clc;

% load constants
const = lib_constants();

% Check if Spice kernels are loaded if not call init_project 
try
    t0_et = cspice_str2et(const.epoch_utc_str);
catch
    init_project();   % Irfan's code call
    t0_et = cspice_str2et(const.epoch_utc_str);
end


% A priori Uncertainites For Stations
sigma_rho_dsn = 0.001 ;            % 1m for DSN Noise
sigma_rho_antarctica = 0.01 ;      % 10m for Antarctica Noise
sigma_rho_dot_dsn = 1e-7 ;         % 0.1mm/s for DSN Noise
sigma_rho_dot_antarctica = 1e-6 ;  % 1mm/s for Antarctica Noise

% A priori Values 
k0_ref    = const.k_SRP_0;            % k_srp = 1.0
lat4_0    = const.stations(4).lat;    % rad, nominal -80 deg
lon4_0    = const.stations(4).lon;    % rad, nominal 0 deg
bias0_ref = 0;                        % km/s, a priori mean of Doppler bias

% A priori Covariance 
sigma_r = 100 ;      % 100km spherical (given)
sigma_v = 0.001 ;    % 1m/s spherical (given)
sigma_ksrp = 1/3 ;   % 3-sigma, 100% uncertainity (given)
sigma_lat4 =  0.5 * const.deg2rad;   
sigma_lon4 =  1 * const.deg2rad;   
sigma_bias =  1e-4;  % Try values 1mm/s , 10mm/s , 100m/s , 1000mm/s


% A Priori Covariance Matrix
P0_bar = diag([ sigma_r^2, sigma_r^2, sigma_r^2, ...   % r_x, r_y, r_z  
                sigma_v^2, sigma_v^2, sigma_v^2, ...   % v_x, v_y, v_z
                sigma_ksrp^2, ...                      % k_SRP
                sigma_lat4^2, ...                      % lat4
                sigma_lon4^2, ...                      % lon4
                sigma_bias^2 ]);                       % bias


% A priori reference state (Sun-centered EMO2000)
r0_ref = const.X0_ref(1:3);   % km
v0_ref = const.X0_ref(4:6);   % km/s
Phi0 = eye(10);         % STM 10x10
X0_ref = [ r0_ref;      % 1–3
           v0_ref;      % 4–6
           k0_ref;      % 7
           lat4_0;      % 8
           lon4_0;      % 9
           bias0_ref ]; % 10


% ---------- 2) Load measurements (0-6 days) ----------
Meas = readtable('ASTE583_Project_LTB_Measurements_0-6D_Truth.csv','VariableNamingRule','preserve');

%Meas.Properties.VariableNames = {'Time (s)', 'Station ID', 'Range (km)'};

Detection = 6*24*3600;                     % 6 days after detection [s]
Meas_6 = Meas(Meas.('Time (s)') <= Detection, :);   % Means take measurements within 6 days after detection 
t_meas   = Meas_6.('Time (s)') ;                    % [s] since detection
N_meas   = length(t_meas);                 % Number of the measurements

% ET time at detection epoch
% t0_et = cspice_str2et(const.epoch_utc_str);  % Time (s) times converted to Ephemeris Time (ET)


% ODE options
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);  % Attention : integration and convergence tolerance should be equal

% First guess for initial state in the iteration, No correction at first iteration      
X0_ref_iter = X0_ref ;     % will be updated after each iteration,  DeltaX_hat is the correction

max_iter  = 12;
tolerance = 1e-10;   % Should be the same with ODE45 Solver's tolerance

for iter = 1:max_iter % Main Iterative Batch Filter Loop Start 

    fprintf('\n=== Iteration %d ===\n', iter);

    % ---- Propagate augmented state (state + STM) ----
    % Augmented initial state: [X0_ref_iter ; Phi0(:)]

    X0_aug = [X0_ref_iter; Phi0(:)];          % 10 + 100 = 110 x 1

    tspan_et = [t0_et, t0_et + Detection];     % from detection to detection + 6days

    % lib_dynamics, will be running in STM mode (110-dimension state)
    [t_prop, X_prop] = ode45(@(t,x) lib_dynamics(t, x, const), tspan_et, X0_aug, opts);

    P0_bar_inverse = inv(P0_bar); 
    deltastar      = P0_bar_inverse;    % Information Matrix Delta Star 
    Nstar          = zeros(10,1);       % Information Vector N star , a priori delta state which is delta x naught bar is zero

     % Preallocation , generating empy matrices to be filled.

    rho_hat        = zeros(N_meas,1);               % model range [km] (won't be used in preliminary OD)
    rhodot_hat       = zeros(N_meas,1);             % model range-rate [km/s]
    prefit_residual_range      = zeros(N_meas,1);   % prefit residual range (won't be used in preliminary OD)
    prefit_residual_rangerate  = zeros(N_meas,1);   % prefit residual range-rate


     % ----  Loop over all measurements ----
    for k = 1:N_meas   % Measurement Loop Start

        % Measurement time in ET
        tk_rel = t_meas(k);               % [s] since detection
        tk_et  = t0_et + tk_rel;          % ET

        % Interpolate propagated augmented state at tk_et
        Xk_aug = interp1(t_prop, X_prop, tk_et);   % 1 x 110
        xk     = Xk_aug(1:10).';                   % 10x1 state at tk
        Phi_k  = reshape(Xk_aug(11:end), 10, 10);  % 10x10 STM at tk

        % Station ID
        station_id = Meas_6.('Station ID')(k);

        % Observed Doppler (range-rate), km/s
        y_obs_dot = Meas_6.('Range (km)')(k); % !!! THIS IS RANGE RATE, NOT RANGE THERE IS A TYPO IN MEASUREMENT FILES!

        % Compute modelled measurement and H at tk
        [Y_comp, H] = lib_measurements(tk_et, xk, station_id, const);
        y_comp_dot  = Y_comp(2);          % take only range-rate
        Hi_bar          = H(2,:);         % 1x10, Doppler row

        % Range-Rate Residual ( Observed - Computed )
         delta_y = y_obs_dot - y_comp_dot;   % scalar

        % Store Prefit Residuals for plotting 
         rhodot_hat(k)     = y_comp_dot;                         % model Doppler at tk
         prefit_residual_rangerate(k) = y_obs_dot - y_comp_dot;  % same as delta_y

        % Sensitivity wrt initial state: H_bar_i = H_i * Phi_k
        H_i = Hi_bar * Phi_k;            % 1x10

        % Measurement noise sigma for this station
        if station_id == 4
            sigma_dop = sigma_rho_dot_antarctica;
        else
            sigma_dop = sigma_rho_dot_dsn;
        end
        Ri = sigma_dop^2;               % scalar variance
        Rinv = inv(Ri);                 % can delete this since it's not being used

        % Update information matrix and vector
        deltastar = deltastar + (H_i.' * H_i) / Ri;    % 10x10
        Nstar  = Nstar  + (H_i.' * delta_y)   / Ri;    % 10x1

    end  % End of Measurement Loop

       % Normal Equation computation
       deltaX0_hat = deltastar \ Nstar;   % 10×1

       % State update
       X0_ref_iter = X0_ref_iter + deltaX0_hat;

       % --- Convergence test (Mahalanobis norm) ---
       J = deltaX0_hat.' * (deltastar * deltaX0_hat);   % scalar

       fprintf('  ||deltaX0_hat|| = %.3e    J = %.3e\n\n', norm(deltaX0_hat), J);

tolerance = 1e-10;   % same as ODE45 RelTol/AbsTol ( to not to get fake perturbations)

if J < tolerance
    fprintf('*** Converged at iteration %d ***\n', iter);
    break;
end
 
end % End Of Main Iterative Batch Filter Loop

% ---------- Postfit propagation with final solution ----------
% After we got final estimation for initial state we will propagate the
% trajectory and calculate the residuals.

X0_hat = X0_ref_iter;   % final estimate of initial state
X0_aug_post = [X0_hat; Phi0(:)];   % can run with or without STM
tspan_et = [t0_et, t0_et + Detection];  % time span for postfit propagation
[t_prop_post, X_prop_post] = ode45(@(t,x) lib_dynamics(t, x, const), tspan_et, X0_aug_post, opts); % propagation

% A Posteriori State Covariance

P0_hat = inv(deltastar);            % Covariance we get from the batch filter 
Phi_6   = reshape(X_prop_post(end,11:end),10,10);   %  extract STM at 6 days
P_6_days  = Phi_6 * P0_hat * Phi_6.';    % postfit ( a posteriori) state covariance

% ---------- Postfit residuals (same style as prefit) ----------
rhodot_hat_post              = zeros(N_meas,1);   % preallocation for model range-rate after estimation 
postfit_residual_rangerate   = zeros(N_meas,1);   % preallocation postfit residuals

for k = 1:N_meas

    % Measurement time (ET)
    tk_rel = t_meas(k);           % [s] since detection
    tk_et  = t0_et + tk_rel;      % ET

    % State at tk from POSTFIT propagation
    Xk_aug_post = interp1(t_prop_post, X_prop_post, tk_et); % interpolating because time steps for propagation and measurements are not the same
    xk_post     = Xk_aug_post(1:10).';  %  gets the state vector from the augmented state vector  % 10x1 

    % Station ID
    station_id = Meas_6.('Station ID')(k);  % get whichever station is getting measurements

    % Observed Doppler
    y_obs_dot = Meas_6.('Range (km)')(k); % get observed range-rate from measurement files
                                          % !!! THIS IS RANGE RATE, NOT RANGE THERE IS A TYPO IN MEASUREMENT FILES!

    % Modelled measurement with FINAL state
    [Y_comp_post, ~] = lib_measurements(tk_et, xk_post, station_id, const); % call the measurement model using the final state
    y_comp_dot_post  = Y_comp_post(2); %  range rate model at time tk

    % Store for plotting (exactly the same format with prefit residuals calculation)
    rhodot_hat_post(k)            = y_comp_dot_post;  
    postfit_residual_rangerate(k) = y_obs_dot - y_comp_dot_post;  % residual = observation - computed 

end




