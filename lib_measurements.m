function [Y_comp, H] = lib_measurements(t, X, station_id, const)
    % LIB_MEASUREMENTS Computes G(X) and H-matrix.
    % Frame: Sun-Centered EMO2000
    
    r_sc = X(1:3);
    v_sc = X(4:6);
    k_SRP = X(7);
    lat4  = X(8);
    lon4  = X(9);
    bias  = X(10);

    
    % 1. Get Earth State (Rotate J2000 -> EMO2000)
    st_earth_J2000 = cspice_spkezr('EARTH', t, 'J2000', 'NONE', 'SUN');
    r_E_Sun_EMO = const.R_EME_EMO * st_earth_J2000(1:3);
    v_E_Sun_EMO = const.R_EME_EMO * st_earth_J2000(4:6);
    
    % 2. Station State (Chain: ECF -> ECI -> EMO)
    % GST relative to J2000
    phi_G = const.phi_G_J2000 + const.we * t;
    c = cos(phi_G); s = sin(phi_G);
    R_ECF_ECI = [c, -s, 0; s, c, 0; 0, 0, 1];
    

    % Station in ECF
    stat = const.stations(station_id);

    % If it's Antarctica, get latitude and longitude from the state.
if station_id == 4
    lat_sta = X(8);   % latitude of station 4 (state)
    lon_sta = X(9);   % longitude of station 4 (state)
else
    % other stations are well known let them come from lib_constants
    lat_sta = stat.lat;
    lon_sta = stat.lon;
end

    cl = cos(lat_sta); 
    sl = sin(lat_sta);
    cd = cos(lon_sta); 
    sd = sin(lon_sta);
    r_sta_ecf = const.R_E * [cl*cd; cl*sd; sl];
    v_sta_ecf = cross([0;0;const.we], r_sta_ecf);


        % Chain Rotation to EMO2000
    % r_emo = R_EME_EMO * R_ECF_ECI * r_ecf
    R_ECF_EMO = const.R_EME_EMO * R_ECF_ECI;
    
    r_sta_emo_earth = R_ECF_EMO * r_sta_ecf;
    v_sta_emo_earth = R_ECF_EMO * v_sta_ecf;
    
    % Shift Origin to Sun
    r_sta = r_E_Sun_EMO + r_sta_emo_earth;
    v_sta = v_E_Sun_EMO + v_sta_emo_earth;
    
    % 3. Compute Measurement
    rho_vec = r_sc - r_sta;
    rho = norm(rho_vec);
    u_rho = rho_vec / rho;
    
    rel_vel = v_sc - v_sta;
    range_rate = dot(rho_vec, rel_vel) / rho;
    
    % Bias Logic (Switch at Day 6)
    t_detect = cspice_str2et('2025 DEC 01 00:00:00.00'); 
    if (t - t_detect) < (6 * 86400)
        y_bias = bias; 
        H_rhodot_bias = 1;
    else
        y_bias = 0; 
        H_rhodot_bias = 0;
    end
    
    Y_comp = [rho; range_rate + y_bias];
    
    %  H Matrix
H_rho_r    = u_rho.';
H_rho_v    = zeros(1,3);
H_rho_srp  = 0;
H_rho_bias = 0;

H_rhodot_r    = (rel_vel' - range_rate * u_rho') / rho;
H_rhodot_v    = u_rho.';
H_rhodot_srp  = 0;

% ---- Sensitivity wrt station-4 latitude and longitude (both range & range-rate) ----
if station_id == 4
    phi4 = X(8);   % lat4 (state)
    lam4 = X(9);   % lon4 (state)

    % 1) Partials in ECEF
    dr_dphi4_ecef = const.R_E * [ -sin(phi4)*cos(lam4);
                                  -sin(phi4)*sin(lam4);
                                   cos(phi4) ];

    dr_dlam4_ecef = const.R_E * [ -cos(phi4)*sin(lam4);
                                   cos(phi4)*cos(lam4);
                                   0 ];

    % 2) Convert to EMO2000
    dr_dphi = R_ECF_EMO * dr_dphi4_ecef;  %dr_dphi is in EMO2000
    dr_dlam = R_ECF_EMO * dr_dlam4_ecef;  %dr_dlam is in EMO2000

    % --- Range partials ---
    H_rho_phi4 = -u_rho.' * dr_dphi;
    H_rho_lam4 = -u_rho.' * dr_dlam;

    % --- Range-rate partials ---
    I3   = eye(3);
    proj = I3 - (u_rho * u_rho.');   % 3x3
    term = proj * rel_vel;           % 3x1

    omega_vec     = [0; 0; const.we];
    dv_dphi4_ecef = cross(omega_vec, dr_dphi4_ecef);
    dv_dlam4_ecef = cross(omega_vec, dr_dlam4_ecef);

    dv_dphi4 = R_ECF_EMO * dv_dphi4_ecef;   % EMO
    dv_dlam4 = R_ECF_EMO * dv_dlam4_ecef;   % EMO

    H_rhodot_phi4 = -(1/rho) * dr_dphi.' * term - u_rho.' * dv_dphi4;
    H_rhodot_lam4 = -(1/rho) * dr_dlam.' * term - u_rho.' * dv_dlam4;

else
    % Measurements from stations 1–3 carry no information about Antarctica's lat/lon
    H_rho_phi4     = 0;
    H_rho_lam4     = 0;
    H_rhodot_phi4  = 0;
    H_rhodot_lam4  = 0;
end


Y_comp = [rho; range_rate + y_bias];

H = [H_rho_r,     H_rho_v ,    H_rho_srp,     H_rho_phi4 ,   H_rho_lam4 ,    H_rho_bias ;
     H_rhodot_r , H_rhodot_v , H_rhodot_srp , H_rhodot_phi4, H_rhodot_lam4 , H_rhodot_bias];


