%% ============================================================
%  EV DISTRICT-LEVEL OPTIMIZATION PIPELINE
%  Script with structured local functions
% ============================================================



%% ============================================================
%  MODEL SUMMARY (District-level EV Charging Optimization)
%
%  Purpose
%  - Coordinate EV charging/discharging across multiple districts to
%    mitigate local peak demand, while respecting individual EV
%    mobility, charging availability, and battery constraints.
%  - This shared version uses DUMMY mobility and demand data for
%    DEMONSTRATION only. Full mobility data cannot be publicly released.
%    Example/demo inputs are provided under: "data/SystemLoad"
%  - EV energy consumption is modeled based on driving distance; While in this
%    shared demo it is simplified using a constant speed and per-step
%    distance because detailed trajectories are not provided.
%
%  Decision variables
%  - p_ch(t,i)   >= 0 : charging power of EV i at time t
%  - p_dis(t,i)  >= 0 : discharging power of EV i at time t
%  - SOC(t,i) ∈ [0,1] : state of charge (t = 1..T+1)
%  - z_k         >= 0 : normalized peak net demand proxy for district k
%
%  Objective
%  - Minimize normalized district-level peak demand:
%      minimize  sum_k  z_k / Dmax_k
%    where Dmax_k is the baseline peak demand of district k.
%
%  Key constraints 
%  - Charging/discharging only when available:
%      0 <= p_ch(t,i), p_dis(t,i) <= p_max * avail(t,i)
%  - SOC dynamics with driving energy consumption:
%      SOC(t+1,i) = SOC(t,i)
%                 + (eta_ch*dt/C)*p_ch(t,i)
%                 - (dt/(eta_dis*C))*p_dis(t,i)
%                 - SOC_drop_from_driving(t,i)
%  - Feasibility guards:
%      SOC at departure times >= SOC_req(t,i)  (backward DP)
%      SOC(T+1,i) >= SOC_init(i)
%  - District peak constraints with time-varying EV locations:
%      z_k >= D_k(t) + sum_{i in district k at time t}(p_ch - p_dis)
%
%  Parameter customization
%  - Time resolution, battery capacity, charging power, efficiencies, and
%    number of vehicles/districts can be modified in init_parameters().
%
%% ============================================================




clear; clc;




%% ============================================================
% [1] PARAMETERS & INPUT DATA
%% ============================================================
Num_veh = 200;
params = init_parameters(Num_veh);


driving = readmatrix("../data/MobilityInputSample_dummy/is_driving.csv");
avail = readmatrix("../data/MobilityInputSample_dummy/is_available_charging.csv");
SOC_init =  readmatrix("../data/MobilityInputSample_dummy/SOC_ini.csv");


driving = driving(:, 1:params.target_N);
avail = avail(:, 1:params.target_N);
SOC_init = SOC_init(1:params.target_N);


SystemDemand = readmatrix("../data/MobilityInputSample_dummy/SystemDemand_dummy.csv");
PV = readmatrix("../data/MobilityInputSample_dummy/PV_dummy.csv");


veh2dist = generate_random_veh2dist(96, params.target_N, 55, 42);



%% ============================================================
% [2] BUILD & SOLVE OPTIMIZATION
%% ============================================================
N = size(driving,2);
fprintf("[2] Building and solving optimization model...\n");

result = build_and_solve_EV_DistrictOpt(SystemDemand, PV, driving, avail, SOC_init, veh2dist, params);

if strcmp(result.status, 'OPTIMAL')
    fprintf("    Optimization solved successfully\n");
else
    warning("Optimization status: %s", result.status);
end


%% ============================================================
% [3] POST-PROCESSING
%% ============================================================
fprintf("[3] postprocess and plot...\n");
postprocess_and_plot(result, SystemDemand, PV, driving, veh2dist, params);


fprintf("=== Pipeline finished ===\n");








function veh2dist = generate_random_veh2dist(T, N, K, seed)
    %GENERATE_RANDOM_VEH2DIST
    %   veh2dist(t,i) ∈ {1..K}, IID over time
    %
    % Inputs:
    %   T     number of time steps
    %   N     number of vehicles
    %   K     number of districts
    %   seed  random seed (integer)
    %
    % Output:
    %   veh2dist [T × N]

    rng(seed, 'twister');   % reproducible

    veh2dist = randi(K, T, N);
end







%% ============================================================
%% ====================== LOCAL FUNCTIONS ======================
%% ============================================================

function params = init_parameters(Num_veh)
    params.T  = 96;
    params.dt = 0.25;

    params.battery_kWh = 71;
    params.p_ch_max    = 7.2;
    params.eta_ch      = 0.87;
    params.eta_dis     = 0.90;

    params.energy_kWh_per_km = 0.195;
    params.speed_kmph = 60;

    params.target_N = Num_veh;
    params.K = 55;
end



% ------------------------------------------------------------
function result = build_and_solve_EV_DistrictOpt(SystemDemand, PV, ...
    driving, avail, SOC_init, veh2dist, params)
    %BUILD_AND_SOLVE_EV_DISTRICTOPT
    %   District-level EV charging optimization (time-varying locations).
    %
    %   min_{p_ch, p_dis, SOC, z}  sum_{k=1}^K z_k / Dmax_k
    %
    %   s.t.  SOC dynamics,
    %         charging/discharging bounds,
    %         departure SOC constraints,
    %         end-of-day SOC >= initial SOC,
    %         z_k >= D_k(t) + sum_{i : loc(t,i)=k} [p_ch_i(t)-p_dis_i(t)]  ∀k,t
    %
    %   Inputs
    %   ------
    %   SystemDemand : [T × K] non-EV demand per district
    %   PV           : [T × K] PV generation per district
    %   driving      : [T × N] 1 = driving, 0 = not driving
    %   avail        : [T × N] 1 = available to charge/discharge
    %   SOC_init     : [N × 1] initial SOC (fraction)
    %   veh2dist     : [T × N] integer in {1..K}, district of EV i at time t
    %   params       : struct with fields:
    %                    .dt, .battery_kWh, .p_ch_max, .eta_ch, .eta_dis
    %
    %   Output
    %   ------
    %   result : Gurobi result struct (includes .x, .status, .objval, .x, ...)

    
    [T, N]  = size(driving);
    [T2, K] = size(SystemDemand);
    assert(T2 == T, 'Time dimension mismatch between driving and SystemDemand');
    assert(all(size(PV) == [T K]), 'PV must be T×K');
    assert(length(SOC_init) == N, 'SOC_init must be N×1');

    [Tv, Nv] = size(veh2dist);
    assert(Tv == T && Nv == N, 'veh2dist must be T×N');

    veh2dist = round(veh2dist);              
    assert(all(veh2dist(:) >= 1 & veh2dist(:) <= K), ...
        'veh2dist entries must be in {1..K}');

    dt     = params.dt;
    c_kWh  = params.battery_kWh;
    eta_ch = params.eta_ch;
    eta_dis= params.eta_dis;


    Dk = SystemDemand;  
    Dmax_k = max(SystemDemand, [], 1)';    
    Dmax_k(Dmax_k <= 0) = 1;              

    cons_kWh_per_km = 0.195;
    v_km_per_h      = 60;
    dist_per_step_km = v_km_per_h * dt;
    e_step_kWh       = cons_kWh_per_km * dist_per_step_km;

    SOC_drop = driving .* (e_step_kWh ./ c_kWh);


    departures = cell(N,1);
    for i = 1:N
        dep = find(driving(:,i)==1 & [0; driving(1:end-1,i)]==0);
        departures{i} = dep(:)';
    end


    n_pch  = T * N;
    n_pdis = T * N;
    n_soc  = (T+1) * N;
    n_z    = K;

    nvar   = n_pch + n_pdis + n_soc + n_z;

    off_pch  = 0;
    off_pdis = n_pch;
    off_soc  = n_pch + n_pdis;
    off_z    = n_pch + n_pdis + n_soc;

    id_pch  = @(t,i) off_pch  + (t-1)*N + i;
    id_pdis = @(t,i) off_pdis + (t-1)*N + i;
    id_soc  = @(t,i) off_soc  + (t-1)*N + i;
    id_z    = @(k)   off_z   + k;


    model.modelsense = 'min';
    model.vtype      = repmat('C', nvar, 1);


    lb = zeros(nvar,1);
    ub = inf(nvar,1);


    for t = 1:T
        for i = 1:N
            ub(id_pch(t,i))  = params.p_ch_max * avail(t,i);
            ub(id_pdis(t,i)) = params.p_ch_max * avail(t,i);
        end
    end


    for t = 1:(T+1)
        for i = 1:N
            lb(id_soc(t,i)) = 0;
            ub(id_soc(t,i)) = 1;
        end
    end

    for k = 1:K
        lb(id_z(k)) = 0;
    end

    model.lb = lb;
    model.ub = ub;


    Aeq_i = [];
    Aeq_j = [];
    Aeq_v = [];
    beq   = [];

    row = 0;
    for i = 1:N

        row = row + 1;
        Aeq_i(end+1) = row;
        Aeq_j(end+1) = id_soc(1,i);
        Aeq_v(end+1) = 1;
        beq(end+1)   = SOC_init(i);

        alpha = eta_ch * dt / c_kWh;
        beta  = (1/eta_dis) * dt / c_kWh;


        for t = 1:T
            row = row + 1;

            Aeq_i = [Aeq_i, row, row, row, row];
            Aeq_j = [Aeq_j, ...
                     id_soc(t+1,i), ...
                     id_soc(t,i),   ...
                     id_pch(t,i),   ...
                     id_pdis(t,i)];
            Aeq_v = [Aeq_v, ...
                     1, ...
                    -1, ...
                    -alpha, ...
                     beta];

            beq(end+1) = -SOC_drop(t,i);
        end
    end

    Aeq         = sparse(Aeq_i, Aeq_j, Aeq_v, row, nvar);
    model.A     = Aeq;
    model.rhs   = beq(:);
    model.sense = repmat('=', row, 1);


    Aineq_all = [];
    bineq_all = [];
    sense_all = [];


    [SOC_req, ~] = compute_SOC_req_backward(driving, avail, params);

    A_i = [];
    A_j = [];
    A_v = [];
    b   = [];

    r_dep = 0;
    for i = 1:N
        for t_dep = departures{i}
            r_dep = r_dep + 1;
            A_i(end+1) = r_dep;
            A_j(end+1) = id_soc(t_dep, i);
            A_v(end+1) = -1;                 
            b(end+1)   = -SOC_req(t_dep,i);  
        end
    end

    if r_dep > 0
        Aineq_dep = sparse(A_i, A_j, A_v, r_dep, nvar);
        Aineq_all = [Aineq_all; Aineq_dep];
        bineq_all = [bineq_all; b(:)];
        sense_all = [sense_all; repmat('<', r_dep, 1)];
    end


    A_i2 = [];
    A_j2 = [];
    A_v2 = [];
    b2   = [];

    r_end = 0;
    for i = 1:N
        r_end = r_end + 1;
        A_i2(end+1) = r_end;
        A_j2(end+1) = id_soc(T+1,i);
        A_v2(end+1) = -1;            
        b2(end+1)   = -SOC_init(i);  
    end

    if r_end > 0
        Aineq_end = sparse(A_i2, A_j2, A_v2, r_end, nvar);
        Aineq_all = [Aineq_all; Aineq_end];
        bineq_all = [bineq_all; b2(:)];
        sense_all = [sense_all; repmat('<', r_end, 1)];
    end

    
    A_i3 = [];
    A_j3 = [];
    A_v3 = [];
    b3   = [];
    
    r_peak = 0;
    for t = 1:T
        for k = 1:K
            r_peak = r_peak + 1;
    
            
            A_i3(end+1) = r_peak;
            A_j3(end+1) = id_z(k);
            A_v3(end+1) = -1;
    
            idx_i = find(veh2dist(t,:) == k);
    
            for i = idx_i(:)'
                
                A_i3(end+1) = r_peak;
                A_j3(end+1) = id_pch(t,i);
                A_v3(end+1) =  1;
    
                
                A_i3(end+1) = r_peak;
                A_j3(end+1) = id_pdis(t,i);
                A_v3(end+1) = -1;
            end
    
            b3(end+1) = -Dk(t,k);
        end
    end

    if r_peak > 0
        Aineq_peak = sparse(A_i3, A_j3, A_v3, r_peak, nvar);
        Aineq_all  = [Aineq_all; Aineq_peak];
        bineq_all  = [bineq_all; b3(:)];
        sense_all  = [sense_all; repmat('<', r_peak, 1)];
    end


    if ~isempty(Aineq_all)
        model.A     = [Aeq; Aineq_all];
        model.rhs   = [beq(:); bineq_all(:)];
        model.sense = [repmat('=', size(Aeq,1), 1); sense_all];
    else
        model.A     = Aeq;
        model.rhs   = beq(:);
        model.sense = repmat('=', size(Aeq,1), 1);
    end


    obj = zeros(nvar,1);
    for k = 1:K
        obj(id_z(k)) = 1.0 / Dmax_k(k);
    end
    model.obj = obj;


    gparams.Method     = 2;
    gparams.OutputFlag = 1;

    result = gurobi(model, gparams);

    if ~isfield(result, 'x')
        fprintf('No primal solution returned. Status: %s\n', result.status);
        return;
    end
end







% ------------------------------------------------------------

function postprocess_and_plot(result, SystemDemand, PV, driving, veh2dist, params)


    if ~strcmp(result.status, 'OPTIMAL')
        disp("Skipping post-processing due to non-optimal solution");
        return;
    end

    x = result.x;
    K = params.K;

    [T, N] = size(driving);
    dt = params.dt;

    time_vec = (0:T-1)' * dt;



    n_pch  = T*N;
    n_pdis = T*N;
    n_soc  = (T+1)*N;
    n_z    = params.K;
    
    pch_idx  = 1:n_pch;
    pdis_idx = (n_pch+1):(n_pch+n_pdis);
    soc_idx  = (n_pch+n_pdis+1):(n_pch+n_pdis+n_soc);
    z_idx    = (n_pch+n_pdis+n_soc+1):(n_pch+n_pdis+n_soc+n_z);
    
    p_ch  = reshape(x(pch_idx),  [N,T])';
    p_dis = reshape(x(pdis_idx), [N,T])';
    SOC   = reshape(x(soc_idx),  [N,T+1])';
    z     = x(z_idx);



    
    P_EV  = sum(p_ch - p_dis, 2);
    D = sum(SystemDemand - PV, 2);
    D_opt = D + P_EV;


    %% ====================== Plot ============================

    
    
    P_ev_each = p_ch - p_dis;   
    P_ev_dist = zeros(T, K);
    for t = 1:T
        k_t = veh2dist(t,:); 
        for k = 1:K
            idx = (k_t == k);
            if any(idx)
                P_ev_dist(t,k) = sum(P_ev_each(t, idx));
            end
        end
    end
    
    fig = figure;
    
    for k = 1:K
        Dk_base = SystemDemand(:,k) - PV(:,k);
        Dk_opt  = Dk_base + P_ev_dist(:,k);
    
        ax = subplot(7, 8, k);
        hold(ax, 'on');
    
        h1 = plot(ax, time_vec, SystemDemand(:,k), "--", 'LineWidth', 1.2);
        h2 = plot(ax, time_vec, PV(:,k),          ":",  'LineWidth', 1.2);
        h3 = plot(ax, time_vec, Dk_base,          "-.", 'LineWidth', 1.2);
        h4 = plot(ax, time_vec, Dk_opt,           "-",  'LineWidth', 1.2);

        if k == 1
            hLeg = [h1; h2; h3; h4];
        end
    
        grid(ax, 'on');
        set(ax, 'FontSize', 8);
        ylim([0, 150])
    
    end
    

    han = axes('Position',[0 0 1 1], 'Visible','off');
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    xlabel(han, 'Time of day (h)');
    ylabel(han, 'Power (kW)');
    sgtitle('District-level: Demand and Optimized Net Demand');
    lgd = legend(han, hLeg, ...
        {'System demand', 'PV output', 'Baseline net demand', 'Optimized net demand (district-level)'}, ...
        'Orientation','horizontal', 'Box','off');
    
    lgd.Units = 'normalized';
    lgd.Position = [0.20 0.02 0.60 0.04];   
    lgd.FontSize = 9;

    
    figure;
    plot(time_vec, sum(SystemDemand, 2), "--", 'LineWidth', 1.5); hold on;
    plot(time_vec, sum(PV, 2), ":", 'LineWidth', 1.5); hold on;
    plot(time_vec, D, "-.", 'LineWidth', 1.5); hold on;
    plot(time_vec, D_opt, 'LineWidth', 1.5);
    xlabel("Time of day (h)");
    ylabel("Net demand (kW)");
    legend("System demand", "PV output", "Baseline net demand", "Optimized net demand (district-level)");
    title("District-level EV Charging Optimization (city-wide demand)");
    grid on;
    



end





function [SOC_req, Eneed_kWh] = compute_SOC_req_backward(driving, avail, params)
    %COMPUTE_SOC_REQ_BACKWARD
    % Backward-feasible required SOC given:
    %   driving(t,i) in {0,1}  (1 = driving at time step t)
    %   avail(t,i)   in {0,1}  (1 = can charge at time step t)
    %
    % The recursion computes the MIN energy that must be in the battery at SOC(t, i)
    % so that there exists a future charging plan meeting all future driving needs,
    % under max charging power constraints.
    %
    % OUTPUT:
    %   SOC_req   : [T x N] minimum SOC required at each time step


    [T, N] = size(driving);

    dt      = params.dt;
    C_kWh   = params.battery_kWh;
    eta_ch  = params.eta_ch;


    if isscalar(params.p_ch_max)
        p_ch_max = params.p_ch_max * ones(1, N);
    else
        p_ch_max = reshape(params.p_ch_max, 1, N);
    end


    e_step_kWh = params.energy_kWh_per_km * params.speed_kmph * dt;

    Edrive_kWh  = driving * e_step_kWh;                 
    EchgMax_kWh = avail .* (eta_ch * dt .* p_ch_max);   


    Eneed_kWh = zeros(T, N);
    E_next = zeros(1, N); 

    for t = T:-1:1
        Eneed_now = Edrive_kWh(t,:) + max(E_next - EchgMax_kWh(t,:), 0);
        Eneed_kWh(t,:) = Eneed_now;
        E_next = Eneed_now;
    end

    SOC_req = Eneed_kWh / C_kWh;
    
    SOC_req = min(max(SOC_req, 0), 1);
end




