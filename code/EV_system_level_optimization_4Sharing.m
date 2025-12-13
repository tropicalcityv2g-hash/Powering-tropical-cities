%% ============================================================
%  EV SYSTEM-LEVEL OPTIMIZATION PIPELINE
%  Script with structured local functions
% ============================================================



%% ============================================================
%  MODEL SUMMARY (System-level EV Charging/Discharging QP)
%
%  Purpose
%  - Optimize aggregate EV charging/discharging over a day to smooth the
%    system net demand (non-EV demand minus PV), while keeping each EV
%    SOC feasible w.r.t. driving needs and charging availability.
%  - This shared version uses DUMMY input data for DEMONSTRATION because the
%    full mobility dataset cannot be publicly released. Example/demo inputs
%    are provided under: "data/SystemLoad"
%  - EV energy consumption is modeled based on driving distance; While in this
%    shared demo it is simplified using a constant speed and per-step
%    distance because detailed trajectories are not provided.
%
%  Decision variables (per time step t and vehicle i)
%  - p_ch(t,i)   >= 0 : charging power
%  - p_dis(t,i)  >= 0 : discharging power
%  - SOC(t,i) in [0,1]: state of charge (defined for t = 1..T+1) %% this
%  can be changed to your own [lb, ub], e.g., [0.2, 0.8] to better protect
%  the battery



%  Objective 
%  - Minimize the squared optimized net demand over time:
%      minimize  sum_t  ( D(t) + sum_i( p_ch(t,i) - p_dis(t,i) ) )^2
%    where D(t) = sum(SystemDemand(t,:) - PV(t,:)) is the baseline net load.
%
%  Key constraints
%  - Availability bounds:
%      0 <= p_ch(t,i)  <= p_max * avail(t,i)
%      0 <= p_dis(t,i) <= p_max * avail(t,i)
%  - SOC dynamics with driving consumption and (dis)charging efficiencies:
%      SOC(t+1,i) = SOC(t,i) + (eta_ch*dt/C)*p_ch(t,i)
%                           - (dt/(eta_dis*C))*p_dis(t,i)
%                           - SOC_drop_from_driving(t,i)
%      SOC(1,i) = SOC_init(i),   0 <= SOC(t,i) <= 1
%  - Feasibility guards:
%      SOC at departure times >= SOC_req(t,i)  (computed by backward DP)
%      SOC(T+1,i) >= SOC_init(i)
%
%  Parameter customization
%  - Time resolution, battery capacity, charging power, efficiencies, and
%    number of vehicles/districts can be modified in init_parameters().
%
%% ============================================================


clear; clc;
rng(1);



%% ============================================================
% [1] BUILD & SOLVE OPTIMIZATION
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



%% ============================================================
% [2] BUILD & SOLVE OPTIMIZATION
%% ============================================================
N = size(driving,2);
fprintf("[4] Building and solving optimization model...\n");

result = build_and_solve_EV_SystemOpt(SystemDemand, PV, driving, avail, SOC_init, params);

if strcmp(result.status, 'OPTIMAL')
    fprintf("    Optimization solved successfully\n");
else
    warning("Optimization status: %s", result.status);
end



%% ============================================================
% [3] POST-PROCESSING
%% ============================================================
fprintf("[5] postprocess and plot...\n");
postprocess_and_plot(result, SystemDemand, PV, driving, params);

fprintf("=== Pipeline finished ===\n");










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
end






function result = build_and_solve_EV_SystemOpt(SystemDemand, PV, driving, avail, SOC_init, params)
    % (Your existing Gurobi model construction goes here unchanged)
    % This function should return "result" with fields .x and .status
    %
    % IMPORTANT: this function must not load data or change N
    %
    % Placeholder (replace with your real model):




    cons_kWh_per_km = 0.195;
    v_km_per_h      = 60;
    dt = params.dt;
    c = params.battery_kWh;
    eta_ch = params.eta_ch;
    eta_dis = params.eta_dis;


    [T, N] = size(driving);


    rng(1);
    


    
    D = sum(SystemDemand - PV, 2);    




    departures = cell(N,1);
    for i = 1:N
        dep = find(driving(:,i)==1 & [0; driving(1:end-1,i)]==0);
        departures{i} = dep(:)';
    end

    n_pch  = T*N;
    n_pdis = T*N;
    n_soc  = (T+1)*N;
    nvar   = n_pch + n_pdis + n_soc;
    
    off_pch  = 0;
    off_pdis = n_pch;
    off_soc  = n_pch + n_pdis;


    
    id_pch  = @(t,i) off_pch  + (t-1)*N + i;
    id_pdis = @(t,i) off_pdis + (t-1)*N + i;
    id_soc  = @(t,i) off_soc  + (t-1)*N + i;



    model.modelsense = 'min';
    
    model.vtype = repmat('C', nvar, 1); 

    

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

    model.lb = lb;
    model.ub = ub;
    
    





    dist_per_step_km = v_km_per_h * dt;             
    e_step_kWh       = cons_kWh_per_km * dist_per_step_km; 
    

    SOC_drop = driving .* (e_step_kWh ./ c');  



    row = 0;
    Aeq_i = [];
    Aeq_j = [];
    Aeq_v = [];
    beq = [];
    
    for i=1:N
        row = row + 1;
        Aeq_i(end+1) = row;
        Aeq_j(end+1) = id_soc(1,i);
        Aeq_v(end+1) = 1;
        beq(end+1) = SOC_init(i);
    
        alpha = eta_ch*dt/c;
        beta  = (1/eta_dis)*dt/c;
    



        for t = 1:T
            row = row + 1;
    
            Aeq_i = [Aeq_i, row, row, row, row];
            Aeq_j = [Aeq_j, ...
                     id_soc(t+1,i), ...
                     id_soc(t,i), ...
                     id_pch(t,i), ...
                     id_pdis(t,i)];
            Aeq_v = [Aeq_v, ...
                      1, ...
                     -1, ...
                     -alpha, ...
                      beta];
    
            beq(end+1) = -SOC_drop(t,i);
        end

    end
    
    Aeq = sparse(Aeq_i, Aeq_j, Aeq_v, row, nvar);
    model.A = Aeq;
    model.rhs = beq';
    model.sense = repmat('=', row, 1);

    

    
    %% ------------------- Inequalities: Departure SOC ----------
    A_i = [];
    A_j = [];
    A_v = [];
    b   = [];
    
    [SOC_req, ~] = compute_SOC_req_backward(driving, avail, params);


    r=0;
    for i=1:N
        for t_dep = departures{i} 
            r=r+1;
            A_i(end+1)=r;
            A_j(end+1)=id_soc(t_dep,i);
            A_v(end+1)=-1;      
            b(end+1) = -SOC_req(t_dep,i);
        end
    end
    
    if r>0
        Aineq = sparse(A_i, A_j, A_v, r, nvar);
        model.A  = [model.A; Aineq];
        model.rhs = [model.rhs; b'];
        model.sense = [model.sense; repmat('<', r, 1)];
    end

    r_departure = r;


    


%{  %}

    %% ------------------- Inequalities: End SOC >= Initial SOC -------------------
    A_i2 = [];
    A_j2 = [];
    A_v2 = [];
    b2   = [];
    
    r2 = 0;
    for i=1:N
        r2 = r2 + 1;
        A_i2(end+1)=r2;
        A_j2(end+1)=id_soc(T+1,i);
        A_v2(end+1)=-1;             
        b2(end+1)   = -SOC_init(i); 

    end
    
    Aineq2 = sparse(A_i2, A_j2, A_v2, r2, nvar);
    
    model.A     = [model.A; Aineq2];
    model.rhs   = [model.rhs; b2(:)];
    model.sense = [model.sense; repmat('<', r2, 1)];

    r_end = N;





    nnzA = T * (2*N);
    Ai = zeros(nnzA,1);
    Aj = zeros(nnzA,1);
    Av = zeros(nnzA,1);
    
    k = 0;
    for t = 1:T

        idx_ch  = off_pch  + (t-1)*N + (1:N);
        idx_dis = off_pdis + (t-1)*N + (1:N);
    

        Ai(k+(1:N)) = t;
        Aj(k+(1:N)) = idx_ch;
        Av(k+(1:N)) = 1;
        k = k + N;
    
        Ai(k+(1:N)) = t;
        Aj(k+(1:N)) = idx_dis;
        Av(k+(1:N)) = -1;
        k = k + N;
    end
    
    Aobj = sparse(Ai, Aj, Av, T, nvar);
    

    Q = 2 * (Aobj' * Aobj);      
    c = 2 * (Aobj' * D(:));   
    
    model.Q   = Q;
    model.obj = c;



    params.Method = 2;
    params.OutputFlag = 1;



    
    result = gurobi(model, params);


end





function postprocess_and_plot(result, SystemDemand, PV, driving, params)

    if ~strcmp(result.status, 'OPTIMAL')
        disp("Skipping post-processing due to non-optimal solution");
        return;
    end

    [T, N] = size(driving);
    dt = params.dt;

    time_vec = (0:T-1)' * dt;




    n_pch  = T*N;
    n_pdis = T*N;
    n_soc  = (T+1)*N;    
    
    x = result.x;
    
    p_ch  = reshape(x(1:n_pch),                 [N, T])';
    p_dis = reshape(x(n_pch+1:n_pch+n_pdis),    [N, T])';
    
    soc_start = n_pch + n_pdis + 1;
    soc_end   = n_pch + n_pdis + n_soc;
    
    SOC = reshape(x(soc_start:soc_end), [N, T+1])'; 

    
    P_EV  = sum(p_ch - p_dis, 2);
    D = sum(SystemDemand - PV, 2);
    D_opt = D + P_EV;


    %% ====================== Plot ============================
    
    figure;
    plot(time_vec, sum(SystemDemand, 2), "--", 'LineWidth', 1.5); hold on;
    plot(time_vec, sum(PV, 2), ":", 'LineWidth', 1.5); hold on;
    plot(time_vec, D, "-.", 'LineWidth', 1.5); hold on;
    plot(time_vec, D_opt, 'LineWidth', 1.5);
    xlabel("Time of day (h)");
    ylabel("Net demand (kW)");
    legend("System demand", "PV output", "Baseline net demand", "Optimized  net demand");
    title("System-level EV Charging Optimization (city-wide demand)");
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
    %   Eneed_kWh : [T x N] same in kWh (useful for debugging)

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



