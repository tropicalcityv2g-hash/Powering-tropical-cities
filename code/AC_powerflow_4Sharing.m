%% AC Power flow modeling
%%% This module implements a lightweight AC power flow model for the Singapore electricity network. 
%%% It computes steady-state AC power flows on all lines. 
%%% An adaptive network reinforcement procedure was implemented to ensure feasibility under the 2050 electricity demand and PV generation conditions
%%% The line constraints and transformers constriants are based on the thesis:
%%% [1] Trpovski, A., 2023. Synthetic power grid development: Optimal power system planning-based approach 
%%% and its application for singapore case study (Doctoral dissertation, Technische Universität München).



clear;


%%% --- Load the constructed reinforced network ---
NetworkModel_AC_PowerFlow = load("../data/PowerGrid/NetworkModel_AC_PowerFlow.mat");
NetworkModel_AC_PowerFlow = NetworkModel_AC_PowerFlow.NetworkModel_AC_PowerFlow;


%%% Here we use the Net demand (= system_demand - PV ) data for 2024, March 22 
Load_each_district_table = readmatrix("../data/Data4Figure/AC_PowerFlow_LoadSample.csv");
Load_each_district = Load_each_district_table;

%%%% P_flow is the power flow on all the lines
[Vm, P_flow, S_flow, line_loss, loss_percent, total_loss] = ...
    AC_Load_Flow_WithConstraints(Load_each_district, NetworkModel_AC_PowerFlow);





function [Vm, P_flow, S_flow, line_loss, loss_percent, total_loss] = ...
    AC_Load_Flow_WithConstraints(loads, mpc_final)


    % AC_Load_Flow_WithConstraints (v2)
    % - Applies 'loads' (MW) to 66 kV buses
    
    % Outputs:
    %   Vm           : nb × 1 voltage magnitudes (p.u.)
    %   Va           : nb × 1 voltage angles (deg)
    %   P_flow       : nl × 3 [F_BUS, T_BUS, Pf_ij (MW)]
    %   S_flow       : nl × 3 [F_Bran, T_Bran, Sf_ij (MVA)]
    %   line_loss    : nl × 1 active loss on each branch (MW)
    %   loss_percent : nl × 1 loss as % of |Pf|+|Pt|
    %   total_loss   : scalar (MW)

    % Working copy
    mpc = mpc_final;

    is66   = (mpc.bus(:,10) == 66);
    bus66  = find(is66);
    nb     = size(mpc.bus,1);
    nl     = size(mpc.branch,1);

    % Sanity: we expect 55 loads for 55 66-kV buses
    if numel(loads) ~= numel(bus66)
        error('Length mismatch: loads has %d entries but system has %d 66-kV buses.', ...
              numel(loads), numel(bus66));
    end

    mpc.bus(is66, 3:4) = 0;

    Pd_bus = zeros(nb,1);
    Qd_bus = zeros(nb,1);
    Pd_bus(bus66) = loads(:);                
    Qd_bus(bus66) = 0.35 * loads(:);         
    pv_mask_bus = Pd_bus < 0;       
    load_mask   = is66 & ~pv_mask_bus;
    mpc.bus(load_mask,3) = Pd_bus(load_mask);
    mpc.bus(load_mask,4) = Qd_bus(load_mask);

    pv_bus_list   = find(is66 & pv_mask_bus);
    first_gen_row = size(mpc.gen,1) + 1;
    for b = pv_bus_list.'
        Pg = -Pd_bus(b);  
        newgen = zeros(1, max(size(mpc.gen,2), 21));
        newgen(1:10) = [mpc.bus(b,1), Pg, 0, 9e9, -9e9, 1.00, mpc.baseMVA, 1, max(Pg,1), 0];
        mpc.gen = [mpc.gen; newgen]; 
        if mpc.bus(b,2) == 1, mpc.bus(b,2) = 2; end  
    end
    last_gen_row      = size(mpc.gen,1);

    mpopt   = mpoption('verbose', 0, 'out.all', 0);
    results = runpf(mpc, mpopt);


    %---------------------------
    % Outputs
    %---------------------------
    Vm = results.bus(:,8);               

    Pf = results.branch(:,14);      
    Pt = results.branch(:,16);      
    line_loss = Pf + Pt;            

    denom = abs(Pf) + abs(Pt);
    denom(denom == 0) = NaN;
    loss_percent = line_loss ./ denom * 100;

    total_loss = nansum(line_loss);
    P_flow = [results.branch(:,1), results.branch(:,2), Pf];


    Qf = results.branch(:,15);     
    Qt = results.branch(:,17);      
    
    % Apparent power at each end (MVA)
    Sf = sqrt(Pf.^2 + Qf.^2);       
    St = sqrt(Pt.^2 + Qt.^2);       
    S_flow = [results.branch(:,1), results.branch(:,2), Sf];    % MVA


end


