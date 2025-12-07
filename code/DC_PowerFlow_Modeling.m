%% DC Power flow modeling
%%% This module implements a lightweight DC power flow model for the Singapore electricity network. 
%%% It computes steady-state active power flows on all lines using a linearized
%%% formulation derived from the full AC power flow equations. 
%{

The DC model assumes:

- Flat voltage magnitudes (|V| = 1 p.u.)
- Small angle differences
- Negligible line resistance (R ≪ X), we set R = 0 
- No reactive power modeling

%}


%% Preparation of sampling data to run a DC power flow simulation
%%% Read Singapore planning area shape file
disp('Spatial mapping ...');
shpdir = '../data/mp14-plng-area-no-sea-planning-area/MP14_PLNG_AREA_NO_SEA_PL.shp';
disp(shpdir);
S = shaperead(shpdir);
info = shapeinfo(shpdir);
p1 = info.CoordinateReferenceSystem;


A = readmatrix("../data/PowerGrid/A_Topology_Singapore.csv");
centroids = readmatrix("../data/PowerGrid/centroids_locations.csv");



index_GasGenerator = [];

for i = 1:length(S)
    PlanningAreaName = S(i).PLN_AREA_N;
    PlanningAreaName = upper(strrep(PlanningAreaName, ' ', ''));
    PlanningAreaName = upper(strrep(PlanningAreaName, '-', ''));


    if PlanningAreaName == "TUAS"
        index_GasGenerator(1) = i;
    elseif PlanningAreaName == "SEMBAWANG"
        index_GasGenerator(2) = i;
    elseif PlanningAreaName == "WESTERNISLANDS"
        index_GasGenerator(3) = i;
    end

end

%%% Here we use the Net demand (= system_demand - PV ) data for 2024, March 22 
Load_each_district_table = readtable("../data/Data4Figure/DC_PowerFlow_LoadSample.csv");
Load_each_district = Load_each_district_table{1, :}';
%%%% P_flow is the power flow on all the lines
P_flow = DC_Load_Flow_Adapted(A, centroids, index_GasGenerator, Load_each_district);





function P_flow = DC_Load_Flow_Adapted(adjacency_matrix, centroids, generator_indices, loads)
    
    n = size(centroids, 1); % Number of nodes
    X_unit = 0.15; %%%%% 0.15 from Andrej Trpovski's thesis (NTU & TUM)
    
    D = zeros(n); % Distance matrix
    B = zeros(n); % Susceptance matrix
    for i = 1:n
        for j = i+1:n
            if adjacency_matrix(i, j) ==1
                distance = haversine(centroids(i,:), centroids(j,:));  
                D(i,j) = distance;
                D(j,i) = distance;
                B(i,j) = -1 / (X_unit * distance); 
                B(j,i) = B(i,j);
            end
        end
    end
    

    for i = 1:n
        B(i,i) = -sum(B(i, ~isinf(B(i,:))));
    end
   
    loads(generator_indices) = -abs(loads(generator_indices));
    
    B_reduced = B(2:end, 2:end); 
    loads_reduced = loads(2:end);
    
    theta = zeros(n,1); 
    theta(2:end) = B_reduced \ loads_reduced; 
    

    P_flow = zeros(n);
    for i = 1:n
        for j = i+1:n
            if D(i,j) > 0 
                P_flow(i,j) = (theta(i) - theta(j)) * B(i,j);
                P_flow(j,i) = -P_flow(i,j);
            end
        end
    end
    
    % Visualization or further processing of P_flow can be done here, with
    % the "ElectricityGridModel_n_Map_4Sharing.m" file
end






function distance = haversine(coord1, coord2)
    % Haversine formula to calculate distances between two points on the Earth
    radius = 6371; % Earth's radius in km
    lat1 = deg2rad(coord1(2));
    lon1 = deg2rad(coord1(1));
    lat2 = deg2rad(coord2(2));
    lon2 = deg2rad(coord2(1));
    deltaLat = lat2 - lat1;
    deltaLon = lon2 - lon1;
    a = sin(deltaLat/2)^2 + cos(lat1) * cos(lat2) * sin(deltaLon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distance = radius * c;
end





