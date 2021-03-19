% main.m -> computes admittance matrix, load-flow and voltage SCs. Includes
% validation of the computed voltage SCs. 
% In it's present form the code can run. If other simulation inputs should
% change, users need to update the code lines beneath:
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
% @author  Sherif Fahmy
% @version 1.0
% @since   2020-10-01

%% Clear all
clear; clc; close all;

%% Set default interpreter
set(0,'defaulttextInterpreter','latex');

%% Add to path
addpath(genpath('3ph Load Flow - IEEE34'));
addpath(genpath('LoadFlow input examples'));

%% Simulation configuration
%% Number of Phases
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
LF.nph = 3; 
%% Create IEEE34 Admittance Matrix 
% (!!!) The excel file format needs to be respected if another network is used (!!!)
% The voltage regulators of the IEEE-34 feeder are ommitted
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
Y_file = 'IEEE34Feeder.xlsx';
[LF.YY, LF.YYL, LF.YYT, LF.Ampacities, LF.Basis] = Ymatrix_3ph(Y_file,LF.nph);
LF.NumNodes = size(LF.YY,1)/LF.nph; % Number of Nodes

%% Node Configuration
% Single-phase indexing
% index of the slack bus
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
LF.idx.slack = 1;
% indices of the PV buses
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
LF.idx.pv = sort([4, 9, 17, 25].'); % With PV nodes
% LF.idx.pv = []; % Without PV nodes
% Indices of the PQ buses
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
LF.idx.loads = [5, 12, 14, 18, 22, 24, 29].';         % nodes with loads (negative active injections)
LF.idx.gen   = [7, 15, 20, 27, 30].';                 % nodes with generation (positive active injections)
LF.idx.zero  = (1:LF.NumNodes).';                     % PQ nodes with no injections (zero-injection nodes)
LF.idx.zero([LF.idx.slack; LF.idx.pv; LF.idx.loads; LF.idx.gen]) = [];
LF.idx.pq = sort([LF.idx.loads; LF.idx.gen; LF.idx.zero]);
% Indices of the control nodes
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
LF.idxCtrl = sort([LF.idx.slack; LF.idx.pv; LF.idx.pq]); % Test All nodes

% Three-phase expanded indexing
% indices of the slack buses
LF.idx_3ph.slack = []; 
for k = 1:length(LF.idx.slack)
    % 3ph indices
    a = (LF.idx.slack(k)-1)*LF.nph+1;
    b = (LF.idx.slack(k)-1)*LF.nph+LF.nph;    
    % Create structures
    LF.idx_3ph.slack = [LF.idx_3ph.slack; (a:b)'];
end
% sort 3ph indices of slack buses
LF.idx_3ph.slack = sort(LF.idx_3ph.slack); 

% indices of the pv buses
LF.idx_3ph.pv = [];
for k = 1:length(LF.idx.pv)
    % 3ph indices
    a = (LF.idx.pv(k)-1)*LF.nph+1;
    b = (LF.idx.pv(k)-1)*LF.nph+LF.nph;   
    % Create structures
    LF.idx_3ph.pv = [LF.idx_3ph.pv; (a:b)'];
end
% sort 3ph indices of pv buses
LF.idx_3ph.pv = sort(LF.idx_3ph.pv);

% indices of the pq buses
LF.idx_3ph.pq = []; LF.idx_num_3ph.pq = [];
for k = 1:length(LF.idx.pq)
    % 3ph indices
    a = (LF.idx.pq(k)-1)*LF.nph+1;
    b = (LF.idx.pq(k)-1)*LF.nph+LF.nph;
    % Create structures   
    LF.idx_3ph.pq = [LF.idx_3ph.pq; (a:b)'];
end
% sort 3ph indices of pq buses
LF.idx_3ph.pq = sort(LF.idx_3ph.pq);

%% Newton-Raphson LoadFlow Settings
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
LF.tol = 1e-8;                           % Newton-Raphson tolerance
LF.maxIter = 1000;                       % Newton-Raphson maximum number of iterations 
% Active power voltage-dependence
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
LF.vdep.alpha = [ones(size(LF.YY,1),1), 0*ones(size(LF.YY,1),1), 0*ones(size(LF.YY,1),1)];  
LF.vdep.lambda = [2*ones(size(LF.YY,1),1), 0*ones(size(LF.YY,1),1), 0*ones(size(LF.YY,1),1)]; 
% Reactive power voltage-dependence
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
LF.vdep.beta = [ones(size(LF.YY,1),1), 0*ones(size(LF.YY,1),1), 0*ones(size(LF.YY,1),1)];  
LF.vdep.omega = [2*ones(size(LF.YY,1),1), 0*ones(size(LF.YY,1),1), 0*ones(size(LF.YY,1),1)]; 

% Initialize values of voltage in the network (line-to-line voltage) to 1
% p.u. with correct phases
E_0=[];
for k=1:LF.NumNodes
    for ph=1:LF.nph
        E_0=[E_0;1*exp(-1i*(ph-1)*(360/LF.nph)/180*pi)];
    end
end

% Create Apparent Power injection setpoins vector
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
% Here you have to input a complex vector in the form S = [S_0^a, S_0^b, S_0^c, ..., S_N^a, S_N^b, S_N^c].
% The vector has to have the same size as the network nodes (so including
% the slack node) but you can put zeros in the slack and zero-injection
% nodes and only put the injections in the PQ nodes
S_star = load('S_star_IEEE34_example.mat','S_star').S_star;
% Create Voltage Magnitude setpoint vector
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
% Here you have to input a complex vector in the form E = [E_0^a, E_0^b, ?E_0^c, ..., E_N^a, E_N^b, E_N^c].
% The vector has to have the same size as the network nodes (so including
% the slack node) but you can put zeros in the slack and zero-injection
% nodes and only put the injections in the PQ nodes
if(~isempty(LF.idx.pv))
    E_star = load('E_star_IEEE34_example_withPQnodes.mat','E_star').E_star;
else
    E_star = load('E_star_IEEE34_example_noPQnodes.mat','E_star').E_star;
end

%% Solve Load-Flow
% This computes the nodal voltages, nodal currents, apparent powers and
% Jacobian
[E_LF,S_LF,S0_LF,I,J] = NR_VoltageDepJacob(LF.YY, S_star, E_star, E_0, LF.idx_3ph, LF.vdep, LF.tol, LF.maxIter); 

%% Compute SCs
% Computes complex, magnitude and angle nodal voltage SCs for all indices in LF.idxCtrl, 
% for all phases and for all automatically identified control variables
[K, Time] = SC_Voltage(LF.YY,S0_LF,E_LF,LF.idx,LF.idx_3ph,LF.idxCtrl,LF.nph,LF.vdep);

%% Validate Sensitivity coefficients
%% Validate PQ nodes mag and angle SCs w.r.t. P0^{\phi}_{0,l} and Q0^{\phi}_{0,l} and PV nodes mag and angle SCs w.r.t. P0^{\phi}_{0,m}
% Reference SCs (Benchmark 1 - Jacobian inversion)
% Invert Jacobian
J_inv = inv(J);
% Extract needed SCs 
K_pq_jacob.mag.P = J_inv(1:length(LF.idx_3ph.pq), 1:(length(LF.idx_3ph.pq)+length(LF.idx_3ph.pv)));
K_pq_jacob.mag.Q = J_inv(1:length(LF.idx_3ph.pq), (length(LF.idx_3ph.pq)+length(LF.idx_3ph.pv))+1:end);
K_pq_jacob.angle.P = J_inv(length(LF.idx_3ph.pq)+1:end, 1:(length(LF.idx_3ph.pq)+length(LF.idx_3ph.pv)));
K_pq_jacob.angle.Q = J_inv(length(LF.idx_3ph.pq)+1:end, (length(LF.idx_3ph.pq)+length(LF.idx_3ph.pv))+1:end);

% Restructure our method's SCs to compare with benchmarck 1 (K_pq_jacob)
K_pq_method.mag.P = []; K_pq_method.angle.P = [];
K_pq_method.mag.Q = []; K_pq_method.angle.Q = [];

for k = 1:size(K,1)
    if( sum( K{k,1} == LF.idx.pq ))
        K_pq_method.mag.P = [K_pq_method.mag.P, K{k,3}{1,1}];
        K_pq_method.mag.Q = [K_pq_method.mag.Q, K{k,3}{2,1}];
        K_pq_method.angle.P = [K_pq_method.angle.P, K{k,4}{1,1}];
        K_pq_method.angle.Q = [K_pq_method.angle.Q, K{k,4}{2,1}];
    elseif( sum( K{k,1} == LF.idx.pv ) )
        K_pq_method.mag.P = [K_pq_method.mag.P, K{k,3}{1,1}];
        K_pq_method.angle.P = [K_pq_method.angle.P, K{k,4}{1,1}];      
    end
end
% Remove zero-derivatives of slack and pv voltages not given by reduced Jacobian
K_pq_method.mag.P([LF.idx_3ph.slack; LF.idx_3ph.pv],:) = []; 
K_pq_method.mag.Q([LF.idx_3ph.slack; LF.idx_3ph.pv],:) = [];
K_pq_method.angle.P(LF.idx_3ph.slack,:) = [];
K_pq_method.angle.Q(LF.idx_3ph.slack,:) = [];

% Computeand Display maximum and RMS error
% Maximum
disp('Benchmarking of nodal voltage sensitivity coefficients (Benchmark 1 - Jacobian Inversion)')
disp(['Maximum error for Magnitude SCs w.r.t P0^{\phi}_{0,l} and w.r.t P0^{\phi}_{0,m} is: ' num2str(max(max(abs( K_pq_method.mag.P - K_pq_jacob.mag.P ))))]);
disp(['Maximum error for Magnitude SCs w.r.t Q0^{\phi}_{0,l} is: ' num2str(max(max(abs( K_pq_method.mag.Q - K_pq_jacob.mag.Q ))))]);
disp(['Maximum error for Angle SCs w.r.t P0^{\phi}_{0,l} and w.r.t P0^{\phi}_{0,m} is: ' num2str(max(max(abs( K_pq_method.angle.P - K_pq_jacob.angle.P ))))]);
disp(['Maximum error for Angle SCs w.r.t Q0^{\phi}_{0,l} is: ' num2str(max(max(abs( K_pq_method.angle.Q - K_pq_jacob.angle.Q ))))]);
% RMS
disp(['RMS error for Magnitude SCs w.r.t P0^{\phi}_{0,l} and w.r.t P0^{\phi}_{0,m} is: ' num2str(sqrt(mean(mean( (K_pq_method.mag.P - K_pq_jacob.mag.P).^2 ))))]);
disp(['RMS error for Magnitude SCs w.r.t Q0^{\phi}_{0,l} is: ' num2str(sqrt(mean(mean( (K_pq_method.mag.Q - K_pq_jacob.mag.Q).^2 ))))]);
disp(['RMS error for Angle SCs for w.r.t P0^{\phi}_{0,l} and w.r.t P0^{\phi}_{0,m} is: ' num2str(sqrt(mean(mean( (K_pq_method.angle.P - K_pq_jacob.angle.P).^2 ))))]);
disp(['RMS error for Angle SCs for w.r.t Q0^{\phi}_{0,l} is: ' num2str(sqrt(mean(mean( (K_pq_method.angle.Q - K_pq_jacob.angle.Q).^2 ))))]);
disp(' ')

%% Validate slack nodes mag and angle SCs w.r.t. |\bar{E}^{\phi}_k| and \angle(\bar{E}^{\phi}_k)
% Reference SCs (Benchmark 2 - delta variations method)
disp('Benchmarking of nodal voltage sensitivity coefficients (Benchmark 2 - Delta Variations)')
% Variation
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
delta_E = [0.01, 0.001, 0.0001, 0.00001, 0.000001]; % p.u.
delta_Eang = [1, 0.1, 0.01, 0.001, 0.0001]*pi/180; % rad

% Initiliaze structures to store max and rms values of errors
Max.slack.mag.E = []; Max.slack.mag.Eang = [];
Max.slack.angle.E = []; Max.slack.angle.Eang = [];
RMS.slack.mag.E = []; RMS.slack.mag.Eang = [];
RMS.slack.angle.E = []; RMS.slack.angle.Eang = [];

for delta = 1:length(delta_E)
% Display
disp(['Delta variations are ' num2str(delta_E(delta)) 'p.u. for slack voltage magnitude and ' ...
    num2str(delta_Eang(delta)) 'rad for slack voltage phase-angle'])
    
% Initialize SCs
K_slack_delta.mag.E = [];  K_slack_delta.mag.Eang = [];
K_slack_delta.angle.E = []; K_slack_delta.angle.Eang = [];

for k1 = 1:length(LF.idx.slack)
    for ctr_var = 1:2
        for ph1 = 1:LF.nph   
        % Initialize values of voltage in the network (line-to-line voltage) to 1
        % p.u. with correct phases
        E_0=[];
        for k2=1:LF.NumNodes
            for ph2=1:LF.nph
                E_0=[E_0;1*exp(-1i*(ph2-1)*(360/LF.nph)/180*pi)];
            end
        end
        if( ctr_var == 1 ) % |\bar{E}^{\phi}_k|
            % Create variation of slack magnitude
            E_0((LF.idx.slack(k1)-1)*LF.nph + ph1) = (1+delta_E(delta))*exp(-1i*(ph1-1)*(360/LF.nph)/180*pi);
        else % \angle(\bar{E}^{\phi}_k)
            % Create variation of slack phase-angle
            E_0((LF.idx.slack(k1)-1)*LF.nph + ph1) = 1*exp((-1i*((ph1-1)*(360/LF.nph)/180*pi+delta_Eang(delta))));
        end
        
        % Create setpoints
        % E_star
        if(~isempty(LF.idx.pv))
            E_star = load('E_star_IEEE34_example_withPQnodes.mat','E_star').E_star;
        else
            E_star = load('E_star_IEEE34_example_noPQnodes.mat','E_star').E_star;
        end
        
        % Recompute LF
        [E_LF_2,~,~,~,~] = NR_VoltageDepJacob(LF.YY, S_star, E_star, E_0, LF.idx_3ph, LF.vdep, LF.tol, LF.maxIter); 

        % Compute SCs
        if( ctr_var == 1 ) % |\bar{E}^{\phi}_k|
            % Compute delta slack magnitude
            K_slack_delta.mag.E = [K_slack_delta.mag.E, (abs(E_LF_2) - abs(E_LF))/delta_E(delta)];
            K_slack_delta.angle.E = [K_slack_delta.angle.E, (angle(E_LF_2) - angle(E_LF))/delta_E(delta)];
        else % \angle(\bar{E}^{\phi}_k)
            % Compute delta slack angle
            K_slack_delta.mag.Eang = [K_slack_delta.mag.Eang, (abs(E_LF_2) - abs(E_LF))/(-delta_Eang(delta))];
            K_slack_delta.angle.Eang = [K_slack_delta.angle.Eang, (angle(E_LF_2) - angle(E_LF))/(-delta_Eang(delta))];
        end
        
        end
    end
end

% Restructure our method's SCs to compare with reference
K_slack_method.mag.E = []; K_slack_method.angle.E = [];
K_slack_method.mag.Eang = []; K_slack_method.angle.Eang = [];

for k = 1:size(K,1)
    if( sum( K{k,1} == LF.idx.slack ))
        K_slack_method.mag.E = [K_slack_method.mag.E, K{k,3}{1,1}];
        K_slack_method.mag.Eang = [K_slack_method.mag.Eang, K{k,3}{2,1}];
        K_slack_method.angle.E = [K_slack_method.angle.E, K{k,4}{1,1}];
        K_slack_method.angle.Eang = [K_slack_method.angle.Eang, K{k,4}{2,1}];    
    end
end

% Compute maximum and RMS errors
% Maximum
Max.slack.mag.E = [Max.slack.mag.E; max(max(abs( K_slack_method.mag.E - K_slack_delta.mag.E )))]; 
Max.slack.mag.Eang = [Max.slack.mag.Eang; max(max(abs( K_slack_method.mag.Eang - K_slack_delta.mag.Eang )))];
Max.slack.angle.E = [Max.slack.angle.E; max(max(abs( K_slack_method.angle.E - K_slack_delta.angle.E )))]; 
Max.slack.angle.Eang = [Max.slack.angle.Eang; max(max(abs( K_slack_method.angle.Eang - K_slack_delta.angle.Eang )))];

% RMS
RMS.slack.mag.E = [RMS.slack.mag.E; sqrt(mean(mean( (K_slack_method.mag.E - K_slack_delta.mag.E).^2 )))]; 
RMS.slack.mag.Eang = [RMS.slack.mag.Eang; sqrt(mean(mean( (K_slack_method.mag.Eang - K_slack_delta.mag.Eang).^2 )))];
RMS.slack.angle.E = [RMS.slack.angle.E; sqrt(mean(mean( (K_slack_method.angle.E - K_slack_delta.angle.E).^2 )))]; 
RMS.slack.angle.Eang = [RMS.slack.angle.Eang; sqrt(mean(mean( (K_slack_method.angle.Eang - K_slack_delta.angle.Eang).^2 )))];

% Display maximum and RMS errors
% Maximum
disp(['Maximum error for Magnitude SCs w.r.t |\bar{E}^{\phi}_k| is: ' num2str(max(max(abs( K_slack_method.mag.E - K_slack_delta.mag.E ))))]);
disp(['Maximum error for Magnitude SCs w.r.t \angle(\bar{E}^{\phi}_k) is: ' num2str(max(max(abs( K_slack_method.mag.Eang - K_slack_delta.mag.Eang ))))]);
disp(['Maximum error for Angle SCs w.r.t |\bar{E}^{\phi}_k| is: ' num2str(max(max(abs( K_slack_method.angle.E - K_slack_delta.angle.E ))))]);
disp(['Maximum error for Angle SCs w.r.t \angle(\bar{E}^{\phi}_k) is: ' num2str(max(max(abs( K_slack_method.angle.Eang - K_slack_delta.angle.Eang ))))]);
% RMS
disp(['RMS error for Magnitude SCs w.r.t |\bar{E}^{\phi}_k| is: ' num2str(sqrt(mean(mean( (K_slack_method.mag.E - K_slack_delta.mag.E).^2 ))))]);
disp(['RMS error for Magnitude SCs w.r.t \angle(\bar{E}^{\phi}_k) is: ' num2str(sqrt(mean(mean( (K_slack_method.mag.Eang - K_slack_delta.mag.Eang).^2 ))))]);
disp(['RMS error for Angle SCs for w.r.t |\bar{E}^{\phi}_k| is: ' num2str(sqrt(mean(mean( (K_slack_method.angle.E - K_slack_delta.angle.E).^2 ))))]);
disp(['RMS error for Angle SCs for w.r.t \angle(\bar{E}^{\phi}_k) is: ' num2str(sqrt(mean(mean( (K_slack_method.angle.Eang - K_slack_delta.angle.Eang).^2 ))))]);

end

%% Validate PV nodes mag and angle SCs w.r.t. |\bar{E}^{\phi}_m|
% Reference SCs (delta method)
% Variation
% (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
delta_E = [0.01, 0.001, 0.0001, 0.00001, 0.000001]; % p.u.

% Initiliaze structures to store max and rms values of errors
Max.pv.mag.E = []; 
Max.pv.angle.E = []; 
RMS.pv.mag.E = []; 
RMS.pv.angle.E = [];

for delta = 1:length(delta_E)
% Display
disp(['Delta variations are ' num2str(delta_E(delta)) 'p.u. for pv voltage magnitude'])

% Initialize SCs
K_pv_delta.mag.E = [];  K_pv_delta.angle.E = [];
% Initialize values of voltage in the network (line-to-line voltage) to 1
% p.u. with correct phases
E_0=[];
for k2=1:LF.NumNodes
    for ph2=1:LF.nph
        E_0=[E_0;1*exp(-1i*(ph2-1)*(360/LF.nph)/180*pi)];
    end
end

for k1 = 1:length(LF.idx.pv)
    for ph1 = 1:LF.nph   

        % Create variation of pv magnitude |\bar{E}^{\phi}_m|
        % (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!) INPUT (!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
        if(~isempty(LF.idx.pv))
            E_star = load('E_star_IEEE34_example_withPQnodes.mat','E_star').E_star;
        else
            E_star = load('E_star_IEEE34_example_noPQnodes.mat','E_star').E_star;
        end
        E_star((LF.idx.pv(k1)-1)*LF.nph + ph1) = E_star((LF.idx.pv(k1)-1)*LF.nph + ph1)+delta_E(delta);

        % Recompute LF
        [E_LF_2,~,~,~,~] = NR_VoltageDepJacob(LF.YY, S_star, E_star, E_0, LF.idx_3ph, LF.vdep, LF.tol, LF.maxIter); 

        % Compute SCs |\bar{E}^{\phi}_m|
        K_pv_delta.mag.E = [K_pv_delta.mag.E, (abs(E_LF_2) - abs(E_LF))/delta_E(delta)];
        K_pv_delta.angle.E = [K_pv_delta.angle.E, (angle(E_LF_2) - angle(E_LF))/delta_E(delta)];
    end
end
       
% Restructure our method's SCs to compare with reference
K_pv_method.mag.E = []; K_pv_method.angle.E = [];

for k = 1:size(K,1)
    if( sum( K{k,1} == LF.idx.pv ))
        K_pv_method.mag.E = [K_pv_method.mag.E, K{k,3}{2,1}];
        K_pv_method.angle.E = [K_pv_method.angle.E, K{k,4}{2,1}];  
    end
end

% Compute maximum and RMS error
% Maximum
Max.pv.mag.E = [Max.pv.mag.E; max(max(abs( K_pv_method.mag.E - K_pv_delta.mag.E )))]; 
Max.pv.angle.E = [Max.pv.angle.E; max(max(abs( K_pv_method.angle.E - K_pv_delta.angle.E )))]; 

% RMS
RMS.pv.mag.E = [RMS.pv.mag.E; sqrt(mean(mean( (K_pv_method.mag.E - K_pv_delta.mag.E).^2 )))]; 
RMS.pv.angle.E = [RMS.pv.angle.E; sqrt(mean(mean( (K_pv_method.angle.E - K_pv_delta.angle.E).^2 )))];

% Display maximum and RMS errors
% Maximum
disp(['Maximum error for Magnitude SCs w.r.t |\bar{E}^{\phi}_m| is: ' num2str(max(max(abs( K_pv_method.mag.E - K_pv_delta.mag.E ))))]);
disp(['Maximum error for Angle SCs for w.r.t |\bar{E}^{\phi}_m| is: ' num2str(max(max(abs( K_pv_method.angle.E - K_pv_delta.angle.E ))))]);
% RMS
disp(['RMS error for Magnitude SCs w.r.t |\bar{E}^{\phi}_m| is: ' num2str(sqrt(mean(mean( (K_pv_method.mag.E - K_pv_delta.mag.E).^2 ))))]);
disp(['RMS error for Angle SCs for w.r.t |\bar{E}^{\phi}_m| is: ' num2str(sqrt(mean(mean( (K_pv_method.angle.E - K_pv_delta.angle.E).^2 ))))]);

end

%% Vizualize errors 
% Difference between delta-variations method and proposed method for SCs
% w.r.t |\bar{E}^{\phi}_m|, |\bar{E}^{\phi}_k| and \angle(\bar{E}^{\phi}_k)
figure; 
subplot(2,2,1) % Max w.r.t pv and slack magnitudes
semilogx(delta_E, Max.slack.mag.E, 'LineWidth', 2.5, 'DisplayName', '$\partial |\bar{E}^{\phi}_i|$ / $\partial |\bar{E}^{\Phi}_k|$'); hold on;
semilogx(delta_E, Max.slack.angle.E, 'LineWidth', 2.5, 'DisplayName', '$\partial \angle (\bar{E}^{\phi}_i)$ / $\partial |\bar{E}^{\Phi}_k|$');
semilogx(delta_E, Max.pv.mag.E, 'LineWidth', 2.5, 'DisplayName', '$\partial |\bar{E}^{\phi}_i|$ / $\partial |\bar{E}^{\Phi}_m|$');
semilogx(delta_E, Max.pv.angle.E, 'LineWidth', 2.5, 'DisplayName', '$\partial \angle (\bar{E}^{\phi}_i)$ / $\partial |\bar{E}^{\Phi}_m|$'); hold off;
l = legend(); set(l,'Interpreter','latex'); xlabel('$\Delta |\bar{E}^{\Phi}_k|$ or $\Delta |\bar{E}^{\Phi}_m|$ [p.u.]'); grid on; axis tight;
title('Max of the error difference'); set(gca,'fontsize', 25);

subplot(2,2,2) % Max w.r.t slack phase-angle
semilogx(delta_Eang, Max.slack.mag.Eang, 'LineWidth', 2.5, 'DisplayName', '$\partial |\bar{E}^{\phi}_i|$ / $\partial \angle (\bar{E}^{\Phi}_k)$'); hold on;
semilogx(delta_Eang, Max.slack.angle.Eang, 'LineWidth', 2.5, 'DisplayName', '$\partial \angle (\bar{E}^{\phi}_i)$ / $\partial \angle (\bar{E}^{\Phi}_k)$');
l = legend(); set(l,'Interpreter','latex'); xlabel('$\Delta \angle (\bar{E}^{\Phi}_k)$ [rad]'); grid on; axis tight;
title('Max of the error difference'); set(gca,'fontsize', 25);

subplot(2,2,3) % RMS w.r.t pv and slack magnitudes
semilogx(delta_E, RMS.slack.mag.E, 'LineWidth', 2.5, 'DisplayName', '$\partial |\bar{E}^{\phi}_i|$ / $\partial |\bar{E}^{\Phi}_k|$'); hold on;
semilogx(delta_E, RMS.slack.angle.E, 'LineWidth', 2.5, 'DisplayName', '$\partial \angle (\bar{E}^{\phi}_i)$ / $\partial |\bar{E}^{\Phi}_k|$');
semilogx(delta_E, RMS.pv.mag.E, 'LineWidth', 2.5, 'DisplayName', '$\partial |\bar{E}^{\phi}_i|$ / $\partial |\bar{E}^{\Phi}_m|$');
semilogx(delta_E, RMS.pv.angle.E, 'LineWidth', 2.5, 'DisplayName', '$\partial \angle (\bar{E}^{\phi}_i)$ / $\partial |\bar{E}^{\Phi}_m|$'); hold off;
l = legend(); set(l,'Interpreter','latex'); xlabel('$\Delta |\bar{E}^{\Phi}_k|$ or $\Delta |\bar{E}^{\Phi}_m|$ [p.u.]'); grid on; axis tight;
title('RMS of the error difference'); set(gca,'fontsize', 25);

subplot(2,2,4) % RMS w.r.t slack phase-angle
semilogx(delta_Eang, RMS.slack.mag.Eang, 'LineWidth', 2.5, 'DisplayName', '$\partial |\bar{E}^{\phi}_i|$ / $\partial \angle (\bar{E}^{\Phi}_k)$'); hold on;
semilogx(delta_Eang, RMS.slack.angle.Eang, 'LineWidth', 2.5, 'DisplayName', '$\partial \angle (\bar{E}^{\phi}_i)$ / $\partial \angle (\bar{E}^{\Phi}_k)$');
l = legend(); set(l,'Interpreter','latex'); xlabel('$\Delta \angle (\bar{E}^{\Phi}_k)$ [rad]'); grid on; axis tight;
title('RMS of the error difference'); set(gca,'fontsize', 25);

%% Timing Results
% Time to compute A
disp(['Time to construct A: ' num2str(Time.A*1e3) 'ms']);

% Time to solve each linear system (A*x(X) = u(X))
figure;
histfit(Time.K*1e3)
xlabel('Computation times [ms]'); title('Histogram Distribution of LSE solving times');

% Total time for proposed method
disp(['Time to construct A + compute ALL ' num2str(length(LF.idxCtrl)*LF.nph*2) ...
    ' (' num2str(length(LF.idxCtrl)) ' control nodes, ' num2str(LF.nph) ' phases and 2 control variables per node) SCs: ' num2str(Time.A*1e3 + sum(Time.K*1e3)) 'ms']);
