function [YY, YYL, YYT, Ampacities, Base] = Ymatrix_3ph(Y_file,nph)
% This function computes the compound admittance matrix in p.u. based on the data
% containted in Y_file (an excel file)
% Notes:
% 1) The code supports multiple basis. The user must carefully input all
% basis in the excel file and be sure that the transformer line paramters
% are in the correct abosulte units!
% 2) If the user wants to output the "non-p.u." matrix he/she should change
% the base values in the excel file to Ab=Vb=1.
%
% Inputs: 
%   - Y_file, Excel file containing line parameters
%   - nph, Number of phases in the network
% Outputs:
%   - YY, Compound Admittance Matrix in p.u.
%   - YYL, A matrix with zero diagonal and only the multi-phase branch admittances on
%         off-diagonal terms in p.u.
%   - YTT: A matrix with zero diagonal and only the shunt admittances on off-diagonal terms 
%          (shunt of PI-circuit of branch corresponding to the node of the
%          row index) in p.u.
%  - Ampacities: A matrix with Branch ampacity limits
%  - Base: Structure containing all base-values needed to get absolute
%          values of all quantities

%% Read Data
% The information of the line parameters are loaded using an excel file 

% Read linedata, amapcity and line lengths
[num,~,~] = xlsread(Y_file,'LineData');

LineData.from       = num(:, 2); % from node
LineData.to         = num(:, 4); % to node
LineData.length     = num(:, 7); % in Miles
LineData.config     = num(:, 9); % configuration number of the branch
LineData.amp        = num(:, 11); % ampacity limits
LineData.basenum_lines    = num(:, 12); % base pertaining to branch

% Read node base numbers
[num,~,~] = xlsread(Y_file,'Original Nodes');
LineData.basenum_nodes    = num(:, 3);

% Read base-values
[num,~,~] = xlsread(Y_file,'Bases');

% Create Basis structure
% Counter
count = 1;

for l = 1:size(num,1)
    % Ab
    Basis{count,1}.Ab = num(l,2);     % Three-Phase Power
    % Vb
    Basis{count,1}.Vb = num(l,3);     % Phase-to-Phase Voltage
    
    % Eb, Ib, Zb & Yb
    Basis{count,1}.Eb = Basis{count,1}.Vb/sqrt(3); % Phase-to-Ground Voltage
    Basis{count,1}.Ib = Basis{count,1}.Ab/(sqrt(3)*Basis{count,1}.Vb); % Nodal Current Base
    Basis{count,1}.Zb = Basis{count,1}.Eb/Basis{count,1}.Ib; % Impedance Base 
    Basis{count,1}.Yb = 1/Basis{count,1}.Zb; % Admittance Base
    
    % Increment counter
    count = count + 1;
end

% Read all the line configurations
[num,~,~] = xlsread(Y_file,'Configurations');

% Create configuration structure with different line parameters
% Counter
count = 1;

for l = 1:nph:size(num,1)
    % Name
    Configurations{count,1}.OrgName = num2str(num(l,1));
    % Longitudinal
    Configurations{count,1}.R = num(l:l+nph - 1, 3:3+nph - 1); % [Ohms]
    Configurations{count,1}.X = num(l:l+nph - 1, 6:6+nph - 1); % [Ohms]
    Configurations{count,1}.Z = (Configurations{count,1}.R + 1i* Configurations{count,1}.X); % [Ohms]
    % Shunt
    Configurations{count,1}.B = 1i*num(l:l+nph - 1, 9:9+nph - 1)*1e-6; % [S]
    % Increment counter
    count = count + 1;
end

n_lines=length(LineData.from);                      % Number of lines
n_nodes=max(max(LineData.from),max(LineData.to));   % Number of nodes

%% Create Admittance Matrix
YY          = cell(n_nodes,n_nodes);
YYL         = cell(n_nodes,n_nodes);
YYT         = cell(n_nodes,n_nodes);
Ampacities  = cell(n_nodes,n_nodes);
Base.Iij    = cell(n_nodes,n_nodes);

% Init Cells
for i = 1:size(YY,1)
    for j = 1:size(YY,2)
        YY{i,j}         = zeros(nph, nph);
        YYL{i,j}        = zeros(nph, nph);
        YYT{i,j}        = zeros(nph, nph);
        Ampacities{i,j} = zeros(nph, 1);
        Base.Iij{i,j}   = zeros(nph, 1);
    end
end

% Create Admittance Matrices
for line = 1:n_lines
    % Get Admittances
    tmp_base_idx = LineData.basenum_lines(line); % Get index of Basis of this branch
    tmp_YL = inv( Configurations{LineData.config(line),1}.Z/Basis{tmp_base_idx,1}.Zb * LineData.length(line) );
    tmp_YT = 0.5*Configurations{LineData.config(line),1}.B/Basis{tmp_base_idx,1}.Yb * LineData.length(line);
    
    % YY
    % off-diagonal
    YY{LineData.from(line), LineData.to(line)} = -tmp_YL;
    YY{LineData.to(line) , LineData.from(line)} = -tmp_YL;
    % diagonal
    YY{LineData.to(line) , LineData.to(line)} = YY{LineData.to(line) , LineData.to(line)} + tmp_YL + tmp_YT;
    YY{LineData.from(line) , LineData.from(line)} = YY{LineData.from(line) , LineData.from(line)} + tmp_YL + tmp_YT;
    
    % YYT
    YYT{LineData.from(line), LineData.to(line)} = tmp_YT;
    YYT{LineData.to(line) , LineData.from(line)} = tmp_YT;    
    
    % YYL
    YYL{LineData.from(line), LineData.to(line)} = tmp_YL;
    YYL{LineData.to(line) , LineData.from(line)} = tmp_YL;  
    
    % Ampacities
    Ampacities{LineData.from(line), LineData.to(line)} = LineData.amp(line)*ones(nph,1)/Basis{tmp_base_idx,1}.Ib;
    Ampacities{LineData.to(line) , LineData.from(line)} = LineData.amp(line)*ones(nph,1)/Basis{tmp_base_idx,1}.Ib;  
    
    % Branch current base
    Base.Iij{LineData.from(line), LineData.to(line)} = Basis{tmp_base_idx,1}.Ib*ones(nph,1);
    Base.Iij{LineData.to(line) , LineData.from(line)} = Basis{tmp_base_idx,1}.Ib*ones(nph,1);     
end

%% Transform Cell 2 Matrix
YY          = cell2mat(YY);
YYL         = cell2mat(YYL);
YYT         = cell2mat(YYT);
Ampacities  = cell2mat(Ampacities);     

% Create Base quantities to be re-used
Base.Iij = cell2mat(Base.Iij);

% Read node base numbers
for n = 1:size(LineData.basenum_nodes,1)
    Base.Ab(n,1) = Basis{LineData.basenum_nodes(n,1),1}.Ab;
    Base.Vb(n,1) = Basis{LineData.basenum_nodes(n,1),1}.Vb;
    Base.Eb(n,1) = Basis{LineData.basenum_nodes(n,1),1}.Eb;
    Base.Ib(n,1) = Basis{LineData.basenum_nodes(n,1),1}.Ib;
end

end
