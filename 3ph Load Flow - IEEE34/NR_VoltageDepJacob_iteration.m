function [Enew,J]=NR_VoltageDepJacob_iteration(Y,E,dPQ0,idx3ph,vdep)
% This function performs one inner iteration of the Newton-Raphson
% algorithm
% Inputs
% - Y       compound admittance matrix in p.u.
% - E       nodal voltage phasor in p.u.
% - dS0     power injection mismatches in p.u.
% - idx3ph      structure containing 3ph indices of all nodes
%    .slack    indices of the slack buses
%    .pq       indices of the PQ buses
%    .pv       indices of the PV buses
% - vdep        structure containing the voltage-dependent weights and
%               exponents
%    .alpha    weights active power voltage-dependency
%    .lambda    exponents active power voltage-dependency
%    .beta     weights reactive power voltage-dependency
%    .omega    exponents reactive power voltage-dependency
% - nph     number of phases

% Auxiliary Variables
Eabs    =   abs(E);
Eangle  =   angle(E);
Yabs    =   abs(Y); 
Yangle  =   angle(Y);
n_nodes =   size(Y,1);

% Jacobian Construction
% Initiliazations
dPdEabs     = zeros(n_nodes,n_nodes); % derivative: P versus E_abs
dPdEangle   = zeros(n_nodes,n_nodes); % derivative: P versus E_arg (theta)
dQdEabs     = zeros(n_nodes,n_nodes); % derivative: Q versus E_abs
dQdEangle   = zeros(n_nodes,n_nodes); % derivative: Q versus E_arg (theta)

for ki=1:n_nodes
    
    % Extra voltage dependent term
    constP = 1/sum((vdep.alpha(ki,:)).*(repmat(abs(E(ki,1)),[1 size(vdep.alpha,2)]).^vdep.lambda(ki,:)),2);
    constQ = 1/sum((vdep.beta(ki,:)).*(repmat(abs(E(ki,1)),[1 size(vdep.beta,2)]).^vdep.omega(ki,:)),2);
    
    dconstP = -(constP^2)*sum((vdep.alpha(ki,:)).*(vdep.lambda(ki,:)).*(repmat(abs(E(ki,1)),[1 size(vdep.alpha,2)]).^(vdep.lambda(ki,:)-1)),2);
    dconstQ = -(constQ^2)*sum((vdep.beta(ki,:)).*(vdep.omega(ki,:)).*(repmat(abs(E(ki,1)),[1 size(vdep.beta,2)]).^(vdep.omega(ki,:)-1)),2);
    
    % Diagonal terms
    dPdEabs(ki,ki)  = constP*2*Yabs(ki,ki)*Eabs(ki)*cos(Yangle(ki,ki)) + ...
                     dconstP*Eabs(ki)^2*Yabs(ki,ki)*cos(Yangle(ki,ki));
    dPdEangle(ki,ki)= 0;
    dQdEabs(ki,ki)  = -constQ*2*Yabs(ki,ki)*Eabs(ki)*sin(Yangle(ki,ki)) - ...
                     dconstQ*Eabs(ki)^2*Yabs(ki,ki)*sin(Yangle(ki,ki));
    dQdEangle(ki,ki)= 0;
    
    for kh=1:n_nodes
        
        theta_ih = Eangle(ki) - Eangle(kh);
        
        if kh~=ki
            dPdEabs(ki,ki)   = dPdEabs(ki,ki) + ...
                constP*Yabs(ki,kh)*Eabs(kh)*cos(theta_ih-Yangle(ki,kh)) + ...
                dconstP*Eabs(ki)*Yabs(ki,kh)*Eabs(kh)*cos(theta_ih-Yangle(ki,kh));
            dPdEabs(ki,kh)   = constP*Yabs(ki,kh)*Eabs(ki)*cos(theta_ih-Yangle(ki,kh));
            
            dPdEangle(ki,ki) = dPdEangle(ki,ki) - ...
                constP*Eabs(ki)*Yabs(ki,kh)*Eabs(kh)*sin(theta_ih - Yangle(ki,kh));
            dPdEangle(ki,kh) = constP*Yabs(ki,kh)*Eabs(ki)*Eabs(kh)*sin(theta_ih - Yangle(ki,kh));
            
            dQdEabs(ki,ki)   = dQdEabs(ki,ki) + ...
                constQ*Yabs(ki,kh)*Eabs(kh)*sin(theta_ih-Yangle(ki,kh)) + ...
                dconstQ*Eabs(ki)*Yabs(ki,kh)*Eabs(kh)*sin(theta_ih-Yangle(ki,kh));
            dQdEabs(ki,kh)   = constQ*Yabs(ki,kh)*Eabs(ki)*sin(theta_ih-Yangle(ki,kh));
            
            dQdEangle(ki,ki) = dQdEangle(ki,ki) + ...
                constQ*Eabs(ki)*Yabs(ki,kh)*Eabs(kh)*cos(theta_ih-Yangle(ki,kh));
            dQdEangle(ki,kh) = -constQ*Yabs(ki,kh)*Eabs(ki)*Eabs(kh)*cos(theta_ih - Yangle(ki,kh));
        end
    end
end

% Remove extra rows (i.e., unnecessary equations)
% slack bus: P & Q, PV buses: Q

dPdEabs(idx3ph.slack,:) = [];
dPdEangle(idx3ph.slack,:) = [];

dQdEabs([idx3ph.pv;idx3ph.slack],:) = [];
dQdEangle([idx3ph.pv;idx3ph.slack],:) = [];

% Remove extra columns (i.e., variables)
% slack bus: E_abs & E_arg, PV nodes: E_abs
    
dPdEabs(:,[idx3ph.pv;idx3ph.slack]) = [];
dQdEabs(:,[idx3ph.pv;idx3ph.slack]) = [];

dPdEangle(:,idx3ph.slack) = [];
dQdEangle(:,idx3ph.slack) = [];

% Combine Jacobian without extra equations/variables
J = [dPdEabs,dPdEangle;dQdEabs,dQdEangle];

% Solution Update
% Solve
dx = J \ dPQ0;

% Reconstruct the solution
dE_abs = zeros(length(Eabs),1);
dE_abs(idx3ph.pq,1) = dx(1:length(idx3ph.pq));

dE_angle = zeros(length(Eangle),1);
dE_angle(sort([idx3ph.pq;idx3ph.pv]),1) = dx((length(idx3ph.pq)+1):end);

% Update
Eabs = Eabs + dE_abs;
Eangle=Eangle + dE_angle;

% Solution
Enew = Eabs.*cos(Eangle)+1i*Eabs.*sin(Eangle); 

end