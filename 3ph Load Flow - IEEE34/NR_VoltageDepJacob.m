function [E,S,S0,I,Jacob] = NR_VoltageDepJacob(Y, S_star, E_star, ...
    E_0, idx3ph, vdep, tol, maxIter)
% This function performs the Newton-Raphson algorithm to compute the
% network state. It is numerically more stable when the inputs are in p.u.
% INPUT
% - Y           nodal admittance matrix in p.u.
% - S_star      apparent power setpoints in p.u. (for all buses):
%               i. active power relevant for pv and pq busses
%               ii. reactive power relevant for pq buses
% - E_star      voltage magnitudes setpoints in p.u. (for all buses), relevant for pv buses
% - E_0         initial voltages (phasors) in p.u.
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
% - tol         tolerance for Newton-Raphson convergence criterion
% - n_max       maximum number of iterations
%
% OUTPUT
% - E           Phase-to-ground nodal voltage phasors at convergence
% - S           Nodal apparent power injections at convergence
% - I           Phase-to-ground nodal voltage currents at convergence
% - J           Jacobian at the solution

% Initialization
abs_E = abs(E_0);
abs_E(idx3ph.pv) = E_star(idx3ph.pv);
angle_E = angle(E_0);
E = abs_E.*exp(1i*angle_E);
 
 for k=1:maxIter
     % Compute nodal currents and power injections
     I=Y*E;
     S=(E.'.*I').'; 
     
     % Compute "Controllable" PQ injections ("considering voltage
     % depeendence)
     Pi0 = real(S)./sum((vdep.alpha).*(repmat(abs(E),[1 size(vdep.alpha,2)]).^vdep.lambda),2);
     Qi0 = imag(S)./sum((vdep.beta).*(repmat(abs(E),[1 size(vdep.beta,2)]).^vdep.omega),2);
     S0 = Pi0 + 1i*Qi0;
    
     % Compute Apparent Power Mismatches for whole network
     dS0 =(S_star - S0);
     dP0 = real(dS0);
     dQ0 = imag(dS0);
     % Keep only relevant mismatches
     dP0(idx3ph.slack) = [];
     dQ0(sort([idx3ph.slack;idx3ph.pv])) = [];
     dPQ0 = [dP0;dQ0];
     
     % Convergence check
     if max(abs(dPQ0)) < tol % tolerance
         break
     elseif( k == maxIter ) % Max number of iteration reached without convergence
         disp('NR Algorithm did not converge in specified maximum number of iterations')
         break
     end
     
     % Compute Newton-Raphson step
     [E,Jacob]=NR_VoltageDepJacob_iteration(Y,E,dPQ0,idx3ph,vdep);
     
 end