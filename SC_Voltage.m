function [K, Time] = SC_Voltage(Y,S0,E,idx1ph,idx3ph,idxCtrl,nph,vdep)
% This function computes the nodal voltage sensitivity coefficients for
% each control variable of the nodes listed in idxCtrl
% INPUT
% - Y           nodal admittance matrix
% - S0          apparent power injections at nominal voltage for all buses
% - idx1ph      structure containing network indices of all nodes
%    .slack    indices of the slack buses
%    .pq       indices of the PQ buses
%    .pv       indices of the PV buses
% - idx3ph      structure containing 3ph expanded indices of all nodes
%    .slack    indices of the slack buses
%    .pq       indices of the PQ buses
%    .pv       indices of the PV buses
% - idxCtrl     vector containing netowrk indices of the nodes for which
%               the user wants to compute the voltage sensitivity
%               coefficients (SCs)
% - nph         number of phases in considered grid
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
% - K           A cell structure with complex, magnitude and angle voltage SCs
%               the size is length(idxCtrl) x 4, where
%   - column 1, has the index of the node
%   - column 2 a 2x1 cell containing each the complex voltage SCs
%   pertaining to the control variables of the node. 
%   . If a node is PQ then the 1st cell contains SCs w.r.t. P_{0,l}^{\phi} 
%   (active power) and the 2nd cell contains SCs w.r.t. Q_{0,l}^{\phi} (reactive power). 
%   . If a node is PV then the 1st cell contains SCs w.r.t. P_{0,m}^{\phi} (active power) and the
%   2nd cell contains SCs w.r.t. |\bar{E}^{\phi}_m| (voltage magnitude). 
%   . If a node is slack then the 1st cell contains SCs w.r.t. |\bar{E}^{\phi}_k| (voltage magnitude) and the
%   2nd cell contains SCs w.r.t. \angle(\bar{E}^{\phi}_k) (phase-angle).
%   (!!!) Each cell contains a nph|N| x nph matrix where |N| is the number of nodes (1ph) and
%   each column corresponds to the SCs pertaining to a phase and a control variable as explained above.
%   - column 3&4 have the same structure as column 2 but contain,
%   respectively, the magnitude and phase-angle nodal voltage SCs.

%% Construct A

% Timing
T = tic;

% Create F_{in}^{\phi\phi'}, P0 and Q0
F = diag(E)*conj(Y)*diag(conj(E));
P0 = real(S0);
Q0 = imag(S0);

% Create each sub-matrix
% A11
A11 = zeros(size(idx3ph.pq,1), size(idx3ph.pv,1));
% A12
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pq,idx3ph.pq);
tmp_2 = -P0(idx3ph.pq).*sum((vdep.alpha(idx3ph.pq,:)).* ...
    (vdep.lambda(idx3ph.pq,:)).* ...
    (repmat(abs(E(idx3ph.pq)),[1 size(vdep.alpha,2)]).^(vdep.lambda(idx3ph.pq,:)-1)),2);
tmp_3 = sum(real(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.pq);
A12 = tmp_1 + diag(tmp_2 + tmp_3);
% A13
tmp_1 = imag(F); tmp_1 = tmp_1(idx3ph.pq,idx3ph.pq);
tmp_2 = -sum(imag(F),2); tmp_2 = tmp_2(idx3ph.pq);
A13 = tmp_1 + diag(tmp_2);
% A14
tmp_1 = imag(F); tmp_1 = tmp_1(idx3ph.pq,idx3ph.pv);
A14 = tmp_1;
% A21
A21 = zeros(size(idx3ph.pv,1), size(idx3ph.pv,1));
% A22
tmp_1 = real(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pv,idx3ph.pq);
A22 = tmp_1;
% A23
tmp_1 = imag(F); tmp_1 = tmp_1(idx3ph.pv,idx3ph.pq);
A23 = tmp_1;
% A24
tmp_1 = imag(F); tmp_1 = tmp_1(idx3ph.pv,idx3ph.pv);
tmp_2 = -sum(imag(F),2); tmp_2 = tmp_2(idx3ph.pv);
A24 = tmp_1 + diag(tmp_2);
% A31
A31 = zeros(size(idx3ph.pq,1), size(idx3ph.pv,1));
% A32
tmp_1 = imag(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pq,idx3ph.pq);
tmp_2 = -Q0(idx3ph.pq).*sum((vdep.beta(idx3ph.pq,:)).* ...
    (vdep.omega(idx3ph.pq,:)).* ...
    (repmat(abs(E(idx3ph.pq)),[1 size(vdep.beta,2)]).^(vdep.omega(idx3ph.pq,:)-1)),2);
tmp_3 = sum(imag(F),2)./abs(E); tmp_3 = tmp_3(idx3ph.pq);
A32 = tmp_1 + diag(tmp_2 + tmp_3);
% A33
tmp_1 = -real(F); tmp_1 = tmp_1(idx3ph.pq,idx3ph.pq);
tmp_2 = sum(real(F),2); tmp_2 = tmp_2(idx3ph.pq);
A33 = tmp_1 + diag(tmp_2);
% A34
tmp_1 = -real(F); tmp_1 = tmp_1(idx3ph.pq,idx3ph.pv);
A34 = tmp_1;
% A41
if( ~isempty(idx3ph.pv) )
    tmp_1 = -sum((vdep.beta(idx3ph.pv,:)).* ...
        (repmat(abs(E(idx3ph.pv)),[1 size(vdep.beta,2)]).^(vdep.omega(idx3ph.pv,:))),2);
    A41 = diag(tmp_1);
else
    A41 = [];
end
% A42
tmp_1 = imag(F)./repmat(abs(E).',[size(F,1) 1]); tmp_1 = tmp_1(idx3ph.pv,idx3ph.pq);
A42 = tmp_1;
% A43
tmp_1 = -real(F); tmp_1 = tmp_1(idx3ph.pv,idx3ph.pq);
A43 = tmp_1;
% A44
tmp_1 = -real(F); tmp_1 = tmp_1(idx3ph.pv,idx3ph.pv);
tmp_2 = sum(real(F),2); tmp_2 = tmp_2(idx3ph.pv);
A44 = tmp_1 + diag(tmp_2);

% Assemble A
A = [ A11 A12 A13 A14; ...
      A21 A22 A23 A24; ...
      A31 A32 A33 A34; ...
      A41 A42 A43 A44 ];

Time.A = toc(T);
  
 %% Compute nodal voltage SCs for each control variable

% Initialize Outputs
K = cell(length(idxCtrl),4);
Time.K = [];

 for id_x = 1:length(idxCtrl)
     
    % Initialize cell entry
    K{id_x,1} = idxCtrl(id_x); % 1ph index
    K{id_x,2} = cell(2,1);     % Nodal Voltage SCs (complex)
    K{id_x,3} = cell(2,1);     % Magnitude Voltage SCs
    K{id_x,4} = cell(2,1);     % Angle Voltage SCs
    
    % Identify node-type of idxCtrl [1: slack, 2: pq, 3: pv]
     node_type = 1*sum( idxCtrl(id_x) == idx1ph.slack ) + ...
                 2*sum( idxCtrl(id_x) == idx1ph.pq ) + ...
                 3*sum( idxCtrl(id_x) == idx1ph.pv );    
    for ctrl_var = 1:2
       for ph = 1:nph
     % Time
     T = tic;   
           
     % Initialize outputs
     K{id_x,2}{ctrl_var,1}(:,ph) = zeros(size(E)); % complex
     K{id_x,3}{ctrl_var,1}(:,ph) = zeros(size(E)); % magnitude
     K{id_x,4}{ctrl_var,1}(:,ph) = zeros(size(E)); % angle

     % Get 3ph index
     tmp_idx = (idxCtrl(id_x)-1)*nph + ph;
    
     % Compute magnitude and angle nodal voltage SCs for slack nodes
     switch( node_type )
         case 1 % Slack node
             if(ctrl_var == 1) % if X = |\bar{E}^{\phi}_k|
                %K{id_x,2}{ctrl_var,1}(tmp_idx,nph) = exp(1i*angle(E(tmp_idx,1))); % complex
                K{id_x,3}{ctrl_var,1}(tmp_idx,ph) = 1; % magnitude
             else % if X = \angle(\bar{E}^{\phi}_k)
                %K{id_x,2}{ctrl_var,1}(tmp_idx,nph) = 1i*E(tmp_idx,1); % complex
                K{id_x,4}{ctrl_var,1}(tmp_idx,ph) = 1; % angle
             end
         otherwise % non-slack node
             % Do nothing as derivatives are all zeros.
     end
     
     % Compute magnitude nodal voltage SCs for PV nodes
     switch( node_type )
         case 3 % pv node
             if(ctrl_var == 2) % if X = |\bar{E}^{\phi}_m|
                K{id_x,3}{ctrl_var,1}(tmp_idx,ph) = 1; % magnitude
             else % if X = P0^{\phi}_{0,m}
                % Do nothing as derivative is zero
             end
         otherwise % non-pv node
             % Do nothing as derivatives are all zeros.
     end
     
     % Compute nodal voltage SCs for PQ & PV nodes
     % Construct u(X)
      switch( node_type )
         case 1 % slack node
             if(ctrl_var == 1) % if X = |\bar{E}^{\phi}_k|
                u1 = -((real(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))^2))*real(conj(E(tmp_idx))*exp(1i*angle(E(tmp_idx)))) + ...
                     (imag(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))^2))*imag(conj(E(tmp_idx))*exp(1i*angle(E(tmp_idx)))) );
                u2 = -((real(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))^2))*real(conj(E(tmp_idx))*exp(1i*angle(E(tmp_idx)))) + ...
                     (imag(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))^2))*imag(conj(E(tmp_idx))*exp(1i*angle(E(tmp_idx)))) );
                u3 = (real(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))^2))*imag(conj(E(tmp_idx))*exp(1i*angle(E(tmp_idx)))) - ...
                     (imag(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))^2))*real(conj(E(tmp_idx))*exp(1i*angle(E(tmp_idx))));
                u4 = (real(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))^2))*imag(conj(E(tmp_idx))*exp(1i*angle(E(tmp_idx)))) - ...
                     (imag(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))^2))*real(conj(E(tmp_idx))*exp(1i*angle(E(tmp_idx))));
             else % if X = \angle(\bar{E}^{\phi}_k)
                u1 = -((real(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))^2))*real(conj(E(tmp_idx))*1i*E(tmp_idx)) + ...
                     (imag(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))^2))*imag(conj(E(tmp_idx))*1i*E(tmp_idx)));
                u2 = -((real(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))^2))*real(conj(E(tmp_idx))*1i*E(tmp_idx)) + ...
                     (imag(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))^2))*imag(conj(E(tmp_idx))*1i*E(tmp_idx)));
                u3 = (real(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))^2))*imag(conj(E(tmp_idx))*1i*E(tmp_idx)) - ...
                     (imag(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))^2))*real(conj(E(tmp_idx))*1i*E(tmp_idx));
                u4 = (real(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))^2))*imag(conj(E(tmp_idx))*1i*E(tmp_idx)) - ...
                     (imag(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))^2))*real(conj(E(tmp_idx))*1i*E(tmp_idx));               
             end
          case 2 % pq node
             if(ctrl_var == 1) % if X = P0^{\phi}_{0,l}
                u1 = zeros(length(idx3ph.pq),1);
                u1(tmp_idx == idx3ph.pq,1) = sum((vdep.alpha(tmp_idx,:)).* ...
                    (repmat(abs(E(tmp_idx)),[1 size(vdep.alpha,2)]).^ ...
                    (vdep.lambda(tmp_idx,:))),2);
                u2 = zeros(length(idx3ph.pv),1);
                u3 = zeros(length(idx3ph.pq),1);             
                u4 = zeros(length(idx3ph.pv),1);
             else % if X = Q0^{\phi}_{0,l}
                u1 = zeros(length(idx3ph.pq),1);
                u2 = zeros(length(idx3ph.pv),1);
                u3 = zeros(length(idx3ph.pq),1);
                u3(tmp_idx == idx3ph.pq,1) = ...
                    sum((vdep.beta(tmp_idx,:)).* ...
                    (repmat(abs(E(tmp_idx)),[1 size(vdep.beta,2)]).^ ...
                    (vdep.omega(tmp_idx,:))),2);                
                u4 = zeros(length(idx3ph.pv),1);               
             end   
          case 3 % pv node
             if(ctrl_var == 1) % if X = P0^{\phi}_{0,m}
                u1 = zeros(length(idx3ph.pq),1);
                u2 = zeros(length(idx3ph.pv),1);
                u2(tmp_idx == idx3ph.pv,1) = sum((vdep.alpha(tmp_idx,:)).* ...
                    (repmat(abs(E(tmp_idx)),[1 size(vdep.alpha,2)]).^ ...
                    (vdep.lambda(tmp_idx,:))),2);
                u3 = zeros(length(idx3ph.pq),1);
                u4 = zeros(length(idx3ph.pv),1);
             else % if X = |\bar{E}^{\phi}_m|
                u1 = -(real(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))));
                u2 = -(real(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx))));
                u2(tmp_idx == idx3ph.pv,1) = ...
                    u2(tmp_idx == idx3ph.pv,1) + ...
                    P0(tmp_idx)*sum((vdep.alpha(tmp_idx,:)).* ...
                    (vdep.lambda(tmp_idx,:)).*(repmat(abs(E(tmp_idx)),[1 size(vdep.alpha,2)]).^ ...
                    (vdep.lambda(tmp_idx,:)-1)),2) - ...
                    sum(real(F(tmp_idx,:)),2)/abs(E(tmp_idx));                
                u3 = -(imag(F(idx3ph.pq,tmp_idx))/(abs(E(tmp_idx))));
                u4 = -(imag(F(idx3ph.pv,tmp_idx))/(abs(E(tmp_idx)))); 
                u4(tmp_idx == idx3ph.pv,1) = ...
                    u4(tmp_idx == idx3ph.pv,1)+ ...
                    Q0(tmp_idx)*sum((vdep.beta(tmp_idx,:)).* ...
                    (vdep.omega(tmp_idx,:)).*(repmat(abs(E(tmp_idx)),[1 size(vdep.beta,2)]).^ ...
                    (vdep.omega(tmp_idx,:)-1)),2) - ... 
                     sum(imag(F(tmp_idx,:)),2)/abs(E(tmp_idx));                 
             end   
          otherwise % "error"
              % Do nothing
     end    
     
     % Solve A*x(X)=u(X) 
     x = linsolve(A,[u1;u2;u3;u4]);
     % Assemble magnitude nodal voltage SCs
     K{id_x,3}{ctrl_var,1}(idx3ph.pq,ph) = ...
         x( (length(idx3ph.pv)+1):(length(idx3ph.pv) + ...
         length(idx3ph.pq)) ); % PQ nodes
     % Assemble phase-angle nodal voltage SCs
     K{id_x,4}{ctrl_var,1}(idx3ph.pq,ph) = ...
         x( (length(idx3ph.pv)+length(idx3ph.pq)+1) : ...
            (length(idx3ph.pv)+length(idx3ph.pq) + ...
            length(idx3ph.pq)) ); % PQ nodes     
     K{id_x,4}{ctrl_var,1}(idx3ph.pv,ph) = ...
         x( (length(idx3ph.pv)+2*length(idx3ph.pq)+1) : ...
            (length(idx3ph.pv)+2*length(idx3ph.pq) + ...
            length(idx3ph.pv)) ); % PV nodes          
     % Infer nodal voltage SCs 
     K{id_x,2}{ctrl_var,1}(:,ph) = ...
        E.*( (1./abs(E)).*K{id_x,3}{ctrl_var,1}(:,ph) + ...
            1i*K{id_x,4}{ctrl_var,1}(:,ph) );
     
     % Time
     Time.K = [Time.K; toc(T)];
        
        end
    end
 end

end