function a_grad = calculate_a_grad(a_vec,P,theta1,del)
% calculate_a_grad calculates the new axial induction vector
%
% Inputs
% - avec: axial induction factor vector current time step
% - P: power output vector current time step
% - theta1: multiplication coefficient update
% - del: bias coefficient update
% Outputs
% - a_grad: axial induction factor vector next time step

%% Calculate change rate first turbine
m = length(a_vec);
m1 = 3;
sum_g1 = 0;
for idxn = 2:m1 %idxn stands for index_neighbors
    % norm   = sqrt( (P(1)-P(idxn)) * (P(1)-P(idxn)) );
    % g1     = (P(1) - P(idxn) ) * (a_vec(1) - a_vec(idxn))/ norm;
    g1     = sign((P(1) - P(idxn))) * (a_vec(1) - a_vec(idxn));
    sum_g1 = sum_g1 + g1;
end
phi = theta1/(m1-1)*sum_g1 + del;

%% Update axial induction all turbines 
a_grad = nan(size(a_vec));
a_grad(1) = a_vec(1) + phi;
a_grad(1) = min(max(a_grad(1),0),1/2);
for idx = 2:m
        a_grad(idx) = a_vec(idx-1);
end

end

% notes
% in this case, we just assume each turbine has only one strategy
% in fact, each turbine has many strategies
% since a^i_k = [ a^i_1, a^i_2, ..., a^i_n ]'; where n = number of
% strategies
% HW-> build a loop where contains strategies each turbine -> try to insert
% in P_i  so e.g. found a^i_3 is best strategies, P_i(a^i_3) is the max
% energy