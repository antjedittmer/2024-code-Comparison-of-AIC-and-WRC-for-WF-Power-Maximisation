function [PT, P, V ] = calculate_PT(a_vec, c, rho, A, Vinf)
% calculate_PT calculates the total power of m wind farms aligned in wind
% direction
%
% Inputs:
% - a_vec: vector of axial induction factor of turbines
% - c: coeffient matrix wake effects (non-zero elements in subdiagonal)
% - rho: air density
% - A: area wind turbine
% - Vinf: free flow wind speed
% Outputs:
% - Pt: wind farm total power
% - P: vector with wind turbine power
% - V: vecotr with wind speeds at turbines

%% Some optimizer expect vector input
if size(a_vec,1) > size(a_vec,2)
    a_vec = a_vec';
end
a_vec = max(min(a_vec,0.5),0);


%% Calculate power of each turbine based on wind and wakes
m = length(a_vec);
c = c(1:m,1:m);
V = zeros(m,1);
P = zeros(m,1);
for idx = 1:m % loop over turbines
    
    %         % interaction to down stream windturbine
    %         for idxc = 1:idx-1
    %             c(idx,idxc) = ( D/ (D + 2* b * (x(idxc)-x(idx)) ) ) ^2;
    %         end
    
    V(idx) = Vinf * (1 - 2 * sqrt( sum( (c(idx,:).*a_vec).^2 ))); %disp(V)
    Cp     = 4*a_vec(idx)*(1-a_vec(idx)).^2; % power coefficient
    
    P(idx) = 1/2* rho * A * Cp * V(idx)^3;
end

%%
PT = sum(P);
end
