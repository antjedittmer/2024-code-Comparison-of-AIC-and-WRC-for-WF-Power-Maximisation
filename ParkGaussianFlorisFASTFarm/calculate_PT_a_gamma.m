function [PT, P, V ] = calculate_PT_a_gamma(ctrlVec, rho, A, Vinf, xPos,b)
% calculate_PT calculates the total power of m wind farms aligned in wind
% direction
%
% Inputs:
% - ctrlVec: vector of axial induction factor and yaw of turbines
% - rho: air density
% - A: area wind turbine
% - Vinf: free flow wind speed
% - xPos: positions of turbines
% 
% Outputs:
% - Pt: wind farm total power
% - P: vector with wind turbine power
% - V: vecotr with wind speeds at turbines

%% Some optimizer expect vector input
if size(ctrlVec,1) > size(ctrlVec,2) % Some optimizer expect vector input
    ctrlVec = ctrlVec';
end

D =  sqrt(A/pi) *2; % A = pi * (D/2)^2: Diameter from area
nT = (length(ctrlVec)+1)/2; % Prepare cntrl inputs
a_vec = ctrlVec(1:nT); % Axial induction
a_vec = max(min(a_vec,0.33),0);
gamma_vec = [ctrlVec(nT+1:end),0]; % Yaw of turbines

if nargin < 6  || isempty(xPos) % no distance given. Assume line of n turbines with 5D dist.
    xPos = -(0:1:(nT-1))* 5*D; %
end

if nargin < 7 || isempty(b)
    b = 0.075; % model parameters wake interaction
end

%% Overlap 
%gamma = gammaVec(1);
xDist = abs((xPos(2)- xPos(1)));
d = abs( xDist * tan(gamma_vec(1))); % Distance center WT2 and center wake WT1
Dwake = D + 2*b*xDist; % Diameter wake perpendicular to rotated WT1
R = 0.5 * Dwake; %R = R1 Diameter of Wake rotated into the plane of WT2
r = 0.5 * D;

if abs((d^2 + r^2 - R^2)/(2*d*r)) > 1 %<= 0.5*(R-r)
    A12 = r^2*pi;
elseif d >= r + R
    A12 = 0;
else
    A12 = r^2 * acos((d^2 + r^2 - R^2)/(2*d*r)) + ...
        R^2 * acos(( d^2 - r^2 + R^2)/(2*d*R)) - ...
        0.5 *sqrt((4*d^2*R^2)-(R^2-r^2+d^2)^2);     %0.5 *sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
end
A12ToA2 = [1, A12/A];

%% Calculate power of each turbine based on wind and wakes
m = length(a_vec);
V = zeros(m,1);
P = zeros(m,1);
c = zeros(m,1);
deltaV = zeros(m,1);
for idx = 1:m % loop over turbines
    
    %interaction to down stream windturbine
    for idxc = 1:idx-1
        c(idx,idxc) = ( D/ (D + 2* b * (xPos(idxc)-xPos(idx)) ) ) ^2;
    end
    
    deltaV(idx) = sqrt( sum( (c(idx,:).* a_vec * A12ToA2(idx)) .^2 ));
    V(idx) = Vinf * (1 - 2 * deltaV(idx)) * cos(gamma_vec(idx)); %disp(V) ToDo: This is only for two turbines
    Cp     = 4*a_vec(idx)*(1-a_vec(idx)).^2; % power coefficient
    
    P(idx) = 1/2* rho * A* Cp * (V(idx))^3;
end

%%
PT = sum(P);
end
