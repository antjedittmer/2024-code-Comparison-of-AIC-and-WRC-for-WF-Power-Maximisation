%% Set the wind farm parameters
rho = 1.225; % air density
D   = 126; % diameter of each turbine assumed to be same
A   = pi * (D/2)^2; % area swept by the rotor
b   = 0.075; %0.035;% % model parameters wake interaction
Vinf  = 8; % wind unaffected by wind farm

% consider we have formulation like this (wind_dir = 0 deg )
% wind -> T01 --- T02 - T03
% wind Vinf1= V1 ->V2  ->V3
m     = 3; % number of turbines
x     = -(0:1:(m-1))* 5*D; %
c     = zeros(m);

% if a is randomly chosen between 0 and 1/3: a = 1/3*rand(1);
a   = 1/3; % opimal induction rate: baseline greedy approach
amax = 1/2; % maximal induction rate
a_vec0 = a * ones(m,1);

%% Wake interaction coefficient matrix c
% Calculate the fix part of the coefficient mat
% Matrix row i corresponds to turbine i. Matrix column j corresponds to
% wake interaction caused by turbine j. Only sub-diagonal elements are
% non-zero as they correspond to wake effects on downstream turbines.

% Example of C for three turbines
%     0   0   0
% c = c21 0   0
%     c31 c32 0

for idx = 1: m % loop over turbines

    % interaction to down stream windturbine
    for idxc = 1:idx-1
        c(idx,idxc) = ( D/ (D + 2* b * (x(idxc)-x(idx)) ) ) ^2;
    end
end

% c2 = c(1:2,1:2);
c12 = c(2,1);

a2 = 1/3;
Cp2 = 4*a2*(1-a2)^2;
K1 = 3 - 24*c12^2*a2*(1-a2);
K0 = -2 + 12*c12*a2*(1-a2);

% For the roots from the derivative
Ka1_square = 12 - 24*c12^3*Cp2;
Ka1_lin = -16 + 24*c12^2*Cp2;
Ka1_const = 4 - 6*c12*Cp2;
rootP = roots([Ka1_square, Ka1_lin, Ka1_const]);
a1 = rootP(rootP<1)

2*a1*(12-2*c12^3 * Cp2) -16 + 24*c12^2*Cp2
4*(3*a2^2-4*a2 +1)*(-24*c12^3*a1^2 + 24*c12^2*a1-6*c12)
4*(6*a2 -4)*(1 -2*a1*c12)^3

%% 

tmp = load('DataOutAICWRC\pWTopt_AIC.mat');
aWF3 = tmp.aCell{3};
a1 = aWF3(1);
a2 = aWF3(2);
a3 = aWF3(3);

c23 = c12;
c13 = c(3,1);
cp3 = 16/27;
fsqra2 = (a2^2*c23^2 + a1^2*c13^2)^(-0.5);
V_wa2 = (1 - 2*fsqra2^(-1));  %V_wa2 = (1 - 2*(a2^2*c23^2 + a1^2*c13^2)^0.5);

% Check first derivatives:
[df_da1, df_da2, df_da3] = partial_derivatives(a1, a2, a3, c12, c13);
fprintf('%2.4f\n', [df_da1, df_da2, df_da3]);

c_vec = 0:0.01:0.5; 
out = nan(length(c_vec),3);

for idx = 1: length(c_vec)  
    [out(idx,1),~, ~] = partial_derivatives(c_vec(idx),a2, a3, c12, c13);
    [~, out(idx,2),~] = partial_derivatives(a1, c_vec(idx), a3, c12, c13); 
    [~, ~, out(idx,3)] = partial_derivatives(a1, a2, c_vec(idx), c12, c13);
end

figure; 

for idx = 1:3
    subplot(3,1,idx)
    plot(c_vec, out(:,idx), c_vec(1:end-1), diff(out(:,idx)),...
        [c_vec(1) c_vec(end)],[0,0],'k--')
    grid on; axis tight;

    if idx == 1
        hold on; plot(a1,0,'b*'); hold off
    elseif idx == 2
        hold on; plot(a2,0,'b*'); hold off
    else
        hold on; plot(a3,0,'b*'); hold off
    end
end

4*(3*a2^2 - 4*a2 + 1)*(1- 2*a1*c12)^3 - 6*cp3 * ...
a2*c23^2*V_wa2^2 *fsqra2


4*(6*a3 - 4) * (1 - 2*(a2^2*c23^2 + a1^2*c13^2)^0.5)^3


4*(6*a2^2-4)*(1 -2*a1*c12)^3 + 6*cp3 * c23 * V_wa2 * f_a2 * ...
(c23 *a2^2 * V_wa2 * fsqra2^2 - 4*c23*a2^2* fsqra2 - V_wa2)




