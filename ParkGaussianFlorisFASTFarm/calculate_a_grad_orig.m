function [a_grad,p,a_best] = calculate_a_grad_orig(a_k,del,p,PT_I,plotOff,c,rho,A,Vinf)
% calculate_a_grad_2 calculates the new axial induction vector
%
% Inputs
% - a_k: axial induction factor vectors current time step
% - del: power output vector current time step
% - theta1: multiplication coefficient update
% - del: bias coefficient update
% Outputs
% - a_grad: axial induction factor vector next time step

if nargin <5
    plotOff = 1;
end

%% Calculate change rate first turbine
m = size(a_k,1);
sum_g1 = zeros(m,1);
l_N = 3;  % neighborhood turbine 1: 2-3

% PT_I = nan(l_N,1);
% PT_I(1) = calculate_PT(a_k(:,1), c, rho, A, Vinf);


for idxn = 2:l_N %idxn stands for index_neighbors
    %PT_I(idxn) = calculate_PT(a_k(:,idxn), c, rho, A, Vinf);
    multPT = sign(PT_I(idxn) - PT_I(1)); % PT_l - PT_L /abs( PT_l - PT_L)
    g1     = multPT * (a_k(:,idxn)-a_k(:,1));
    sum_g1 = sum_g1 + g1;
end

PI_ub = max(PT_I)*(1 + eps); %2* 10^6 * m; %max(PT_I)+1; 

f_PT_I = (PT_I - PI_ub); %f_pvec = (PT_I - PI_ub) .* p;
DforL = diag([2,1,1]);
AforL = [0,1,1;1,0,0; 1,0,0];
L = DforL - AforL;

%% Calculate trajectory pdot
% ydot = L *[f1*p1, f2*p2, f3*p3]';
% ydot = [l11 * f1, l12 * f2, l13 * f3;
%         l21 * f1, l22 * f2, l23 * f3;
%         l31 * f1, l32 * f2, l33 * f3] * [p1,p2,p3]'
% xdot(t) = A*x(t); x(0) = x0;
% x(t) = C exp(A * t) = x0 *exp(A*t);

L_f_PTI = L .* repmat(f_PT_I',3,1);

tspan = [0 0.01];

y0 = p;

if plotOff
    [~,y] = ode45(@(t,y) L_f_PTI *y, tspan, y0);
else
    [t,y] = ode45(@(t,y) L_f_PTI *y, tspan, y0);
    figure(1)
    plot(t,y(:,1),t,y(:,2),'--',t,y(:,3),'-.');
    legend('p1','p2','p3');
    title('Probabilities strategies with p''(t) = L f(p(t)) ')
    xlabel('time (s)'); ylabel('Probability [-]');
    axis tight; grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1.5)
end


% y = y0 * exp(L_f_PTI*0)
% Theta_k = 1 - p(tau) -> p* = p(tau) = 1;
p = y(end,:);
Theta = (1 - p); %prob. 1,2,3

phi = Theta(1) *sum_g1 + del; %prob 1: Bad (p(1 -> Big change 

%% Update axial induction all turbines 
a_grad = nan(size(a_k));
a_grad(:,1) = a_k(:,1) +  phi;
for idx = 2: l_N
        a_grad(:,idx) = a_k(:,idx-1);
end

%a_best = a_grad(:,max(PT_I) == PT_I);
% calculate_PT(a_vec, c, rho, A, Vinf)
PTNew = calculate_PT(a_grad(:,1), c, rho, A, Vinf);
if PTNew > max(PT_I)
    a_best = a_grad(:,1);
else
    i1 = find(PT_I == max(PT_I));
    a_best =  a_k(:,i1(1));
    a_grad(:,2) =   a_best; % a_k(:,i1(1));

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