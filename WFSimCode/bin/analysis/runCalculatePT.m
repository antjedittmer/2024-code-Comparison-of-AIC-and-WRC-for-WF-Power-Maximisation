clear; close all; clc

fs = 12;

Vinf = 8; % free-stream wind speed
D = 126; % diameter of turbines
nT = 2; % number of turbines
gammaVec = 0*ones(1,nT); % yaw angle
a1Vec = 1/3*ones(1,nT); % induction of turbines
xDistT = 5*D*(0:(nT-1)); 

%% Two turbines: Test whether this gives a similar curve to WFSim
gammaEnd = 30;
gammaVec1 = 0:gammaEnd;
lenP = length(gammaVec1);
PT = nan(lenP,1);
P = nan(lenP,2);
for idx = 0:30
    gammaVec(1)= idx*pi/180;
    [PT(idx+1), P(idx+1,:)] = calculatePT(gammaVec,a1Vec,xDistT,Vinf,D);
end
Psum = sum(P,2);
[Pmax,idxPmax]= max(Psum);
nF = figure;
pos = get(nF,'Position');
set(nF,'Position',[pos(1:3),pos(4)*0.75]);

plot(gammaVec1,P(:,1)/10^6,'bs--',gammaVec1,P(:,2)/10^6,'rs--',...
    gammaVec1,Psum/10^6,'ks--',gammaVec1(idxPmax),Pmax/10^6,'gs')
legend('P WT1','P WT2', 'P T','max(P T)','Location','eastoutside', 'Interpreter','Latex','FontSize',fs);
axis tight;  posAxis = axis; axis([posAxis(1:2),0,posAxis(4)*1.01]);
grid on;

set(gca,'TickLabelInterpreter','Latex','FontSize',fs)
xlabel('yaw $\gamma_1$ [deg]','Interpreter','Latex','FontSize',fs);
ylabel('Power [MW]','Interpreter','Latex','FontSize',fs);


%% Two turbines: Get optimal angles
[gammaOpt,PTopt]  = fminsearch(@(gammaVec) ...
    calculatePT(gammaVec,a1Vec,xDistT,Vinf,D), gammaVec);
PT = abs(PT);
PTopt = abs(PTopt);
disp('Two turbines: Angle and power increase')
disp(gammaOpt * 180/pi)
fprintf('P greedy %2.2f MW, P opt %2.2f MW: %2.2f%% \n\n', min(PT)/10^6, PTopt/10^6, [100*(PTopt - min(PT))/min(PT)]); %, 'perCent'])

%% Three turbines: Reset parameter and get optimal angles for two turbine
nT = 3; % number of turbines
gammaVec = 0*ones(1,nT); % yaw angle

% gammaVec = [26.2510   18.6058    0.0000]
a1Vec = 1/3*ones(1,nT); % induction of turbines
xDistT = 5*D*(0:(nT-1)); 
[gammaOpt,PTopt]  = fminsearch(@(gammaVec) ...
    calculatePT(gammaVec,a1Vec,xDistT,Vinf,D), gammaVec);

PT = abs(calculatePT(gammaVec,a1Vec,xDistT,Vinf,D));

PTopt = abs(PTopt);
disp('Three turbines: Angle and power increase')
disp(gammaOpt * 180/pi)
fprintf('P greedy %2.2f MW, P opt %2.2f MW: %2.2f%% \n', min(PT)/10^6, PTopt/10^6, [100*(PTopt - min(PT))/min(PT)]); %, 'perCent'])


%% Ten turbines: Reset parameter and get optimal angles for two turbine
nT = 4; % number of turbines
gammaVec = 0*ones(1,nT); % yaw angle
a1Vec = 1/3*ones(1,nT); % induction of turbines
xDistT = 5*D*(0:(nT-1)); 
[gammaOpt,PTopt]  = fminsearch(@(gammaVec) ...
    calculatePT(gammaVec,a1Vec,xDistT,Vinf,D), gammaVec);

PT = abs(calculatePT(gammaVec,a1Vec,xDistT,Vinf,D));

PTopt = abs(PTopt);
disp('Ten turbines: Angle and power increase')
disp(gammaOpt * 180/pi)
fprintf('P greedy %2.2f MW, P opt %2.2f MW: %2.2f%% \n', min(PT)/10^6, PTopt/10^6, [100*(PTopt - min(PT))/min(PT)]); %, 'perCent'])

%% Five turbines: Reset parameter and get optimal angles for two turbine
nT = 5; % number of turbines
gammaVec = 0*ones(1,nT); % yaw angle
a1Vec = 1/3*ones(1,nT); % induction of turbines
xDistT = 5*D*(0:(nT-1)); 
[gammaOpt,PTopt]  = fminsearch(@(gammaVec) ...
    calculatePT(gammaVec,a1Vec,xDistT,Vinf,D), gammaVec);

PT = abs(calculatePT(gammaVec,a1Vec,xDistT,Vinf,D));

PTopt = abs(PTopt);
disp('Ten turbines: Angle and power increase')
disp(gammaOpt * 180/pi)
fprintf('P greedy %2.2f MW, P opt %2.2f MW: %2.2f%% \n', min(PT)/10^6, PTopt/10^6, [100*(PTopt - min(PT))/min(PT)]); %, 'perCent'])








