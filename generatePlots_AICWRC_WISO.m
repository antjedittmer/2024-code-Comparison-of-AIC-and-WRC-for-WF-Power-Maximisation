%% Code for 'Data driven Decentralized Algorithm for Wind Farm Control'
% Code illustrating ideas in 'Data driven Decentralized Algorithm for Wind
% Farm Control with Population-Games Assistance'
function generatePlots_AICWRC_WISO(plotAll)

if ~nargin
    clc; clear; close all;
    plotAll = 0; % also generate plots not used in paper  
end

% figDir
figDir ='figDirETC_1';
if ~isfolder(figDir)
    mkdir(figDir)
end

loadAIC = 1;
loadWRC = 1;


mainDir = fileparts(mfilename('fullpath'));
dirCode = fullfile(mainDir,'ParkGaussianFlorisFASTFarm');
addpath(dirCode)
dirFigures = fullfile(mainDir,figDir);
dirProject = fileparts(fileparts(fileparts(dirCode))); % this is my main project folder
dirFloris = fullfile(dirProject,'floris','examples');

dirFASTFarm = fullfile(dirCode,'BinduFASTFarm');
dirFASTFarmFigDir = fullfile(dirFASTFarm, 'Figures');
plotsFASTFarm = {'Farm1Steady85DmTwrBTinOne'};

% Set plotting properties
dAFS = get(0,'DefaultAxesFontSize');
dTFS = get(0,'DefaultTextFontSize');
dLLw = get(0,'DefaultLineLineWidth');

%% Define the parameters
rho = 1.225; % air density
D   = 126; % diameter of each turbine assumed to be same
A   = pi * (D/2)^2; % area swept by the rotor
b   = 0.075; %0.035;% % model parameters wake interaction
gamma = 0.99; %tunable parameter for all strategies
n = 3;
tau = 0.01;
epsilon = eps;
Vinf  = 8; % wind unaffected by wind farm

% consider we have formulation like this (wind_dir = 0 deg )
% wind -> T01 --- T02 - T03
% wind Vinf1= V1 ->V2  ->V3
m     = 100; % number of turbines
x     = -(0:1:(m-1))* 5*D; %
c     = zeros(m);

% if a is randomly chosen between 0 and 1/3: a = 1/3*rand(1);
a   = 1/3 * 0.95; % opimal induction rate: baseline greedy approach
amax = 1/2; % maximal induction rate
a_vec0 = a * ones(m,1)* 0.95;

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

%% Calculate wind farm yield from axial induction vectors sweep
% Sweep all possible combinations between 0 and 0.5 of two wind turbines.

%PTtest2 = nan(combA,1);
mtest2     = 2; % number of turbines
xtest2     = x(1:mtest2); %-(0:1:(mtest2-1))* 5*D; %
ctest2     = c(m-(mtest2-1):m,m-(mtest2-1):m);

mtest3     = 3; % number of turbines
xtest3     = x(1:mtest3); %-(0:1:(mtest3-1))* 5*D; %
ctest3     = c(m-(mtest3-1):m,m-(mtest3-1):m);

mtest4     = 4; % number of turbines
xtest4     = x(1:mtest4); %-(0:1:(mtest3-1))* 5*D; %
ctest4     = c(m-(mtest4-1):m,m-(mtest4-1):m);

% if a is randomly chosen between 0 and 1/2: a = 1/3*rand(1);
a_vecTest2 = 0: 0.01: amax * ones(mtest2,1); % sweep of induction factors a
lenA = length(a_vecTest2); % length induction vector
a_repmat = repmat(a_vecTest2,lenA,1);
a1_vec = a_repmat(:);
a2_vec = repmat(a_vecTest2,1,lenA);

a_vecTest3 = 0: 0.01: amax * ones(3,1); % sweep of induction factors a
lenA = length(a_vecTest2); % length induction vector

a_repmat = repmat(a_vecTest2,lenA^2,1);
a1_vec_4 = a_repmat(:);
a2_repvec = repmat(a_vecTest2,lenA,lenA);
a2_vec_4 = a2_repvec(:);
a3_vec_4 = repmat(a_vecTest2,1,lenA^2);

combA = lenA * lenA; %number of combinations of a1 and a2

PTtest2 = nan(combA,1);
for idx = 1: combA
    [PTtest2(idx)] = calculate_PT([a1_vec(idx),a2_vec(idx)], ctest2, rho, A, Vinf);
end

%% if gamma is chosen between 0 and 30: a = 1/3*rand(1);
gamma_vecTest2 = -5: 0.5: 35* ones(mtest2,1); % sweep of induction factors a
lenGamma = length(gamma_vecTest2); % length induction vector
gamma_repmat = repmat(gamma_vecTest2,lenGamma,1);
gamma1_vec = gamma_repmat(:);
gamma2_vec = repmat(gamma_vecTest2,1,lenGamma);

combGamma = lenGamma * lenGamma; %number of combinations of a1 and a2
PTtestGamma = nan(lenGamma,1);
for idx = 1: combGamma
    % Calculate power of array of three turbines
    [Psum,P] = maxPowerYawSeveralWT1([gamma1_vec(idx),gamma2_vec(idx)]');
    PTtestGamma(idx) = sum(P(1:2)); % Use the first two turbine
end
PTtestGammaGreedy = - maxPowerYawSeveralWT1(0);

%% design the c matrix for three turbines
PTtest3 = nan(combA,1);

for idx = 1: combA
    [PTtest3(idx)] = calculate_PT([a1_vec(idx),a2_vec(idx),1/3], ctest3, rho, A, Vinf);
end

PTalpha = nan(combA,1);

alpha_vec1 = a1_vec*35/max(a1_vec);
alpha_vec2 = a2_vec*35/max(a2_vec);

for idx = 1: combA
    [PTalpha(idx)] = maxPowerYawSeveralWT1([alpha_vec1(idx),alpha_vec2(idx)]');
end

tic
PTtest4 = nan(combA*lenA,1);
for idx = 1: combA*lenA
    [PTtest4(idx)] = calculate_PT([a1_vec_4(idx),a2_vec_4(idx),a3_vec_4(idx), 1/3], ctest4, rho, A, Vinf);
end
toc

%% Visualize wind farm yield from different axial induction sweep
% Plot power yield of possible combinations of axial induction between 0
% and 0.5 of two wind turbines.
cp_additional = abs(maxPowerYawSeveralWT1(0)/calculate_PT([1/3,1/3], ctest2, rho, A, Vinf));

[X,Y] = meshgrid(a_vecTest2,a_vecTest2);
PTmatrix = reshape(PTtest2/10^6,lenA,lenA)  * cp_additional;
idxMax = find(max(PTtest2) == PTtest2);
PTMax_MW = PTtest2(idxMax)/10^6  * cp_additional;

PTmatrix3 = reshape(PTtest3/10^6,lenA,lenA)  * cp_additional;
idxMax3 = find(max(PTtest3) == PTtest3);
PTMax_MW3 = PTtest3(idxMax3)/10^6 * cp_additional;

PTmatrix4 = reshape(PTtest4/10^6,lenA,lenA,lenA) * cp_additional;
idxMax4 = find(max(PTtest4) == PTtest4);
PTMax_MW4 = PTtest4(idxMax4)/10^6 * cp_additional;

[X1,Y1] = meshgrid(a_vecTest2 * 35/max(a_vecTest2),a_vecTest2*35/max(a_vecTest2));
PTalphamatrix = reshape(PTalpha/10^6,lenA,lenA) * cp_additional;
idxMaxA = find(max(PTalphamatrix) == PTalphamatrix);
PTMax_MWalpha = PTalphamatrix(idxMaxA)/10^6 * cp_additional;

[PTGreedy,Pgreedy,vGreedy] = calculate_PT([1/3,1/3], ctest2, rho, A, Vinf);
PT_Greedy_MW = PTGreedy/10^6 * cp_additional;
P_Greedy_MW = Pgreedy/10^6 * cp_additional;
perInc2AIC = P_Greedy_MW * cp_additional;

[PTGreedy3,Pgreedy3,vGreedy3] = calculate_PT([1/3,1/3,1/3], ctest3, rho, A, Vinf);
PT_Greedy_MW3 = PTGreedy3/10^6 * cp_additional;
P_Greedy_MW3 = Pgreedy3/10^6 * cp_additional;

[PTGreedy4,Pgreedy4,vGreedy4] = calculate_PT([1/3,1/3,1/3,1/3], ctest4, rho, A, Vinf);
PT_Greedy_MW4 = PTGreedy4/10^6 * cp_additional;
P_Greedy_MW4 = Pgreedy3/10^6 * cp_additional;

[~,Pmax2,vmax2] = calculate_PT([a1_vec(idxMax),a2_vec(idxMax)], ctest2, rho, A, Vinf);
P_max2_MW = Pmax2/10^6  * cp_additional;

% for gamma
[Xg,Yg] = meshgrid(gamma_vecTest2,gamma_vecTest2);
PTmatrixG = reshape(PTtestGamma/10^6,lenGamma,lenGamma);
[tmp, idxMaxG] = max(PTtestGamma);
PTMaxG_MW = tmp/10^6;

%% Two turbines: Plot farm power vs. axial induction factor
% sGr = sprintf('%2.2f, a = [%2.2f,%2.2f]',[PT_Greedy_MW ,1/3,1/3]);
sGr = sprintf('%2.2f',PT_Greedy_MW);
sMax = sprintf('%2.2f',PTMax_MW);

sGrG = sprintf('%2.2f',PT_Greedy_MW);
sMaxG = sprintf('%2.2f',PTMaxG_MW);

strMaxA = sprintf('[%2.2f, %2.2f]', a1_vec(idxMax),a2_vec(idxMax));
legCell = {'P_{WF} (MW)',...
    sprintf('Greedy P_{WF}: %s MW',sGr),...
    sprintf('Max P_{WF}: %s MW', sMax)};
s1 = sprintf('Greedy P_{WF} = %2.3f at [%2.2f,%2.2f]',[PT_Greedy_MW ,1/3,1/3]);
s2 = sprintf('Max P_{WF} = %2.3f at [%2.2f,%2.2f]',[PTMax_MW, a1_vec(idxMax),a2_vec(idxMax)]);
Phvec = [1:0.1:PTMax_MW,1:0.1:PTMax_MW];

setTitle = 0;
nF = 1;

if plotAll == 1
    figure(nF); nF = nF +1;

    contour(X,Y,PTmatrix,Phvec,'ShowText','on')
    hold on;
    plot(1/3,1/3,'*','MarkerSize',8,...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    plot(a1_vec(idxMax),a2_vec(idxMax),'o','MarkerSize',8,...
        'MarkerEdgeColor','k','MarkerFaceColor',[0,0.9,0]);

    grid on; hold off;
    xlabel('Axial Induction Turbine 1 a_1 (-)')
    ylabel('Axial Induction Turbine 2 a_2 (-)')
    if setTitle == 1
        title({'Two Turbines: Wind Farm Power in MW vs. Induction Factors',[s1,', ',s2]})
        legend('P_{WF} (MW)','Greedy P_{WF} (MW)','Max P_{WF} (MW)','Location','Southeast');
    else
        legend(legCell,'Location','Southeast');
        text(1/3 +0.01,1/3 +0.03,'[0.33,0.33]')
        text(a1_vec(idxMax) +0.01,a2_vec(idxMax) +0.03,strMaxA,'color',[0,0.9,0]);

    end

    figname = 'twoTurbinesGreedyOpt_2d';
    print(fullfile(figDir,figname), '-dpng');
    print(fullfile(figDir,figname), '-depsc');
end

%% Three turbines:
s1_3 = sprintf('Greedy P_{WF} = %2.3f at [%2.2f,%2.2f]',[PT_Greedy_MW3 ,1/3,1/3]);
s2_3 = sprintf('Max P_{WF} = %2.3f at [%2.2f,%2.2f]',[PTMax_MW3, a1_vec(idxMax3),a2_vec(idxMax3)]);

sGr_3 = sprintf('%2.2f',PT_Greedy_MW3);
sMax_3 = sprintf('%2.2f',PTMax_MW3);
strMaxA_3 = sprintf('[%2.2f, %2.2f, 0.33]', a1_vec(idxMax3),a2_vec(idxMax3));
legCell_3 = {'P_{WF} (MW)',...
    sprintf('Greedy P_{WF}: %s MW',sGr_3),...
    sprintf('Max P_{WF}: %s MW', sMax_3)};
Phvec = [1:0.1:PTMax_MW3,1:0.1:PTMax_MW3];

if plotAll == 1
    figure(nF); nF = nF +1;

    contour(X,Y,PTmatrix3,Phvec,'ShowText','on')
    xlabel('Axial Induction Turbine 1 a_1 (-)')
    ylabel('Axial Induction Turbine 2 a_2 (-)')
    grid on; hold on;
    plot(1/3,1/3,'*','MarkerSize',8,...
        'MarkerEdgeColor','k','MarkerFaceColor','k');
    plot(a1_vec(idxMax3),a2_vec(idxMax3),'o','MarkerSize',8,...
        'MarkerEdgeColor','k','MarkerFaceColor',[0,0.9,0]);

    if setTitle == 1
        title({'Three Turbines (a_3 = 1/3): Farm Power [MW] vs. Induction',[s1_3,', ',s2_3]})
        legend('P_{WF} (MW)','Greedy P_{WF} (MW)','Max P_{WF} (MW)','Location','Southeast');
    else
        legend(legCell_3,'Location','Southeast');
        text(1/3 +0.005,1/3 +0.05,'[0.33,0.33,0.33]')
        text(a1_vec(idxMax3) +0.01,a2_vec(idxMax3) +0.03,strMaxA_3,'color',[0,0.9,0]);
    end

    figname = 'twoTurbinesGreedyOpt_3d';
    print(fullfile(figDir,figname), '-dpng');
    print(fullfile(figDir,figname), '-depsc');

end

%% Two and three turbines
pos0 = get(0,'defaultFigurePosition');

yLabelAIC = 'Axial Induction a_2 (-)';
xLabelAIC = 'Axial Induction a_1 (-)';

settingOrig = 1; % this is the original setting

if settingOrig ==1
    posMult = 1.1;
    locLegend = 'Southwest';
    tVec1 = [0.01, -0.03, + 0.01, +0.07];
    tVec2 = [0.005,-0.05, 0.01, -0.01];
else
    posMult = 1.25;
    locLegend = 'EastOutside';
    tVec1 = [-0.04, -0.05, - 0.1, +0.08];
    tVec2 = [0.005,-0.05, 0.01, -0.01];

end

set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',0.7);

figure(nF); nF = nF +1;

tiledlayout(2,1,'TileSpacing','compact')
set(gcf,'position',[pos0(1:3),pos0(4)*posMult])
nexttile
contour(X,Y,PTmatrix,Phvec,'ShowText','on')
hold on;
plot(1/3,1/3,'*','MarkerSize',8,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(a1_vec(idxMax),a2_vec(idxMax),'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[0,0.8,0]);
hold off;
if settingOrig == 0, axis equal; end

ylabel(yLabelAIC)
legend(legCell,'Location',locLegend);
text(1/3 + tVec1(1), 1/3 + tVec1(2),'[0.33,0.33]')
text(a1_vec(idxMax) + tVec1(3),a2_vec(idxMax) + tVec1(4),strMaxA,'color',[0,0.8,0]);

nexttile
contour(X,Y,PTmatrix3,Phvec,'ShowText','on')
xlabel(xLabelAIC), ylabel(yLabelAIC)
grid on;
hold on;
plot(1/3,1/3,'*','MarkerSize',8,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(a1_vec(idxMax3),a2_vec(idxMax3),'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[0,0.8,0]);
hold off;
if settingOrig == 0, axis equal;end

legend(legCell_3,'Location', locLegend);
text(1/3 + tVec2(1), 1/3 + tVec2(2),'[0.33,0.33,0.33]')
text(a1_vec(idxMax3) + tVec2(3),a2_vec(idxMax3) + tVec2(4),strMaxA_3,'color',[0,0.8,0]);

set(0,'DefaultAxesFontSize',dAFS);
set(0,'DefaultTextFontSize',dTFS);
set(0,'DefaultLineLineWidth',dLLw);

figname = 'twoAndThreeTurbinesGreedyOpt_2d';
print(fullfile(figDir,figname), '-dpng');
print(fullfile(figDir,figname), '-depsc');

%% AIC and WRC
settingOrig = 1; % this is the original setting

sGr = sprintf('%2.2f',PT_Greedy_MW);
sMax = sprintf('%2.2f',PTMax_MW);
strMaxA = sprintf('[%2.2f, %2.2f]', a1_vec(idxMax),a2_vec(idxMax));
legCell = {'P_{WF} (MW)',...
    sprintf('Greedy P_{WF}: %s MW',sGr),...
    sprintf('Max P_{WF}: %s MW', sMax)};

Phvec = [0.5:0.1:4,0.5:0.1:4];

if settingOrig ==1
    posMult = 1.1;
    locLegend = 'Southwest';
    tVec1 = [0.01, -0.03, + 0.01, +0.07];
    tVec2 = [0.005,-0.05, 0.01, -0.01];

else
    posMult = 1.25;
    locLegend = 'EastOutside';
    tVec1 = [-0.04, -0.05, - 0.1, +0.08];
    tVec2 = [0.005,-0.05, 0.01, -0.01];

end
xLabelWRC = 'Yaw angle \gamma_1 (deg)';
yLabelWRC = 'Yaw angle \gamma_2 (deg)';

set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',0.7);

%colorbar = gray;

figure(nF); nF = nF +1;

tiledlayout(2,1,'TileSpacing','compact')
set(gcf,'position',[pos0(1:3),pos0(4)*posMult])
nexttile
contour(X,Y,PTmatrix,Phvec,'ShowText','on');
hold on;
plot(1/3,1/3,'*','MarkerSize',8,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(a1_vec(idxMax),a2_vec(idxMax),'o','MarkerSize',8,...
    'MarkerEdgeColor','k','MarkerFaceColor',[0,0.8,0]);
hold off;
if settingOrig == 0, axis equal; end

ylabel(yLabelAIC),xlabel(xLabelAIC)
legend(legCell,'Location',locLegend);
text(1/3 + tVec1(1), 1/3 + tVec1(2),'[0.33,0.33]')
text(a1_vec(idxMax) + tVec1(3),a2_vec(idxMax) + tVec1(4),strMaxA,'color',[0,0.8,0]);

nexttile
contour(Xg,Yg,PTmatrixG,Phvec,'ShowText','on');
xlabel(xLabelWRC), ylabel(yLabelWRC)
grid on;
hold on;
plot(0,0,'*','MarkerSize',8,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
plot(gamma1_vec(idxMaxG),gamma2_vec(idxMaxG),'o','MarkerSize',8,...
    'MarkerEdgeColor','k','MarkerFaceColor',[0,0.8,0]);
hold off;
if settingOrig == 0, axis equal;end

legend(legCell_3,'Location', locLegend);
text(1/3 + tVec2(1), 1/3 + tVec2(2),'[0,0]')
text(a1_vec(idxMax3) + tVec2(3),a2_vec(idxMax3) + tVec2(4),strMaxA_3,'color',[0,0.8,0]);

set(0,'DefaultAxesFontSize',dAFS);
set(0,'DefaultTextFontSize',dTFS);
set(0,'DefaultLineLineWidth',dLLw);

figname = 'twoTurbinesGreedyOptAICWRC_2d';
print(fullfile(figDir,figname), '-dpng');
print(fullfile(figDir,figname), '-depsc');


%% m turbines optimization dP/dt
mEnd = 15;

if exist('pWTopt_AIC.mat','file') == 2 && loadAIC
    load('pWTopt_AIC.mat','PT_opt','PT_greedy','aCell');
else
    dPsum = 0;

    aCell = cell(mEnd,1);
    PT_opt = nan(mEnd,1);
    PT_greedy = nan(mEnd,1);

    for m = 2 : mEnd
        c = zeros(m);
        x = -(0:1:(m-1))* 5*D; %
        for idx = 1 : m % loop over turbines
            % interaction to down stream windturbine
            for idxc = 1:idx-1
                c(idx,idxc) = ( D/ (D + 2* b * (x(idxc)-x(idx)) ) ) ^2;
            end
        end
        a_vec0 = a * ones(m,1);

        a_opt = fmincon(@(a_vec) - calculate_PT(a_vec,c,rho,A,Vinf), a_vec0,[],[],[],[], zeros(size(a_vec0)), amax *ones(size(a_vec0)));
        [PT_opt(m), P_opt, Vopt] = calculate_PT(a_opt, c, rho, A, Vinf);
        [PT_greedy(m), P_greedy] = calculate_PT(1/3*ones(size(a_vec0)), c, rho, A, Vinf);
        %ratioP(m) = (PT_opt(m) - PT_greedy(m))/PT_greedy(m);
        aCell{m} = a_opt;
    end
    save('pWTopt_AIC.mat','PT_opt','PT_greedy','aCell');
end
ratioP = (PT_opt - PT_greedy)./PT_greedy;

if exist('pWTopt_WRC.mat','file') == 2 && loadWRC
    load('pWTopt_WRC.mat','Popt','P0','pWTopt')

else
    nT = 10;
    perMore = NaN(nT,1);
    pWTopt = cell(nT,1);

    Popt = NaN(nT,1);
    P0 = NaN(nT,1);
    x0 = 0;
    for idx = 2:nT
        [pWTopt{idx}, P0a] = patternsearch(@maxPowerYawSeveralWT1,x0);
        x0 =  [pWTopt{idx}(:); pWTopt{idx}(idx-1)];
        Popt(idx) = abs(P0a);
        P0(idx) = abs(maxPowerYawSeveralWT1(0 * pWTopt{idx}));
    end
    save('pWTopt_WRC.mat','Popt','P0','pWTopt')


end

figure(nF); nF = nF +1;

set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.3);

subplot(2,1,1); plot(2:mEnd,PT_greedy(2:mEnd)/10^6,2:mEnd,PT_opt(2:mEnd)/10^6)
title('AIC: WT in line in main wind direction: 1 x 2,3,...100');
grid on;
legend('P_{greedy}','P_{opt}','Location','SouthEast'); ylabel('P_{T} (MW)');
subplot(2,1,2); plot(2:mEnd,100*ratioP(2:mEnd));
grid on; %axis tight
ylabel('Power increase (%)');
xlabel('No. wind turbines (-)')

strName = 'AICImprovement';

print(gcf,fullfile(strName), '-dpng');
print(gcf,fullfile(strName), '-depsc');

% Two turbines: Angle and power increase
%    18.6065    0.0000
%
% P greedy 2.89 MW, P opt 3.10 MW: 7.34%
%
% Three turbines: Angle and power increase
%    26.2510   18.6058    0.0000
%
% P greedy 3.03 MW, P opt 3.49 MW: 15.00%
% Ten turbines: Angle and power increase
%    30.4818   26.2538   18.6094    0.0000
%
% P greedy 3.07 MW, P opt 3.71 MW: 20.71%

% save('pWTopt_WRC.mat','Popt','P0','pWTopt');

perMoreWRC = 100*(Popt - P0)./P0;
ratioP = (PT_opt - PT_greedy)./PT_greedy;
cl = lines;

WRC_floris_unused= [7.39100391 15.053566   18.54092809 20.90854975 23.10266198 24.97313607,...
    26.56170468 28.05012005 29.50473681];

WRC_floris = [ 7.46446676 16.4140829  20.43789747 23.01275746 25.24222557 27.10788875,...
    28.59286803 29.99725078 31.40874593];

mEnd = length(perMoreWRC);
figure; plot(2:mEnd,ratioP(2:mEnd)*100,'o-',2:mEnd, perMoreWRC(2:mEnd),'*-');
hold on; plot(2:mEnd,WRC_floris,'+--', 'color',cl(2,:))

% title('AIC vs WRC: WT in line in main wind direction: 1 x 2,3,4');
grid on;
ylabel('Power increase (%)');
xlabel('No. wind turbines (-)')
legend('AIC, Jensen Model','WRC, Bastankhah Model', 'WRC, FLORIS GCH Optimization','Location','SouthEast')

figname = 'AICWRCImprovement';
print(gcf,fullfile(figDir,figname), '-dpng');
print(gcf,fullfile(figDir,figname), '-depsc');


%% FLORIS plots
% plotAll = 0; % Some plots were created to test how they look like keeping this for now

whatFloris = what(dirFloris);
tmp = load(fullfile(dirFloris,whatFloris.mat{2}));

if plotAll  == 1
    idxX = 1:size(tmp.x1_mesh,1);
    idxY = 1:size(tmp.x1_mesh,2);
    ldyy = tmp.x1_mesh;
    ldxx2 = tmp.x2_mesh;
    sol.u = tmp.vecU_0;
    u_Inf = max(sol.u(:));
    fs = 12;

    % subplot(2,2,[1 3]);
    figure;
    subplot(2,1,1);
    contourf(ldyy(1,idxY),ldxx2(idxX,1),sol.u(idxX,idxY),(0:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(summer); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; colorbar;
    %set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight;

    sol.u = tmp.vecU_opt;
    subplot(2,1,2);
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(0:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; colorbar;
    %set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight;
end

%%
%savemat(mafilename, {'vecU_0': vecU_0, 'vec_V': vecV_0, 'vecU_opt': vecU_opt, 'vecV_opt': vecV_opt,
% 'x1_mesh': x1_mesh, 'x2_mesh': x2_mesh, 'layout_x': layout_x,'layout_y': layout_y1,
% 'farm_power_baseline': farm_power_baseline, 'farm_power_opt': farm_power_opt, 'gamma_turbine': gamma_turbine})

if plotAll  == 1

    idx2 = contains(whatFloris.mat,'2') & contains(whatFloris.mat,'cc');
    mat2 = whatFloris.mat{idx2};
    tmp = load(fullfile(dirFloris,mat2));
    tmp2 = tmp; % save for later;

    idxX = 1:size(tmp.x1_mesh,2);
    idxY = 1:size(tmp.x1_mesh,1);
    ldyy = tmp.x2_mesh';
    ldxx2 = tmp.x1_mesh';
    sol.u = tmp.vecU_0';
    u_Inf = max(sol.u(:));
    u_min = min(sol.u(:));

    figure(nF); nF = nF + 1;
    set(gcf,'position',[pos0(1:3),pos0(4)*1.9])

    set(0,'DefaultAxesFontSize', 12);
    set(0,'DefaultTextFontSize',12);
    set(0,'DefaultLineLineWidth',0.7);

    tiledlayout(2,3,'TileSpacing','compact')

    yStr = 'y (m)';
    nexttile
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; colorbar;
    xlabel(yStr); ylabel('x (m)')
    % set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight;

    sol.u = tmp.vecU_opt';
    nexttile
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; colorbar;
    xlabel(yStr);
    %set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight;

    % Get the mat file containing the 3-turbine array
    idx3 = contains(whatFloris.mat,'3') & contains(whatFloris.mat,'cc');
    mat3 = whatFloris.mat{idx3};
    tmp = load(fullfile(dirFloris,mat3));

    idxX = 1:size(tmp.x1_mesh,2);
    idxY = 1:size(tmp.x1_mesh,1);
    ldyy = tmp.x2_mesh';
    ldxx2 = tmp.x1_mesh';
    sol.u = tmp.vecU_0';
    u_Inf = max(sol.u(:));

    sol.u = tmp.vecU_opt';nexttile
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
    xlabel(yStr);

    % Get the mat file containing the 4-turbine array
    idx4 = contains(whatFloris.mat,'4') & contains(whatFloris.mat,'cc');
    mat4 = whatFloris.mat{idx4};
    tmp = load(fullfile(dirFloris,mat4));

    idxX = 1:size(tmp.x1_mesh,2);
    idxY = 1:size(tmp.x1_mesh,1);
    ldyy = tmp.x2_mesh';
    ldxx2 = tmp.x1_mesh';
    sol.u = tmp.vecU_0';
    u_Inf = max(sol.u(:));

    sol.u = tmp.vecU_opt';nexttile
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
    xlabel(yStr);

    idx5 = contains(whatFloris.mat,'5') & contains(whatFloris.mat,'cc');
    mat5 = whatFloris.mat{idx5};
    tmp = load(fullfile(dirFloris,mat5));

    idxX = 1:size(tmp.x1_mesh,2);
    idxY = 1:size(tmp.x1_mesh,1);
    ldyy = tmp.x2_mesh';
    ldxx2 = tmp.x1_mesh';
    sol.u = tmp.vecU_0';
    u_Inf = max(sol.u(:));

    sol.u = tmp.vecU_opt';nexttile
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
    xlabel(yStr);

    idx6 = contains(whatFloris.mat,'6') & contains(whatFloris.mat,'cc');
    mat6 = whatFloris.mat{idx6};
    tmp = load(fullfile(dirFloris,mat6));

    idxX = 1:size(tmp.x1_mesh,2);
    idxY = 1:size(tmp.x1_mesh,1);
    ldyy = tmp.x2_mesh';
    ldxx2 = tmp.x1_mesh';
    sol.u = tmp.vecU_0';
    u_Inf = max(sol.u(:));

    sol.u = tmp.vecU_opt';nexttile
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
    xlabel(yStr);

    ta = annotation('textarrow'); fs = 12; % plot the air speed arrow
    ta.FontSize = 11;
    ta.Position = [0.20 0.14 0 0.035];
    tmpStr = sprintf('%s%2.1f',' $V_{\infty}$= 8 m/s'); %'$V_{\infty}$'; %
    ta.Text.String = strrep(tmpStr, '$','');
    ta.Text.Interpreter = 'tex';

    figname = 'FlorisWakes2';
    print(gcf,fullfile(figDir,figname), '-dpng');
    print(gcf,fullfile(figDir,figname), '-depsc');
end

%% Floris wake plot



idx2 = contains(whatFloris.mat,'2') & contains(whatFloris.mat,'cc');
mat2 = whatFloris.mat{idx2};
tmp = load(fullfile(dirFloris,mat2));
tmp2 = tmp; % save for later;

idxX = 1:size(tmp.x1_mesh,2);
idxY = 1:size(tmp.x1_mesh,1);
ldyy = tmp.x2_mesh';
ldxx2 = tmp.x1_mesh';
sol.u = tmp.vecU_0';
u_Inf = max(sol.u(:));
u_min = min(sol.u(:));

figure(nF); nF = nF + 1;

set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',0.7);

tiledlayout(1,3,'TileSpacing','compact')
set(gcf,'position',[pos0(1:3),pos0(4)*1.05])

yStr = 'y (m)';
nexttile
contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; colorbar;
xlabel(yStr); ylabel('x (m)')
text(min(ldyy(1,idxY))*0.95, max(ldxx2(idxX,1))*0.95, sprintf('P_{WF} %2.2f MW',tmp2.farm_power_baseline/10^6),'Fontsize',9)

% set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
axis equal; axis tight;

sol.u = tmp.vecU_opt';
nexttile
contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; colorbar;
xlabel(yStr);
text(min(ldyy(1,idxY))*0.95, max(ldxx2(idxX,1))*0.95, sprintf('P_{WF}: %2.2f MW',tmp2.farm_power_opt/10^6),'Fontsize',9)
%set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
axis equal; axis tight;

% Get the mat file containing the 3-turbine array
idx3 = contains(whatFloris.mat,'3') & contains(whatFloris.mat,'cc');
mat3 = whatFloris.mat{idx3};
tmp = load(fullfile(dirFloris,mat3));
idxX = 1:size(tmp.x1_mesh,2);
idxY = 1:size(tmp.x1_mesh,1);
ldyy = tmp.x2_mesh';
ldxx2 = tmp.x1_mesh';
sol.u = tmp.vecU_0';
u_Inf = max(sol.u(:));

sol.u = tmp.vecU_opt';nexttile
contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(u_min:0.1:u_Inf*1.05),'Linecolor','none');
colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
xlabel(yStr);
text(min(ldyy(1,idxY))*0.95, max(ldxx2(idxX,1))*0.95, sprintf('%2.2f MW',tmp.farm_power_opt/10^6),'Fontsize',9)
axis equal; axis tight;

ta = annotation('textarrow'); fs = 12; % plot the air speed arrow
ta.FontSize = 11;
ta.Position = [0.20 0.17 0 0.035];
tmpStr = sprintf('%s%2.1f',' $V_{\infty}$= 8 m/s'); %'$V_{\infty}$'; %
ta.Text.String = strrep(tmpStr, '$','');
ta.Text.Interpreter = 'tex';

figname = 'FlorisWakes';
print(gcf,fullfile(figDir,figname), '-dpng');
print(gcf,fullfile(figDir,figname), '-depsc');


%%

if plotAll  == 1
    idxX = 1:size(tmp.x1_mesh,2);
    idxY = 1:size(tmp.x1_mesh,1);
    ldyy = tmp.x2_mesh';
    ldxx2 = tmp.x1_mesh';
    sol.u = tmp.vecU_0';
    u_Inf = max(sol.u(:));

    figure;
    subplot(1,4,1);
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(0:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
    %set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight;

    sol.u = tmp.vecU_opt';
    subplot(1,4,2);
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(0:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
    %set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight;

    tmp = load(fullfile(dirFloris,whatFloris.mat{4}));
    idxX = 1:size(tmp.x1_mesh,2);
    idxY = 1:size(tmp.x1_mesh,1);
    ldyy = tmp.x2_mesh';
    ldxx2 = tmp.x1_mesh';
    sol.u = tmp.vecU_0';
    u_Inf = max(sol.u(:));

    subplot(1,4,3);
    contourf(ldyy(1,idxY),ldxx2(idxX,1)',sol.u(idxX,idxY),(0:0.1:u_Inf*1.05),'Linecolor','none');
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
    %set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight;
end


%% FASTFarm plots

% Check if all plots exist, otherwise sugest running FASTFarm plot script
% This is not done automatically as it would close all other plots

fileFormats = {'.eps','.png'};

for idx = 1: length(plotsFASTFarm)
    fileFASTFarm = fullfile(dirFASTFarmFigDir,[plotsFASTFarm{idx},'.eps']);

    if exist(fileFASTFarm,'file') == 2

        for idxF = 1: length(fileFormats)
            % Copy png and eps
            aFileFASTFarm = strrep(fileFASTFarm,'.eps',fileFormats{idxF});
            aFilePaper = strrep(aFileFASTFarm, dirFASTFarmFigDir, dirFigures);
            copyfile(fileFASTFarm,aFilePaper,'f')

            % Open fig file. If this does not work catch this
            try
                aFigFile = strrep(fileFASTFarm,'.eps','.fig');
                open(aFigFile); %nF = nF +1;
            catch
            end
        end
    else
        warning( [plotsFASTFarm{idx} 'does not exist. Consider running testVisualizeDifferentYawAngles2.m in folder BinduFASTFarm'])
    end

end

