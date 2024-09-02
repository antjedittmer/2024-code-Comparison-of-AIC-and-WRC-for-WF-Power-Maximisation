function generatePlots_WFSim(plotAll)

if ~nargin
    clc; clear; close all;
    plotAll = 0; % also generate plots not used in paper  
end

dirFig = 'figureTimeSeriesYawWeight'; % directory for plots with all signals
if ~isfolder(dirFig)
    mkdir(dirFig)
end

paperDir = fullfile(pwd,'figPaperWISO');
if ~isfolder(paperDir)
    mkdir(paperDir)
end


% Load or run open-loop WFSim simulation
dirDataWFSim = 'DataInOutWfSim'; %WFSim run data
dirDataWFSimSubDir1 = 'eDMDresults_UasOutput';
dirDataWFSimDir1 = fullfile(dirDataWFSim,dirDataWFSimSubDir1); %for AIC-WRC
userLoadDataWFSim = 1; % Load Koopman models if available
loadDataWFSim = isdir(dirDataWFSimDir1)&& userLoadDataWFSim; %#ok<ISDIR>

% Load or calculate Koopman sys ID results
dirDataKoopman = 'DataT2OLWFSim'; %Koopman models
dirDataKoopmanSubDir1 = 'Vinf8dot0_CL_K06_P0';
dirDataKoopmanSubDir2 = 'Vinf8dot0_OL_Ct';
dirDataKoopman1 = fullfile(dirDataKoopman,dirDataKoopmanSubDir1); %for AIC-WRC
dirDataKoopman2 = fullfile(dirDataKoopman,dirDataKoopmanSubDir2); %for AIC-WRC
userLoadDataKoopman = 1; % Load Koopman models if available
loadDataKoopman = isdir(dirDataKoopman1) && isdir(dirDataKoopman2) && userLoadDataKoopman; %#ok<ISDIR>

% Load or run closed-loop WFSim simulation
WFSimData = fullfile(dirDataWFSim, 'WFSimSimulationData.mat');
userLoadData = 1; % load WFSim data if available
loadData = exist(WFSimData,'file') == 2 && userLoadData;

%% 1. Create or load open loop test data with WFSim

% Inputs
R = 1e-6; % Weights on control input changed (J = sum(e'Qe + dU'Rdu)
refstairs = 0; %refstairs: Set reference to stairs
measured = 1; %measured (for Koopman model): Use measured values as feedback
KoopmanStates = 6; %KoopmanStates (for Koopman model):  Number of Koopman states
PolyLiftingFunction = 0; %PolyLiftingFunction(for Koopman model): Higher order lifting function 2,4,6,10,12,14,18,24
controller = 0; %controller: integer switch: Open-loop: 0,  Wfsim NMPC: 1, KMPC: 2
Vinf = 8;
vinfstr = '';

if ~loadDataWFSim
    % Run WFSim_demo in open loop for identification set
    ControlSetStr = 'sowfa_2turb_yaw_noise_step';
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf,vinfstr);
    fig = gcf; fig.Name = sprintf('OL_%s',ControlSetStr);

    % RunWFsim demo for validation set
    ControlSetStr = 'sowfa_2turb_yaw_steps_Ct_comb';
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf,vinfstr);
    fig = gcf; fig.Name = sprintf('OL_%s',ControlSetStr);

    %Run WFSimdemo for open loop yaw step
    ControlSetStr = 'sowfa_2turb_yaw_steps';
    WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,controller,ControlSetStr,Vinf,vinfstr);
    fig = gcf; fig.Name = sprintf('OL_%s',ControlSetStr);

    runWFSimAnimationPlots;
end
%% 2. Generate Koopman Sys Id for windfarm controller

if ~loadDataKoopman
    % Inputs
    yawmode = 1; %0 for ct (pitch/torque) control, 1 for additional yaw control
    filenameId = 'Vinf8dot0_sowfa_2turb_yaw_noise_step';
    filenameVal = 'Vinf8dot0_sowfa_2turb_yaw_steps'; %'Vinf8dot0_sowfa_2turb_alm_turbl.mat'; %
    noStates = 6; %number of Koopman states
    useVal = 1; % use extra data for validation
    percentTrain = .6;
    NTMstr = ''; % normal turbulent wind (if empty constant free flow wind)NTMstr
    PolyVec = 0;
    t0 = 240;

    % Run Koopman main function for WFSim simulation environent
    %MainWfSim(yawmode,filenameId,filenameVal,noStates,useVal,percentTrain);
    MainWfSimVinYawPowerOut(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0);
    MainWfSimVinYaw(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)% Code based on
    yawmode = 0; %0 for ct (pitch/torque) control, 1 for additional yaw control

    filenameId = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb';
    filenameVal = 'Vinf8dot0_sowfa_2turb_alm_turbl_AllComb';
    MainWfSimVin(yawmode,filenameId,filenameVal,noStates,PolyVec,useVal,percentTrain,NTMstr,t0)% Code based on
end

%% 3. Evaluate quality of AIC Koopman Sys ID in WFsim in closed loop with MPC

% Inputs
R = 1e-6; % Weights on control input changed (J = sum(e'Qe + dU'Rdu)
refstairs = -1; %refstairs: Set reference to stairs. If -1
measured = 0; %measured (for Koopman model): Use measured values as feedback
controller = 2; %controller: integer switch: Open-loop: 0,  Wfsim NMPC: 1, KMPC: 2
Vinf = 8;
ControlSetStr = 'sowfa_2turb_alm_turbl';
PolyLiftingFunction = 0; %PolyLiftingFunction(for Koopman model): Higher order lifting function 2,4,6,10,12,14,18,24Force
KoopmanStates = 6; %KoopmanStates (for Koopman model):  Number of Koopman states
vinfStr = '';% '_Vinf'; %_Vinf'; % alternativ ''; % Use
kTimestep =  2;
sol_update = 0;

% Load data or run WFSim_demo with AIC
OptYaw = 0; % Keep the yaw angle at zero
KoopAIC = 1; % Use Koopman matrix for AIC
idx = 1;
time0 = 240;

if loadData
    load(WFSimData,'CT_prime*', 'Phi*', 'Force*','Power*', 'time0' , 'JR2', 'JQ2','mpc_Pref');
else
    % run simulation
    [sol_array,JR2AIC(idx),JQ2AIC(idx),fileName,mpc]  = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
        controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC);

    % assign the simulated signals from the AIC run
    CT_prime = cell(1); Phi = cell(1); Force = cell(1); Power = cell(1);
    CT_prime{idx} = cell2mat(arrayfun(@(x)(x.turbine.CT_prime),sol_array,'UniformOutput', false));
    Phi{idx} = cell2mat(arrayfun(@(x)(x.turbine.Phi),sol_array,'UniformOutput', false));
    Force{idx} = cell2mat(arrayfun(@(x)(x.turbine.force),sol_array,'UniformOutput', false));
    Power{idx} = cell2mat(arrayfun(@(x)(x.turbine.power),sol_array,'UniformOutput', false));
end

% calculate the mean values
PhiVecMeanAIC(idx) = mean(mean(abs(Phi{idx}(1,time0:end)')));
ForceMeanAIC(idx) = mean(mean(abs(Force{idx}(1,time0:end)')));

% print to figure folder
fig = gcf;
fig.Name = sprintf('TimeSeries_OptYaw%d_KoopAIC%d',OptYaw,KoopAIC);
filenamepng = fig.Name;
print(gcf,fullfile(dirFig,filenamepng), '-dpng');

%% 4. Run WFSim_demo with AIC + WRC with power estimates
OptYaw = 1; % Yaw angle as control input
KoopAIC = 0; % Use Koopman matrix for AIC
R = 10^-6; % R weight
Rvec = [1,1,0.1]; % Weights on control inputs Ct1, Ct2, and gamma1
RUvec = [0, 10, 100, 1000];
lRU = length(RUvec);


% initialize values for for-loop
if ~loadData
    CT_prime_RU = cell(lRU,1); %cells forsimulated signals from AIC-WRC run
    Phi_RU = cell(lRU,1);
    Force_RU = cell(lRU,1);
    Power_RU = cell(lRU,1);


    PhiVecMean = nan(lRU,1);
    ForceMean = nan(lRU,1);
    JR2 = nan(lRU,1);
    JQ2 = nan(lRU,1);
end

for idx = 1 :length(RUvec)
    RU = RUvec(idx);

    if ~loadData
        [sol_array,JR2(idx),JQ2(idx),fileName,mpc]  = WFSim_demo(R,refstairs,measured,KoopmanStates,PolyLiftingFunction,...
            controller,ControlSetStr,Vinf, vinfStr,kTimestep,sol_update,OptYaw,KoopAIC,Rvec,RU);
        mpc_Pref = mpc.Pref;

        CT_prime_RU{idx} = cell2mat(arrayfun(@(x)(x.turbine.CT_prime),sol_array,'UniformOutput', false));
        Phi_RU{idx} = cell2mat(arrayfun(@(x)(x.turbine.Phi),sol_array,'UniformOutput', false));
        Force_RU{idx} = cell2mat(arrayfun(@(x)(x.turbine.force),sol_array,'UniformOutput', false));
        Power_RU{idx} = cell2mat(arrayfun(@(x)(x.turbine.power),sol_array,'UniformOutput', false));
    end

    PhiVecMean(idx) = mean(mean(abs(Phi_RU{idx}(1,time0:end)')));
    ForceMean(idx) = mean(mean(abs(Force_RU{idx}(1,time0:end)')));

    % print to figure folder
    fig = gcf;
    fig.Name = sprintf('TimeSeries_RU%02d',RU);
    filenamepng =  fig.Name;
    print(gcf,fullfile(dirFig,filenamepng), '-dpng');
    print(gcf,fullfile(dirFig,filenamepng), '-depsc');

end

%% 5.a) Settings for plots
dAFS = get(0,'DefaultAxesFontSize');
dTFS = get(0,'DefaultTextFontSize');
dLLw = get(0,'DefaultLineLineWidth');

set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1.3);

lc = lines;
ls = {'-.','--',':'};
TCCell = {'TC3 $r_{\gamma0}=10^{-2}$','TC2 $r_{\gamma0}=10^{-3}$','TC1 $r_{\gamma0}$=0'};
TCCell = fliplr(TCCell);
CriteriaCell = {'Yaw $\gamma$ [$^{\circ}$]', 'Actuator activity $\Delta U$ $\cdot$ 1000 [-]', ...
    'Power tracking error $\bar{e}$ $\cdot$ 100 [kW]'};
CriteriaCellShort = {'$\gamma$ [$^{\circ}$]', ' $\Delta U$ $\cdot$ 1000 [-]', ...
    '$e$ $\cdot$ 100 [kW]'};

XCriteriaCell = categorical(CriteriaCell);
idxTC = [1,3,4];

%% 5.b. Comparison time series
titleON = 0;
fs = 12;
default_position = get(0, 'DefaultFigurePosition');

timeIdx = 500: length(Power_RU{idxTC(1)});
titleStr = 'Different $R_{\gamma}$ for Cost: $J = e^T Q e + \Delta U^T R \Delta U + U ^T R_\gamma U $';
yStrP = '$P_{WF} (MW)$';
yStrGamma = '$\gamma_1$ (deg)'; % (${\circ})$]';
legStr = ['$P_{ref}$','BL AIC', TCCell];
interStr = 'tex';
fs2 = fs;
if strcmp(interStr ,'tex')
    titleStr = strrep(titleStr,'$','');
    yStrP = strrep(yStrP,'$','');
    yStrGamma = strrep(yStrGamma,'$','');
    legStr = strrep(legStr,'$','');
end

figure(100)
set(gcf, 'Position', [default_position(1:3), 1.35*default_position(4)]);
tiledlayout(5,1,'TileSpacing','compact')
nexttile([2 1])

plot(timeIdx,mpc_Pref(timeIdx)/10^6,'k'); hold on;
tmpPAIC = sum(Power{1});
plot(timeIdx, tmpPAIC(timeIdx)/10^6,'color',[0.5, 0.5,0.5],'LineWidth',1.5,'LineStyle','--'); hold on;

eRMS = nan(length(idxTC)+1,1);
errorPVec = (mpc_Pref(timeIdx)/10^6 - tmpPAIC(timeIdx)'/10^6);
eRMS(1) = sqrt(errorPVec' * errorPVec / length(errorPVec));

for idx = 1: length(idxTC)
    tmpP = sum(Power_RU{idxTC(idx)});
    plot(timeIdx, tmpP(timeIdx)/10^6,'color',lc(idx,:),'LineStyle',ls{idx}); hold on;
    axis tight; grid on; ylabel(yStrP,'Interpreter', interStr,'FontSize',fs2)

    errorPVec = (mpc_Pref(timeIdx)/10^6 - tmpP(timeIdx)'/10^6);
    eRMS(idx+1) = sqrt(errorPVec' * errorPVec/ length(errorPVec));
end

if titleON == 1
    title(titleStr,'Interpreter', interStr,'FontSize',fs2)
end
legend_handle = legend(legStr,'interpreter',interStr,'FontSize',fs-1.5,...
    'Location','South','NumColumns',2);
legend_position = legend_handle.Position; % get the current position of the legend
legend_position(1) = legend_position(1) - 0.05; % move left by 0.05 normalized units
legend_handle.Position = legend_position; % re-assign

% second subfigure
nexttile([1 1])
phiRMS = zeros(length(idxTC)+1,1);
diffPhiRMS = zeros(length(idxTC)+1,1);

for idx = 1: length(idxTC)
    tmpPhi = Phi_RU{idxTC(idx)};
    tmpPhi1 = tmpPhi(1,timeIdx);

    plot(timeIdx, tmpPhi1,'color',lc(idx,:),'LineStyle',ls{idx}); hold on;
    axis tight; grid on; ylabel(yStrGamma,'Interpreter',interStr,'FontSize',fs)

    phiRMS(idx+1) = sqrt(tmpPhi1 * tmpPhi1'/ length(tmpPhi1));
    diffPhi = diff(tmpPhi1);
    diffPhiRMS(idx+1) = sqrt( diffPhi *  diffPhi'/ length( diffPhi));

end

% Plot the CT curves
diffCTRMS = zeros(length(idxTC)+1,2); %

for idx0 = 1:2
    str0 = [sprintf('$$C_{T%d}', idx0),'^{\prime}$$ (-)']; % $C_{T1}^{\prime}$
    if strcmp(interStr ,'tex'), str0 = strrep(str0,'$',''); end

    nexttile([1 1])
    tmpCTAIC = CT_prime{1}(idx0,:);
    plot(timeIdx, tmpCTAIC(timeIdx),'color',[0.5, 0.5,0.5],'LineWidth',1.5,'LineStyle','--'); hold on;

    diffCt1 = diff(tmpCTAIC);
    diffCTRMS(1,idx0) = sqrt( diffCt1 * diffCt1'/ length( diffCt1));

    for idx = 1: length(idxTC)
        tmpCT = CT_prime_RU{idxTC(idx)};
        plot(timeIdx, tmpCT(idx0,timeIdx),'color',lc(idx,:),'LineStyle',ls{idx}); hold on;

        diffCt1 = diff(tmpCT); % Difference Ct
        diffCTRMS(idx+1,idx0) = sqrt( diffCt1 * diffCt1'/ length( diffCt1)); %RMS

        axis tight; grid on; ylabel(str0,'Interpreter', interStr,'FontSize',fs2)
    end
    if idx0 == 2

        xlabel('time (s)','Interpreter', interStr,'FontSize',fs)
    end
end

if titleON == 1, filenamepng = sprintf('timeseriesDifferentRgammaScaledCt');
else, filenamepng = sprintf('timeseriesDifferentRgammaScaledCt0');end

print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');

%% Table values
diffU = sum(diffCTRMS,2) + 0.1 * diffPhiRMS; % weighted sum delta U
eRMSKW = eRMS * 10^3; % error in MW
disp(['eRMSKW &', sprintf('%2.2f &', eRMSKW)])
disp(['phiRMS &', sprintf('%2.2f &', phiRMS)])
disp(['diffU &', sprintf('%2.2f &', diffU)])

%% bar plot

if plotAll == 1
    figure(101);

    X = categorical(CriteriaCellShort);
    X = reordercats(X,CriteriaCellShort);
    if strcmp(interStr,'tex')
        TCCell = strrep(TCCell,'$','');
    end

    b = bar(X,[PhiVecMean(idxTC), JR2(idxTC)*1000, JQ2(idxTC) *100]);
    set(gca,'TickLabelInterpreter',interStr,'FontSize',fs)


    legend(TCCell,'interpreter',interStr,'FontSize',fs) % 'Location','South'
    grid on;

    filenamepng = sprintf('barDifferentRgammaScaled');
    print(gcf,fullfile(dirFig,filenamepng), '-dpng');
    print(gcf,fullfile(dirFig,filenamepng), '-depsc');
end

%% Copy png and eps in one folder

% Manually set file names for paper plots
origDir = fullfile(pwd,dirFig); %DataT2OLWFSim','Vinf8dot0_CL_K06_P0');

% Manually set file names for paper plots
filenames = {'timeseriesDifferentRgammaScaledCt0'};
filetype  = {'.png','.eps'};

for idx = 1: length(filenames)
    for idxT = 1: length(filetype)
        SOURCE = fullfile(origDir,[filenames{idx},filetype{idxT}]);
        copyfile(SOURCE, paperDir)
    end
end



%% Unused code

% % Change color of bar plot
% cl1 = lines;
% cl2 = [0.5,0.5,0.5; cl1(4,:); [0.8,0,0]; cl1(1,:); cl1(5,:)];
% cl =[0,0,0.9; 0,0.9,0; 0.9,0,0];
% for k = 1:3, b(k).FaceColor = cl(k,:); end
%
% title({'Normalized average yaw, force and power trackin error', ...
%     'Cost function: $J = e^T Q e + \Delta U ^T R \Delta U + U^T R_U U$'},...
%     'Interpreter', 'Latex')

% % Possible different labels and legends for plots
% legend('Yaw $\gamma$ [$^{\circ}$]', 'Actuator activity $\Delta U$ $\cdot$ 1000 [-]', ...
%    'Power tracking error $\bar{e}$ $\cdot$ 100 [kW]', 'Location','NorthWest','Interpreter', interStr,'FontSize',fs)
%
% ylabel('Scaled Minimization Criteria','Interpreter', interStr,'FontSize',fs)
% set(gca,'TickLabelInterpreter',interStr,'FontSize',fs)
%
% filenamepng = sprintf('barDifferentRgammaScaled');
% print(gcf,fullfile(dirFig,filenamepng), '-dpng');