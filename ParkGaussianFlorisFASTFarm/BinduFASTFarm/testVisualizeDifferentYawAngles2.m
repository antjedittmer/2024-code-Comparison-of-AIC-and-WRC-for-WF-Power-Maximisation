%% Script to visualize open loop FASTFarm runs
clc; clear; close all;

% Folder with FAST.Farm outdata
outdataName = 'FAST.Farm1Steady85Dm0to35'; % 5D

% get the file name removing FAST from front and 0to35 from from back
idxDotFast = strfind(outdataName,'.');
idxm = strfind(outdataName,'m');
outmidName = outdataName(idxDotFast+1:idxm(2));

matfileName = 'genPwrMat.mat';


%% Load data or read from out files
if exist(matfileName,'file') == 2
load(matfileName ,'noYaw','timeVec','genPwrFarm','genPwrMatrix','towerBM*Matrix','towerM*Matrix');

else

    %% Set path for Matlab plot files and outputs
plotDir = fileparts(mfilename('fullpath'));
matlabDir = fileparts(plotDir);
parentDir = fileparts(matlabDir);
projectDir = fileparts(fileparts(parentDir));
FASTDir = fullfile(projectDir,'230510_FastFarm','matlab-toolbox');

addpath(genpath(FASTDir));
outdataPath = fullfile(plotDir,outdataName);
disp(['FASTFarm output files in dir: ', outdataPath]);
disp(['Analysis script in: ',plotDir])


% Get information for *.out data (information also used for plotting
oudataDir = dir(outdataPath);
outDirTable = struct2table(oudataDir); % convert to table for easier handling
outdataTable = outDirTable(contains(outDirTable.name,'FAST'),:);
disp(head(outdataTable))

% Find which out files contains turbine and which farm data
outT1 = contains(outdataTable.name, 'T1'); % Turbine 1
outT2 = contains(outdataTable.name, 'T2'); % Turbine 2
outF = ~ (outT1|outT2); % Farm

outFtable = outdataTable(outF,:); % table with files for farm
outT1table = outdataTable(outT1,:);
outT2table = outdataTable(outT2,:);

nY =  sum(outF); % Number yaw angles
noYaw = nan(nY,1); % Initialize vector with yaw angles

%% Read data from file

% Write mat file for farm tables
OutDataFarmAll = struct;
for idx = 1:nY
    outFName = outFtable.name{idx}; % Get file name from table
    idxDot = strfind(outFName,'.'); % Find dots in file name
    gamma1 = str2double(strrep(outFName(idxDot(1)+1: idxDot(2)-1),outmidName,'')); % Use dots to get yaw angle

    fileFarm = fullfile(outdataPath, outFName);  % full file name
    [outDataFarm,outDataChan] = ReadFASTtext(fileFarm); % Read out data
    varname = sprintf('outDataFarmTable%02d', gamma1); % Var name with yaw angle
    OutDataFarmAll.(varname) = array2table(outDataFarm,'VariableNames',outDataChan); % Build table in struct
end
OutDataFarmAll = orderfields(OutDataFarmAll); % order struct fields alphabetically
%save('OutDataFarmAll.mat','OutDataFarmAll')

% Write mat file for turbine table 1
OutDataT = struct('T1All',struct,'T2All',struct);
for idx = 1: nY % Write struct for turbine table 1
    outT1Name = outT1table.name{idx}; % Get file name from table
    idxDot = strfind(outT1Name,'.');
    gamma1 = str2double(strrep(outT1Name(idxDot(1)+1: idxDot(2)-1),outmidName,'')); % Use dots to get yaw angle

    fileT1 = fullfile(outdataPath, outT1Name); % Find dots in file name
    [outDataT1,outDataT1Chan] = ReadFASTtext(fileT1); % Read out data
    varname = sprintf('outDataT1Table%02d', gamma1); % Var name with yaw angle
    OutDataT.T1All.(varname) = array2table(outDataT1,'VariableNames',outDataT1Chan); % Build table in struct
end
OutDataT.T1All = orderfields(OutDataT.T1All); % order struct fields alphabetically

for idx = 1: nY % Write struct for turbine table 2
    outT1Name = outT2table.name{idx}; % Get file name from table
    idxDot = strfind(outT1Name,'.');
    gamma1 = str2double(strrep(outT1Name(idxDot(1)+1: idxDot(2)-1),outmidName,'')); % Use dots to get yaw angle


    fileT1 = fullfile(outdataPath, outT1Name); % Find dots in file name
    [outDataT1,outDataT1Chan] = ReadFASTtext(fileT1); % Read out data
    varname = sprintf('outDataT2Table%02d', gamma1); % Var name with yaw angle
    OutDataT.T2All.(varname) = array2table(outDataT1,'VariableNames',outDataT1Chan); % Build table in struct
end
OutDataT.T2All = orderfields(OutDataT.T2All);

%% Consolidate turbine data into power output and moments table

genPwrTable = table;
rootMomentsTable = table;
towerSSTable = table;
towerFATable = table;
towerMxTable = table;
towerMyTable = table;
towerMzTable = table;
towerBMxTable = table;
towerBMyTable = table;
towerBMzTable = table;
towerYawPznTable = table;
nT = 2;
tmp = regexp(num2str(1:nT),' ','split');
idxWT = tmp(1:2:end);

for iWT = 1:nT % loop over number turbine: 1:2

    iWTa = idxWT{iWT};
    iWTb = ['T',iWTa,'All'];

    structName = OutDataT.(iWTb);
    fTables = fieldnames(structName);

    for idxY = 1:nY % loop over yaw angles

        nameData = fTables{idxY};
        aTable = structName.(nameData);
        yawStr = nameData(end-1:end);
        strYaw = ['Yaw',yawStr];
        noYaw(idxY) = str2double(yawStr);

        varnames1 = aTable.Properties.VariableNames;
        idxPwr = contains(varnames1,'GenPwr'); % Get lines with GenPwr
        idxMOoP = contains(varnames1,'RootMOoP'); % Get lines with OoP moments
        idxYMx = contains(varnames1,'YawBrMxp'); % Get lines with tower-top / yaw bearing roll moment 
        idxYMy = contains(varnames1,'YawBrMyp'); % Get lines with tower-top / yaw bearing pitch moment 
        idxYMz = contains(varnames1,'YawBrMzp'); % Get lines with tower-top / yaw bearing yaw moment 
        idxSS = contains(varnames1,'TTDspSS'); % Get lines with TTower-top side by side deflection
        idxFA = contains(varnames1,'TTDspFA'); % Get lines with Tower FA deflection
        idxBMx = contains(varnames1,'TwrBsMxt'); % Get lines with moment caused by side-to-side forces
        idxBMy = contains(varnames1,'TwrBsMyt'); % Get lines with moment caused by fore-aft forces
        idxBMz = contains(varnames1,'TwrBsMzt'); % Get lines with Tower base yaw (or torsional) moment 
        %idxYaw = contains(varnames1,'YawPzn'); % Get lines with Yaw position 

        genPwrTable.([iWTb,strYaw,'Pwr']) = aTable.(varnames1{idxPwr});
        rootMomentsTable.([iWTb,strYaw,'MOoP']) = mean(aTable{:,varnames1(idxMOoP)},2);
        towerSSTable.([iWTb,strYaw,'SS']) = aTable.(varnames1{idxSS});
        towerFATable.([iWTb,strYaw,'FA']) = aTable.(varnames1{idxFA});
        towerMxTable.([iWTb,strYaw,'SS']) = aTable.(varnames1{idxYMx});
        towerMyTable.([iWTb,strYaw,'SS']) = aTable.(varnames1{idxYMy});
        towerMzTable.([iWTb,strYaw,'SS']) = aTable.(varnames1{idxYMz});
        towerBMxTable.([iWTb,strYaw,'SS']) = aTable.(varnames1{idxBMx});
        towerBMyTable.([iWTb,strYaw,'SS']) = aTable.(varnames1{idxBMy});
        towerBMzTable.([iWTb,strYaw,'SS']) = aTable.(varnames1{idxBMz});
        %towerYawPznTable.([iWTb,strYaw,'SS']) = aTable.(varnames1{idxYaw});
    end
end

%% Further data wrangling: Generate matrices from tables to sum up power
genPwrMatrix = nan( size(genPwrTable,1), size(genPwrTable,2)/nT,nT);
rootMomentsMatrix = nan( size(genPwrTable,1), size(genPwrTable,2)/nT,nT);
towerSSMatrix = nan( size(towerSSTable,1), size(towerSSTable,2)/nT,nT);
towerFAMatrix = nan( size(towerFATable,1), size(towerFATable,2)/nT,nT);
towerMxMatrix = nan( size(towerMxTable,1), size(towerMxTable,2)/nT,nT);
towerMyMatrix = nan( size(towerMyTable,1), size(towerMyTable,2)/nT,nT);
towerMzMatrix = nan( size(towerMzTable,1), size(towerMzTable,2)/nT,nT);
towerBMxMatrix = nan( size(towerBMxTable,1), size(towerBMxTable,2)/nT,nT);
towerBMyMatrix = nan( size(towerBMyTable,1), size(towerBMyTable,2)/nT,nT);
towerBMzMatrix = nan( size(towerBMzTable,1), size(towerBMzTable,2)/nT,nT);

for iWT = 1:nT
    iWTa = idxWT{iWT};
    iWTb = ['T',iWTa];
    idxTurbine = contains(genPwrTable.Properties.VariableNames,iWTb);
    genPwrMatrix(:,:,iWT) = genPwrTable{:,idxTurbine};
    rootMomentsMatrix(:,:,iWT) = rootMomentsTable{:,idxTurbine};
    towerSSMatrix(:,:,iWT) = towerSSTable{:,idxTurbine};
    towerFAMatrix(:,:,iWT) = towerFATable{:,idxTurbine};
    towerMxMatrix(:,:,iWT) = towerMxTable{:,idxTurbine};
    towerMyMatrix(:,:,iWT) = towerMyTable{:,idxTurbine};
    towerMzMatrix(:,:,iWT) = towerMzTable{:,idxTurbine};
    towerBMxMatrix(:,:,iWT) = towerBMxTable{:,idxTurbine};
    towerBMyMatrix(:,:,iWT) = towerBMyTable{:,idxTurbine};
    towerBMzMatrix(:,:,iWT) = towerBMzTable{:,idxTurbine};

end
genPwrFarm = sum(genPwrMatrix,3);
t0 = 160;
timeVec = (aTable.Time > t0);
save(matfileName,'noYaw','timeVec','genPwrFarm','genPwrMatrix','towerBM*Matrix','towerM*Matrix');
end


%% Create figure folder and save default plotting properties
dirFig = 'Figures';
if ~isdir(dirFig) %#ok<ISDIR> 
    mkdir(dirFig);
end
% Set plotting properties
dAFS = get(0,'DefaultAxesFontSize');
dTFS = get(0,'DefaultTextFontSize');
dLLw = get(0,'DefaultLineLineWidth');

% Create labels and title string for plots
xlabelStr = 'Yaw $\gamma$ (deg)';
ylabelStr.P = '$P$ (kW)';
ylabelStr.M = '$M_{TwrY,B}$ (kNm)';

xlabelStrTex = strrep(xlabelStr,'$','');
ylabelStr.MTex  = strrep(ylabelStr.M ,'$','');
ylabelStr.PTex  = strrep(ylabelStr.P ,'$','');

cl = lines;
titleStr = ['Mean values power and moments vs. yaw',': WF{\color[rgb]{',num2str(cl(1,:)),'} WT1 ',...
    '\color[rgb]{',num2str(cl(2,:)),'}WT2}'];

%% Figure 1: Set plotting properties  and plot side-by-side moments
set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',0.7);
fs = 12; 
noF = 1;
noF = plotFigurePwrMntX(xlabelStr,ylabelStr,noYaw, timeVec, genPwrFarm,....
    genPwrMatrix, towerBMxMatrix, outmidName,dirFig,fs,noF);


%% Figure 2: Reset plotting properties to default and plot all tower base moments
set(0,'DefaultAxesFontSize',dAFS);
set(0,'DefaultTextFontSize',dTFS);
set(0,'DefaultLineLineWidth',dLLw);
noF = plotFigurePwrMnts1(xlabelStrTex,ylabelStr, titleStr, noYaw, timeVec,...
    genPwrFarm,genPwrMatrix, towerBMxMatrix,towerBMyMatrix, towerBMzMatrix, ...
    outmidName,dirFig,noF);

%%  Figure 3: Reset plotting properties to default and plot all tower base moments
noF = plotFigurePwrMntsBT(xlabelStrTex,ylabelStr, titleStr, noYaw, timeVec,...
    genPwrFarm,genPwrMatrix, towerBMxMatrix,towerBMyMatrix, towerBMzMatrix,...
    towerMxMatrix,towerMyMatrix, towerMzMatrix,outmidName,dirFig,noF);

%% Figure 4: Plot all tower base moments in seven subplots

