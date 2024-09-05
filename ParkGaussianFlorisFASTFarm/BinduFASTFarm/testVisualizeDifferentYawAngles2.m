%% Script to visualize open loop FASTFarm runs
clc; clear; close all;

% Folder with FAST.Farm outdata
outdataName = 'FAST.Farm1Steady85Dm0to35'; % 5D

% get the file name removing FAST from front and 0to35 from from back
idxDotFast = strfind(outdataName,'.');
idxm = strfind(outdataName,'m');
outmidName = outdataName(idxDotFast+1:idxm(2));

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

%% Load data or read from out files

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

%% %
xlabelStr = 'Yaw $\gamma$ (deg)';
ylabelStrP = '$P$ (kW)';
ylabelStrM = '$M_{TwrY,B}$ (kNm)';



%% 
dirFig = 'Figures';
if ~isdir(dirFig) %#ok<ISDIR> 
    mkdir(dirFig);
end
% Set plotting properties
dAFS = get(0,'DefaultAxesFontSize');
dTFS = get(0,'DefaultTextFontSize');
dLLw = get(0,'DefaultLineLineWidth');

set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',0.7);
 fs = 12; 

figure(1);
t0 = 160;
timeVec = (aTable.Time > t0);

subplot(2,1,1)
plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
%xlabel(xlabelStr,'interpreter','latex'); 
ylabel(ylabelStrP,'interpreter','latex');
legend('WT1','WT2','WF','Location','SouthEast','Orientation','Horizontal','interpreter','latex')
%title(sprintf('Mean values power vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))
set(gca,'TickLabelInterpreter','Latex','FontSize',fs)

subplot(2,1,2)
plot(noYaw,abs(squeeze(mean(towerBMxMatrix(timeVec,:,1:nT)))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrM,'interpreter','latex');xlabel(xlabelStr,'interpreter','latex');
legend('WT1','WT2','Location','East','Orientation','Horizontal','interpreter','latex')
%title(sprintf('Tower base moment caused by side-to-side forces: %2.1f to %2.1f s', t0,aTable.Time(end)))
set(gca,'TickLabelInterpreter','Latex','FontSize',fs)

filenamepng = matlab.lang.makeValidName(outmidName);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');

%% Reset plotting properties to default and plot all tower base moments
set(0,'DefaultAxesFontSize',dAFS);
set(0,'DefaultTextFontSize',dTFS);
set(0,'DefaultLineLineWidth',dLLw);

figure(2);

% save('genPwrMat.mat','genPwrMatrix', 'genPwrFarm', 'timeVec')
cl = lines;

titleStr = ['Mean values power and moments vs. yaw',': WF{\color[rgb]{',num2str(cl(1,:)),'} WT1 ',...
    '\color[rgb]{',num2str(cl(2,:)),'}WT2}'];

t0 = 160;

xlabelStrTex = strrep(xlabelStr,'$','');
ylabelStrPTex = strrep(ylabelStrP,'$','');

nAx = 3;

subplot(nAx+1,1,1)
FarmPwrAtYaw = mean(genPwrFarm(timeVec,:));
powerIncreaseInPerc = 100*(max(FarmPwrAtYaw) - min(FarmPwrAtYaw))/min(FarmPwrAtYaw);
% 2.30 MW, 2.52 MW,  increase 9.7 %

plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrPTex); %xlabel(xlabelStrTex);

title(titleStr)
strCell = {'x','y','z'};

for idx = 1:nAx

    ylabelStrMtex = strrep(strrep(ylabelStrM,'$',''),'Y',upper(strCell{idx}));
    eval(['temp = towerBM',strCell{idx},'Matrix;'])

    subplot(nAx+1, 1,1+idx)
    plot(noYaw,abs(squeeze(mean(temp(timeVec,:,1:nT)))));
    grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.01])
    ylabel(ylabelStrMtex);
end

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.1]);

xlabel(xlabelStrTex);

filenamepng = matlab.lang.makeValidName([outmidName, 'TwrB']);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');

%% 
figure(3)

subplot(nAx+1,1,1)
FarmPwrAtYaw = mean(genPwrFarm(timeVec,:));
powerIncreaseInPerc = 100*(max(FarmPwrAtYaw) - min(FarmPwrAtYaw))/min(FarmPwrAtYaw);

% 2.30 MW, 2.52 MW, increase 9.7 %

plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrPTex); %xlabel(xlabelStrTex);
% legend('WT1','WT2','WF','Location','SouthEast','Orientation','Horizontal','interpreter','latex')
% title(sprintf('Mean values power vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))
title({titleStr,'Tower base: solid lines, top: dashed lines'})

strCell = {'x','y','z'};

for idx = 1:nAx

    ylabelStrMtex = strrep(strrep(ylabelStrM,'$',''),'Y',upper(strCell{idx}));
    eval(['temp = towerBM',strCell{idx},'Matrix;'])
    tempVec = abs(squeeze(mean(temp(timeVec,:,1:nT))));

    eval(['tempT = towerM',strCell{idx},'Matrix;']);
    tempTVec = abs(squeeze(mean(tempT(timeVec,:,1:nT))));

    subplot(nAx+1, 1,1+idx)
    plot(noYaw,tempVec(:,1),'color',cl(1,:))
    hold on;
    plot(noYaw,tempVec(:,2),'color',cl(2,:))
    plot(noYaw,tempTVec(:,1),'--','color',max(cl(1,:) -0.2,0),'linewidth',2);
    plot(noYaw,tempTVec(:,2),'--','color',max(cl(2,:) -0.2,0),'linewidth',2);
    hold off;
    
    ylabelStrMTtex = strrep(ylabelStrMtex,',B','');

    grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.01])
    ylabel(ylabelStrMTtex);
end

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.3]);

xlabel(xlabelStrTex);

filenamepng = matlab.lang.makeValidName([outmidName, 'TwrBTinOne']);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');


%%

figure(4)
subplot(4+3,1,1)
plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrPTex); %xlabel(xlabelStrTex);
% legend('WT1','WT2','WF','Location','SouthEast','Orientation','Horizontal','interpreter','latex')
% title(sprintf('Mean values power vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))
title(titleStr)

strCell = {'x','y','z'};

for idx = 1:3

    ylabelStrMtex = strrep(strrep(ylabelStrM,'$',''),'Y',upper(strCell{idx}));
    eval(['temp = towerBM',strCell{idx},'Matrix;'])

    subplot(4+3,1,1+idx)
    plot(noYaw,abs(squeeze(mean(temp(timeVec,:,1:nT)))));
    grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
    ylabel(ylabelStrMtex);

    eval(['tempT = towerM',strCell{idx},'Matrix;']);
    ylabelStrMTtex = strrep(ylabelStrMtex,',B',',T');

    subplot(4+3,1,4+idx)
    plot(noYaw,abs(squeeze(mean(tempT(timeVec,:,1:nT)))));
    grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
    ylabel(ylabelStrMTtex);

end

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*2.4]);

xlabel(xlabelStrTex);

filenamepng = matlab.lang.makeValidName([outmidName, 'TwrBTwrT']);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');



return;
%% Plots Power vs. yaw

% Set plotting properties
dAFS = get(0,'DefaultAxesFontSize');
dTFS = get(0,'DefaultTextFontSize');
dLLw = get(0,'DefaultLineLineWidth');

set(0,'DefaultAxesFontSize', 14);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultLineLineWidth',1.4);

t0 = 160;
timeVec = (aTable.Time > t0);

hf = figure(1);
xlabelStr = 'yaw (deg)';
ylabelStrP = 'power (W)';

plot(noYaw,squeeze(genPwrMatrix(end,:,1:nT)),noYaw,genPwrFarm(end,:),'k');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrP);xlabel(xlabelStr);
legend('WT1','WT2','WF','Location','SouthEast')
title(sprintf('Last sample power at %2.1f s vs. yaw',aTable.Time(end)))
filenamepng = sprintf('PowerVsYawFASTFarmExample');
print(gcf,filenamepng, '-dpng');

hf = figure(hf.Number + 1);
plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
xlabel(xlabelStr); ylabel(ylabelStrP);
legend('WT1','WT2','WF','Location','SouthEast')
title(sprintf('Mean values power vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))

%% Plot MMoP
hf = figure(hf.Number + 1);
ylabelStrP = 'MMoP (kNm)';

plot(noYaw,abs(squeeze(rootMomentsMatrix(end,:,1:nT))));
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrP);xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('Last sample MOoP at %2.1f s vs. yaw',aTable.Time(end)))

hf = figure(hf.Number + 1);
plot(noYaw,abs(squeeze(mean(rootMomentsMatrix(timeVec,:,1:nT)))));
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
xlabel(xlabelStr); ylabel(ylabelStrP);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('Mean values MOoP vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))

%% Plot TwSS side-to-side tower-top displacement
hf = figure(hf.Number + 2);
ylabelStrP = 'TwSS (m)';

plot(noYaw,squeeze(towerSSMatrix(end,:,1:nT)));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrP);xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('Last tower-top SS deflection at %2.1f s vs. yaw',aTable.Time(end)))

hf = figure(hf.Number + 2);
plot(noYaw,abs(squeeze(mean((towerSSMatrix(timeVec,:,1:nT))))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
xlabel(xlabelStr); ylabel(ylabelStrP);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('tower-top SS deflection vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))

%% Plot TwFA fore-aft tower-top displacement
hf = figure(hf.Number + 2);
ylabelStrP = 'TwFA (m)';

plot(noYaw,squeeze(towerFAMatrix(end,:,1:nT)));
grid on; axis tight;% pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrP);xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('Last fore-aft tower-top displacement at %2.1f s vs. yaw',aTable.Time(end)))

hf = figure(hf.Number + 2);
plot(noYaw,abs(squeeze(mean((towerFAMatrix(timeVec,:,1:nT))))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
xlabel(xlabelStr); ylabel(ylabelStrP);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('tower fore-aft tower-top displacement vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))

%% tower-top / yaw bearing roll moment
hf = figure(hf.Number + 3);
ylabelStrP = 'TwrSS (kNm)';

subplot(3,1,1)
plot(noYaw,abs(squeeze(mean(towerMxMatrix(timeVec,:,1:nT)))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrP);xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('tower-top / yaw bearing roll moment vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))

subplot(3,1,2)
plot(noYaw,abs(squeeze(mean(towerMyMatrix(timeVec,:,1:nT)))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel('TwrFA (kNm)');xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('tower-top / yaw bearing pitch moment vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))

subplot(3,1,3)
plot(noYaw,abs(squeeze(mean(towerMzMatrix(timeVec,:,1:nT)))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel('TwrTorque (kNm)');xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('tower-top / yaw bearing yaw moment vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))



%% moment caused by side-to-side, foreaft and torsional (Tower base)
hf = figure(hf.Number + 3);
ylabelStrP = 'TwrSS (kNm)';

subplot(3,1,1)
plot(noYaw,abs(squeeze(mean(towerBMxMatrix(timeVec,:,1:nT)))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStrP);xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('Tower base moment caused by side-to-side forces: %2.1f to %2.1f s', t0,aTable.Time(end)))

subplot(3,1,2)
plot(noYaw,abs(squeeze(mean(towerBMyMatrix(timeVec,:,1:nT)))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel('TwrFA (kNm)');xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('Tower base moment caused by fore-aft forces: %2.1f to %2.1f s', t0,aTable.Time(end)))
subplot(3,1,3)
plot(noYaw,abs(squeeze(mean(towerBMzMatrix(timeVec,:,1:nT)))));
grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel('TwrTorque (kNm)');xlabel(xlabelStr);
legend('WT1','WT2','Location','SouthEast')
title(sprintf('Tower base yaw (or torsional) moment  vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))

%% Reset plotting properties to default
set(0,'DefaultAxesFontSize', dAFS);
set(0,'DefaultTextFontSize',dTFS);
set(0,'DefaultLineLineWidth',dLLw);

