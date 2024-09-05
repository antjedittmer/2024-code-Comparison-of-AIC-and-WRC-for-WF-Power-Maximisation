function noF = plotFigurePwrMnts1(xlabelStrTex,ylabelStr, titleStr, noYaw, timeVec, genPwrFarm,genPwrMatrix, towerBMxMatrix,towerBMyMatrix, towerBMzMatrix, outmidName,dirFig,noF) %#ok<INUSL> 
% plotFigurePwrMnts1 plots the steady state power of the two turbines and on
% farm level and all tower base moments


if nargin < 9
    noF = 2;
end

figure(noF); noF = noF + 1;

% save('genPwrMat.mat','genPwrMatrix', 'genPwrFarm', 'timeVec')
nT = size(genPwrMatrix,3);
nAx = 3;

subplot(nAx+1,1,1)
FarmPwrAtYaw = mean(genPwrFarm(timeVec,:));
powerIncreaseInPerc = 100*(max(FarmPwrAtYaw) - min(FarmPwrAtYaw))/min(FarmPwrAtYaw); %#ok<NASGU> 
% 2.30 MW, 2.52 MW,  increase 9.7 %

plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStr.PTex); %xlabel(xlabelStrTex);

title(titleStr)
strCell = {'x','y','z'};

for idx = 1:nAx

    ylabelStr.Mtex = strrep(strrep(ylabelStr.M,'$',''),'Y',upper(strCell{idx}));
    eval(['temp = towerBM',strCell{idx},'Matrix;'])

    subplot(nAx+1, 1,1+idx)
    plot(noYaw,abs(squeeze(mean(temp(timeVec,:,1:nT)))));
    grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.01])
    ylabel(ylabelStr.Mtex);
end

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.1]);

xlabel(xlabelStrTex);

filenamepng = matlab.lang.makeValidName([outmidName, 'TwrB']);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');
