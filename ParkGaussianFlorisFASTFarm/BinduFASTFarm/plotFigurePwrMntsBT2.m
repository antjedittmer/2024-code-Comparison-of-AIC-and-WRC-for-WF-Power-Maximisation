function noF = plotFigurePwrMntsBT2(xlabelStrTex,ylabelStr, titleStr, noYaw, timeVec, genPwrFarm,genPwrMatrix, towerBMxMatrix,towerBMyMatrix, towerBMzMatrix,...
    towerMxMatrix,towerMyMatrix, towerMzMatrix,outmidName,dirFig,noF) %#ok<INUSL> 
%plotFigurePwrMntsBT2 plots the steady state power of the two turbines and on
% farm level and all tower base and top moments in seven subplots

% save('genPwrMat.mat','genPwrMatrix', 'genPwrFarm', 'timeVec')
nT = size(genPwrMatrix,3);


figure(noF); noF = noF + 1;


subplot(4+3,1,1)
plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStr.PTex); %xlabel(xlabelStrTex);
% legend('WT1','WT2','WF','Location','SouthEast','Orientation','Horizontal','interpreter','latex')
% title(sprintf('Mean values power vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))
title(titleStr)

strCell = {'x','y','z'};

for idx = 1:3

    ylabelStr.Mtex = strrep(strrep(ylabelStr.M,'$',''),'Y',upper(strCell{idx}));
    eval(['temp = towerBM',strCell{idx},'Matrix;'])

    subplot(4+3,1,1+idx)
    plot(noYaw,abs(squeeze(mean(temp(timeVec,:,1:nT))))); %#ok<USENS> 
    grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
    ylabel(ylabelStr.Mtex);

    eval(['tempT = towerM',strCell{idx},'Matrix;']);
    ylabelStr.MTtex = strrep(ylabelStr.Mtex,',B',',T');

    subplot(4+3,1,4+idx)
    plot(noYaw,abs(squeeze(mean(tempT(timeVec,:,1:nT)))));
    grid on; axis tight; %pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
    ylabel(ylabelStr.MTtex);

end

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*2.4]);

xlabel(xlabelStrTex);

filenamepng = matlab.lang.makeValidName([outmidName, 'TwrBTwrT']);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');
