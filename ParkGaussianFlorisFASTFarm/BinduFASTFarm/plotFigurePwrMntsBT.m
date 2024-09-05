function noF = plotFigurePwrMntsBT(xlabelStrTex,ylabelStr, titleStr, noYaw, timeVec, genPwrFarm,genPwrMatrix, towerBMxMatrix,towerBMyMatrix, towerBMzMatrix,...
    towerMxMatrix,towerMyMatrix, towerMzMatrix,outmidName,dirFig,noF) %#ok<INUSL> 
%plotFigurePwrMntsBT plots the steady state power of the two turbines and on
% farm level and all tower base and top moments



% save('genPwrMat.mat','genPwrMatrix', 'genPwrFarm', 'timeVec')
nT = size(genPwrMatrix,3);
nAx = 3;
cl = lines;


figure(noF); noF = noF + 1;

subplot(nAx+1,1,1)
FarmPwrAtYaw = mean(genPwrFarm(timeVec,:));
powerIncreaseInPerc = 100*(max(FarmPwrAtYaw) - min(FarmPwrAtYaw))/min(FarmPwrAtYaw); %#ok<NASGU>  Fod debugging
% 2.30 MW, 2.52 MW, increase 9.7 %

plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStr.PTex); %xlabel(xlabelStrTex);
% legend('WT1','WT2','WF','Location','SouthEast','Orientation','Horizontal','interpreter','latex')
% title(sprintf('Mean values power vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))
title({titleStr,'Tower base: solid lines, top: dashed lines'})

strCell = {'x','y','z'};

for idx = 1:nAx

    ylabelStr.Mtex = strrep(strrep(ylabelStr.M,'$',''),'Y',upper(strCell{idx}));
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
    
    ylabelStr.MTtex = strrep(ylabelStr.Mtex,',B','');

    grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.01])
    ylabel(ylabelStr.MTtex);
end

posDefault = get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.3]);

xlabel(xlabelStrTex);


filenamepng = matlab.lang.makeValidName([outmidName, 'TwrBTinOne']);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');
