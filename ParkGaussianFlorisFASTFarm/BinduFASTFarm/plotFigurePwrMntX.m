function noF = plotFigurePwrMntX(xlabelStr,ylabelStr,noYaw, timeVec, genPwrFarm,genPwrMatrix, towerBMxMatrix, outmidName,dirFig,fs,noF) 
% plotFigurePwrMntX plots the steady state power of the two turbines and on
% farm level and the side-by-side tower base moments

%% Information for plots from data
nT = size(genPwrMatrix,3);

%% Generate figure
figure(noF); noF = noF+ 1;

subplot(2,1,1)
plot(noYaw,squeeze(mean(genPwrMatrix(timeVec,:,1:nT))),'-',noYaw, mean(genPwrFarm(timeVec,:)),'k-');
grid on; axis tight; pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStr.P,'interpreter','latex');
legend('WT1','WT2','WF','Location','SouthEast','Orientation','Horizontal','interpreter','latex')
set(gca,'TickLabelInterpreter','Latex','FontSize',fs)

subplot(2,1,2)
plot(noYaw,abs(squeeze(mean(towerBMxMatrix(timeVec,:,1:nT)))));
grid on; axis tight; % pos1 = axis; axis([pos1(1:2),0,pos1(4)*1.1])
ylabel(ylabelStr.M,'interpreter','latex');
xlabel(xlabelStr,'interpreter','latex');
legend('WT1','WT2','Location','East','Orientation','Horizontal','interpreter','latex')
set(gca,'TickLabelInterpreter','Latex','FontSize',fs)

%% Print figure to png
filenamepng = matlab.lang.makeValidName(outmidName);
saveas(gcf,fullfile(dirFig,filenamepng), 'fig')
print(gcf,fullfile(dirFig,filenamepng), '-dpng');
print(gcf,fullfile(dirFig,filenamepng), '-depsc');

%% Unused code
% % Set title with time information
% title(sprintf('Mean values power vs. yaw : %2.1f to %2.1f s', t0,aTable.Time(end)))
%title(sprintf('Tower base moment caused by side-to-side forces: %2.1f to %2.1f s', t0,aTable.Time(end)))