function hfig = WFSim_animation_compare_h( Wp,solCell,hfig,deltaCell,sigmaYCell,xDistVecCell)
% Import local variables from large structs


if nargin < 5
    plotCenterline = 0;
else
    plotCenterline = 1;
end
%plotCenterline = 0;

Dr     = Wp.turbine.Drotor;
ldyy   = Wp.mesh.ldyy;
ldxx2  = Wp.mesh.ldxx2;
yline  = Wp.mesh.yline;

u_Inf  = Wp.site.u_Inf;

N      = Wp.turbine.N;
Cry    = Wp.turbine.Cry;
Crx    = Wp.turbine.Crx;


sol = solCell{1};

% time   = sol.time;
% k      = sol.k;

%input  = sol.turbInput;
%turb_coord = .5*Dr*exp(1i*input.phi*pi/180);  % Yaw angles
fs = 12;

if N == 2 %plotCenterline &&
    clT = {'b','r'};
    clP = clT;
else
    clT = repmat({'k'},N,1);
    clP = repmat({'r'},N,1);
end

%% Create figure window, if not yet available
if nargin <= 2
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end

set(0,'CurrentFigure',hfig);

pCell = {[1 2],[3,4]};
for idxP = 1: 2
    sol = solCell{idxP};
    delta = deltaCell{idxP};
    sigmaY = sigmaYCell{idxP};
    xDistVec = xDistVecCell{idxP};
       
    idxX = 5:size(sol.u,1)-4;
    idxY = 5:size(sol.u,2)-4;
    
    input  = sol.turbInput;
    turb_coord = .5*Dr*exp(1i*input.phi*pi/180);  % Yaw angles
    
    %% Plot u velocity flow component
    sol_u = sol.u';
    subplot(2,2,pCell{idxP});
    contourf(ldxx2(idxX,1),ldyy(1,idxY ),sol_u(idxY,idxX),(0:0.1:u_Inf*1.05),'Linecolor','none');
    %contourf(ldyy(1,idxY ),ldxx2(idxX,1)',sol.u(idxX,idxY),(0:0.1:u_Inf*1.05),'Linecolor','none');
    
    set(gca,'TickLabelInterpreter','Latex','FontSize',fs);
    colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.05]);  hold on; hc = colorbar;
    set(hc,'TickLabelInterpreter','Latex','FontSize',fs);
    axis equal; axis tight;
    
    % Plot the turbines in the field
    for kk=1:N
        Qy     = (Cry(kk)-real(turb_coord(kk))):1:(Cry(kk)+real(turb_coord(kk)));
        Qx     = linspace(Crx(kk)-imag(turb_coord(kk)),Crx(kk)+imag(turb_coord(kk)),length(Qy));
        plot(Qx,Qy,'k','linewidth',1)
        str = strcat('$WT_{',num2str(kk),'}$');
        text(Crx(kk)-1,Cry(kk)+150,str,'FontSize',fs+2,'interpreter','latex',...
            'Color',clT{kk})
        
        if kk == 1 && plotCenterline == 1
            mIdxQ = round(length(Qy)/2);
            Qy1 = Qy(mIdxQ);
            Qx1 = Qx(mIdxQ);
        end
    end
    
        % text(-600,ldxx2(end,end),['$t=~$ ', num2str(time,'%.1f'), '[s]'],'FontSize',fs,'interpreter','latex');
    
    if idxP == 1 
        ta = annotation('textarrow');
        ta.FontSize = fs;
        ta.Position = [0.24 0.633 0.05 0.0]; %[0.175595238095238 0.229166666666667] [0.649000000000001 0.65]
        ta.Text.String = sprintf('%s%2.1f ',' $V_{\infty}$=',sol.u(1,1)); %'$V_{\infty}$'; %
        ta.Text.Interpreter = 'latex';
    elseif idxP == 2 
        ta1 = annotation('textarrow');
        ta1.FontSize = fs;
        ta1.Position = [0.24 0.153 0.05 0.0]; %[0.175595238095238 0.229166666666667] [0.649000000000001 0.65]
        ta1.Text.String = sprintf('%s%2.1f ',' $V_{\infty}$=',sol.u(1,1)); %'$V_{\infty}$'; %
        ta1.Text.Interpreter = 'latex';
    end
    
    if plotCenterline == 1
        
        if any(isnan(sigmaY)), offset = 0; else
            offset = 36;
        end
        
        % plot(-delta + Qy1, xDistVec' + Qx1, 'b--');
        plot(xDistVec' + Qx1, -delta + Qy1, 'b--');
        %text(Qy1+50,Qx1+100,'$\delta$','FontSize',fs+2,'interpreter','latex',...
        text(Qx1+100,Qy1+50,'$\delta$','FontSize',fs+2,'interpreter','latex',...
            'Color','b')
        
        %plot(-delta(1:end-offset) + Qy1 + sigmaY(1:end-offset), xDistVec(1:end-offset) + Qx1+ offset, 'b-.');
        plot( xDistVec(1:end-offset) + Qx1+ offset, -delta(1:end-offset) + Qy1 + sigmaY(1:end-offset), 'b-.');
        %plot(-delta(1:end-offset) + Qy1 - sigmaY(1:end-offset), xDistVec(1:end-offset) + Qx1+ offset, 'b-.');
        plot( xDistVec(1:end-offset) + Qx1+ offset, -delta(1:end-offset) + Qy1 - sigmaY(1:end-offset), 'b-.');

        sigmaVecX = -delta(1:end-offset) + Qy1 - sigmaY(1:end-offset);
        sigmaVecY = xDistVec(1:end-offset) + Qx1+ offset;
        idx1 = find(~isnan(sigmaVecX),1);
        sigma0X = sigmaVecX(idx1);
        sigma0Y = sigmaVecY(idx1);
        
        % text(sigma0X(1)+100, sigma0Y(1),'$\sigma_y$','FontSize',fs+2,'interpreter','latex',...
        text(sigma0Y(1),sigma0X(1)-50,'$\sigma_y$','FontSize',fs+2,'interpreter','latex',...
            'Color','b')
        
    end
    
    xlabel('$y$ [m]','FontSize',fs+2,'interpreter','latex')
    ylabel('$x$ [m]','FontSize',fs+2,'interpreter','latex');
    set(gca,'TickLabelInterpreter','Latex','FontSize',fs)
    
    
    gamma = sol.turbInput.phi(1);
    %titleStr = {['$k=~$', num2str(round(time),'%d'),', $\gamma_1=~$',num2str(round(gamma),'%d')],'Long. wind $u$ [m/s]'};
    %     titleStr = {'Long. wind $u_x$ [m/s],', ['yaw angle $\gamma_1=~$',num2str(round(gamma),'%d'),' deg']};
    %     title(titleStr,'FontSize',fs,'interpreter','latex');
end