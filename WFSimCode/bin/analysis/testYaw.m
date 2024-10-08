close all; clear; clc;

angleVec = 0:5:35;

dAFS = get(0,'DefaultAxesFontSize');
dTFS = get(0,'DefaultTextFontSize');
dLLw = get(0,'DefaultLineLineWidth');

set(0,'DefaultAxesFontSize', 12);
set(0,'DefaultTextFontSize',12);
set(0,'DefaultLineLineWidth',1);
fs = 12;

%% Define parameters
%a1Vec = repmat(1/3,1, length(gammaVec))*0.9; %
% ctVec = [0,15,30,45;0.79, 0.76,0.69,0.58]; %temp = interp1([0,15,30,45],[0.79, 0.76,0.69,0.58],0:45);
% cpVec = [0,15,30,45; 0.5, 0.48,0.42,0.26]; %temp = interp1([0,15,30,45],[0.5, 0.48,0.4,0.26],0:45);

gammaVec = (0:5:35); % change
alphaStar =  2.32; %32; % 4*alpha = 4* 0.58 = 2.32 2.8
betaStar = 0.154; % 2*beta = 0.154 from jetstreams
b   = 0.022; % 0.075; % model parameters wake interaction 0.022
I0 = 0.05; %Wp.site.v_Inf/Wp.site.u_Inf; % Intensity: 0 for our case

rho = 1.225; %Wp.site.Rho; % air density
D   = 126; %Wp.turbine.Drotor;  % rotor diameter
r   = D/2; %rotor radius
A   = r^2 *pi;% rotor area

Vinf  =  8; %Wp.site.u_Inf; % wind unaffected by wind farm
xDist = 5*D; % diff(Wp.turbine.Crx); %  Distance between turbines
%c12 = ( D/ (D + 2* b * xDist ) )^2;

a2 = 1/3*0.95; % Induction factor WT2 for maximal power  1/9; %
Ct2 = 4*a2 *(1- a2); % eq. 6.2 thrust coefficient for max. power
a1_0 = 1/3*0.95; %
gammaCorVec = 0:35;
aCorVec = linspace(1/3, 1/3.5, length(gammaCorVec)) / (1/3);


A12Vec = nan(length(gammaVec),1);
A12ToA1 = nan(length(gammaVec),1);
Ct1Vec = nan(length(gammaVec),1);
P = nan(length(gammaVec),2);

% For centerline and wake expansion
tempXDist = xDist; %this is the currently used distance
nD = ceil(10*D);
xDistVec = 0 : nD;
deltaVec = nan(nD,length(gammaVec));
sigmaYVec = nan(nD,length(gammaVec));
sigmaZVec = D/sqrt(8)*ones(nD,length(gammaVec));
x0Vec = nan(length(gammaVec),1);

%figure; plot(diff(vtempT))
Ivec = nan(length(gammaVec),1);
Vvec = nan(length(gammaVec),2);
ThetaCoVec = nan(length(gammaVec),1);

for idx = 1 : length(gammaVec)
    gammaDeg = gammaVec(idx); % current yaw angle WT1
    a1 = a1_0 * interp1(gammaCorVec,aCorVec,gammaDeg);
    gamma =  gammaDeg/180*pi;
    
    % Calculate variables based on paper
    % eq.6.7 Thrust coefficient (approximation of eq. 6.1)
    % Ct1 = 4*a1 *(1- a1*cos(gamma)); %* interp1(ctVec(1,:),ctVec(2,:),gammaDeg);
    Ct1 = 4*a1 * sqrt(1- a1*(2*cos(gamma)-a1)); %* eq. 6.1 interp1(ctVec(1,:),ctVec(2,:),gammaDeg)/ctVec(2,1);
    % eq. 6.12: Wake angle at rotor 30 -> 6 deg
    ThetaCo = - 0.3*gamma/cos(gamma) *(1- sqrt(1-Ct1*cos(gamma)));
    
    % Question phi = - sin(gamma)*Ct1/(2*(1+ b * xDist/D)^2): 14 deg >> 6 deg
    I = I0 +   0.0150 * gammaDeg/5; %(idx-1); %%abs(sol_array(kPlotVec(idx)).v(nTx,nTy)/sol_array(kPlotVec(idx)).u(nTx,nTy));
    Ivec(idx) = I; 
    % end of near wake area x0
    sqrt1minCt = sqrt(1-Ct1);
    x0 = cos(gamma)* (1 + sqrt1minCt) / ...
        (sqrt(2)*(alphaStar*I + betaStar*(1 - sqrt1minCt))) * D; % eq. 6.16/ 7.3
    x0Vec(idx) = x0;
    
    % Calculate centerline and wake expansion
    for idx2 = 1: length(xDistVec)
        xDist =  xDistVec(idx2);
        % Wake width
        if xDist <= x0
            deltaVec(idx2,idx) = ThetaCo * xDist;
            sigmaYVec(idx2,idx) = D * cos(gamma)/sqrt(8);
        else
            
            % wake defelection delta
            term1 = ThetaCo * x0;
            % 7.2 eq. for sigmaY and sigmaZ
            sigmaY = b * (xDist - x0) + D * cos(gamma)/sqrt(8);
            sigmaZ = b * (xDist - x0) + D/sqrt(8);
            
            sqrtCt = sqrt(Ct1);
            tempTerm = 1.6 * sqrt((8*sigmaY*sigmaZ)/(D^2*cos(gamma)));
            
            term2 = ThetaCo * D/14.7 * sqrt(cos(gamma)/(b^2*Ct1)) *...
                (2.9 + 1.3 *sqrtCt - Ct1) * ...
                log( (1.6 + sqrtCt)*(tempTerm - sqrtCt)/...
                ((1.6 - sqrtCt)*(tempTerm + sqrtCt)));
            
            deltaVec(idx2,idx) = term1 + term2;
            sigmaYVec(idx2,idx) = sigmaY;
            sigmaZVec(idx2,idx) = sigmaZ;
        end
    end
    
    % Get wake deflection and expansion at second turbine
    xDist = tempXDist;
    [~,idxD] = min(abs(xDistVec - xDist));
    d = deltaVec(idxD,idx);
    
    sigmaY = sigmaYVec(idxD,idx); % wake widths at xDist;
    sigmaZ = sigmaZVec(idxD,idx);
    
     %Eq. 7.1 far wake turbine deficit: sigma
    fracCT = (Ct1*cos(gamma)) / (8*(sigmaY*sigmaZ/D^2));
    
    yvec = -r : r; %d - R : d + R;
    zvec = -r : r;
    dUVec = nan(length(yvec),length(zvec));
    for idxZ = 1: length(zvec)
        zDist = zvec(idxZ);
        for idxY = 1: length(yvec)
            yDist = yvec(idxY);
            if (zDist^2 + yDist^2 < r^2)
                dUVec(idxY,idxZ) = (1- sqrt(1 - fracCT ))*exp(-0.5*( ((yDist - d)/sigmaY)^2 + (zDist/sigmaZ)^2));
            end
        end
    end
   
    dUvecCell{idx} = dUVec;
    
    dURWake = mean(dUVec(:),'omitnan');
    dWT2 = (1 - dURWake);
    
    deltaV = [1, dWT2]; % deltaV =  [0,(c12* a1 * A12ToA1Vec)];
    U = Vinf .* deltaV; %(1 - 2 * deltaV).* %[cos(gamma), 1]; %disp(V) ToDo: This is only for two turbines
    Cp1 = 1; %interp1(cpVec(1,:),cpVec(2,:),gammaDeg)/cpVec(2,1);

    a12 = [a1,a2];
    gamma12 = [gamma,0];
    Cp = 4*a12.*(sqrt(1- a12.*(2*cos(gamma12)-a12))).^2; %Cp = 4*[a1,a2].*([sqrt(1- a1*(2*cos(gamma)-a1)),(1- a2)]).^2;
    Cp = [Cp(1)*Cp1, Cp(2)]; %*0.6
    
    P(idx,:) = 1/2* rho * A* Cp .* (U .* [cos(gamma), 1]).^3;%  - [0, 5.4798e+05];
    Vvec(idx,:) = U;
    ThetaCoVec(idx) = ThetaCo;
end

% Comparison WFSim/ Gaussian Wake model
Psum = sum(P,2);
[Pmax,idxPmax]= max(Psum);

%% Figure PowerVsYawWFSim and Gaussian Bastankah Wake model
nF = figure;
pos = get(nF,'Position');
set(nF,'Position',[pos(1:3),pos(4)*0.75]);

plot(gammaVec,P(:,1)/10^6,'bs--',gammaVec,P(:,2)/10^6,'rs--',...
    gammaVec,Psum/10^6,'ks--',gammaVec(idxPmax),Pmax/10^6,'gs')
legend('P WT1','P WT2', 'P T','max(P T)','Location','eastoutside', 'Interpreter','Latex','FontSize',fs);
axis tight;  posAxis = axis; axis([posAxis(1:2),0,posAxis(4)*1.01]);

set(gca,'TickLabelInterpreter','Latex','FontSize',fs)
xlabel('yaw $\gamma_1$ [deg]','Interpreter','Latex','FontSize',fs);
ylabel('Power [MW]','Interpreter','Latex','FontSize',fs);



%sprintf('max P: $\\gamma$ = %d [deg]', phi1K(vecIdx)),...
strName = sprintf('PowerVsYaw');
set(gcf,'Name', strName);
% print(gcf,fullfile(dirFig,strName), '-dpng');
% print(gcf,fullfile(dirFig,strName), '-depsc');



set(0,'DefaultAxesFontSize',dAFS);
set(0,'DefaultTextFontSize',dTFS);
set(0,'DefaultLineLineWidth',dLLw);
