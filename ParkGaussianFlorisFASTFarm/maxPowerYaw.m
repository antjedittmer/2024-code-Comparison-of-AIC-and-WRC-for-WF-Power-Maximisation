function [Psum,P] = maxPowerYaw(gammaVec)
if ~nargin
    gammaVec = [];
end

%% Define parameters
alphaStar =  2.32; %32; % 4*alpha = 4* 0.58 = 2.32 2.8
betaStar = 0.154; % 2*beta = 0.154 from jetstreams
b   = 0.022; % 0.075; % model parameters wake interaction 0.022
I0 = 0.05; %Wp.site.v_Inf/Wp.site.u_Inf; % Intensity: 0 for our case

rho = 1.225; % air density
D   = 126;  % rotor diameter
r   = D/2; %rotor radius
A   = r^2 *pi;% rotor area
n  = length(gammaVec); % number upstream turbines

Vinf  =  8;% Wp.site.u_Inf; % wind unaffected by wind farm
xDistWT = 5*D*(1:n); %diff(Wp.turbine.Crx); % 5*D; % Distance between turbines
xDistWTinv = fliplr(xDistWT);
%c12 = ( D/ (D + 2* b * xDist ) )^2;

a2 = 1/3*0.95; % Induction factor WT2 for maximal power  1/9; %
% Ct2 = 4*a2 *(1- a2); % eq. 6.2 thrust coefficient for max. power

gammaVec1 = (0:5:35);
a1Vec = linspace(1/3, 1/3.5, length(gammaVec1))*0.95; %


if isempty(gammaVec)
    Cp = 4*a2*(1- a2)^2;
    P = 1/2* rho * A* Cp * Vinf^3;
    Psum = P;
    return;
end

% For centerline and wake expansion
nD = ceil(5*D*(n+1));
xDistVec = 0 : nD;
deltaVec = nan(nD,length(gammaVec));
sigmaYVec = nan(nD,length(gammaVec));
sigmaZVec = D/sqrt(8)*ones(nD,length(gammaVec));
%x0Vec = nan(length(gammaVec),1);

% b = 0.022; % model parameters wake expansion 0.022 (from paper), 0.075
% commonly used
dUvecCell = cell(n,1);
dURWakeMatrix = zeros(n,1);

a1WT = nan(n,1);

gammaRad = gammaVec/180*pi;

for idx = 1 : n

    gammaDeg = gammaVec(idx); % current yaw angle WT1
    a1 = interp1(gammaVec1,a1Vec,gammaDeg);
    gamma =  gammaRad(idx);

    % Calculate variables based on paper
    % eq.6.7 Thrust coefficient (approximation of eq. 6.1
    % Ct1 = 4*a1 *(1- a1*cos(gamma)); %* interp1(ctVec(1,:),ctVec(2,:),gammaDeg);
    Ct1 = 4*a1 * sqrt(1- a1*(2*cos(gamma)-a1)); %* interp1(ctVec(1,:),ctVec(2,:),gammaDeg)/ctVec(2,1);
    % eq. 6.12: Wake angle at rotor 30 -> 6 deg
    ThetaCo = - 0.3*gamma/cos(gamma) *(1- sqrt(1-Ct1*cos(gamma)));

    % Question phi = - sin(gamma)*Ct1/(2*(1+ b * xDist/D)^2): 14 deg >> 6 deg
    idx1 = interp1(gammaVec1,1 : length(gammaVec1),gammaDeg);
    I = I0 +   0.0150 *(idx1-1); %%abs(sol_array(kPlotVec(idx)).v(nTx,nTy)/sol_array(kPlotVec(idx)).u(nTx,nTy));
    % end of near wake area x0
    sqrt1minCt = sqrt(1-Ct1);
    x0 = cos(gamma)* (1 + sqrt1minCt) / ...
        (sqrt(2)*(alphaStar*I + betaStar*(1 - sqrt1minCt))) * D; % eq. 6.16/ 7.3
    %x0Vec(idx) = x0;

    % Calculate centerline and wake expansion
    anXDistVec = xDistVec(1:5*D*(n+2-idx));
    for idx2 = 1: length(anXDistVec)
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

    % Get wake deflection and expansion at downwind turbines
    % idxT = idx;
    xWTDist = xDistWTinv(n); %only for the next turbine

    %for idx3 = 1: length(xWTDist)
    anXDistT = xWTDist; % distance up- to downwind turbine
    [~,idxD] = min(abs(anXDistVec -  anXDistT));
    d = deltaVec(idxD,idx);

    sigmaY  = sigmaYVec(idxD,idx); % wake widths at xDist;
    % R = sigmaY;% wake widths at xDist;
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
                dUVec(idxY,idxZ) = (1- sqrt(1 - fracCT ))*exp(-0.5*(((yDist - d)/sigmaY)^2 + (zDist/sigmaZ)^2));
            end
        end
    end

    dUvecCell{idx} = dUVec;
    dURWakeMatrix(idx) = mean(dUVec(:),'omitnan');
    %end
    a1WT(idx) = a1;
end

dURWake =  dURWakeMatrix';
dWT2 = cumprod(1 - dURWake);

deltaV = [1, dWT2]; % deltaV =  [0,(c12* a1 * A12ToA1Vec)];
U = Vinf .* deltaV; %(1 - 2 * deltaV).* %[cos(gamma), 1]; %disp(V) ToDo: This is only for two turbines
Cp1 = ones(n,1); %interp1(cpVec(1,:),cpVec(2,:),gammaDeg)/cpVec(2,1);
Cp = 4*[a1WT;a2].*([sqrt(1- a1WT.*(2*cos(gammaRad)-a1WT));(1- a2)]).^2;
Cp = [Cp(1:n).*Cp1; Cp(n+1)]; %*0.6

P = 1/2* rho * A* Cp .* (U' .* [cos(gammaRad); 1]).^3;%  - [0, 5.4798e+05];


% Comparison WFSim/ Gaussian Wake model
Psum = - sum(P);
%[Pmax,idxPmax] = max(Psum);


