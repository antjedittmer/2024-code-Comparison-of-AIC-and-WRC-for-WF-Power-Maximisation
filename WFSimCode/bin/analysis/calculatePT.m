function [PT, P] = calculatePT(gammaVec,a1Vec,xDistT,Vinf,D)

% [gammaOpt,PTopt]  = fminsearch(@(gammaVec) ...
% calculatePT(gammaVec,a1Vec,xDistT,Vinf,D), gammaVec)

if ~nargin
    Vinf = 8;
    gammaVec = [0,0,0];
    a1Vec = [1/3,1/3,1/3];
    D = 126; 
    xDistT = 5*D*(0:(length( a1Vec)-1)); % Distance turbines in that row
end
% nD = ceil(max(xDistT));
%xDistVec = 0 : nD;
alphaStar = 2.32; %32; % 4*alpha = 4* 0.58 = 2.32 2.8
betaStar = 0.154; % 2*beta = 0.154 from jetstreams
I0 = 0.05; % Wp.site.v_Inf/Wp.site.u_Inf; % Intensity: 0 for our case
b = 0.022; % 0.075; % model parameters wake interaction 0.022
r = D/2;

gammaCorVec = (0:35)/180 *pi; %gamma in radian
aCorVec = linspace(1/3, 1/3.5, length(gammaCorVec)) / (1/3);

a1 = a1Vec*0.95 .* interp1(gammaCorVec,aCorVec,gammaVec); % reduction for all induction factors

% Calculate variables based on paper
% eq.6.1 thrust coefficient (approximated in equation 6.7) 
Ct1 = 4*a1 .* sqrt(1 - a1.*(2*cos(gammaVec)-a1)); 
% eq. 6.12: Wake angle at rotor 30 -> 6 deg
ThetaCo = - 0.3*gammaVec./cos(gammaVec) .*(1-sqrt(1-Ct1.*cos(gammaVec)));

% Question phi = - sin(gamma)*Ct1/(2*(1+ b * xDist/D)^2): 14 deg >> 6 deg
Ivec = I0 + 0.0030*gammaVec *180/pi;%0.0150 *(idx-1); %%abs(sol_array(kPlotVec(idx)).v(nTx,nTy)/sol_array(kPlotVec(idx)).u(nTx,nTy));

% end of near wake area x0
sqrt1minCt = sqrt(1-Ct1);
x0Vec = cos(gammaVec).* (1 + sqrt1minCt) ./ ...
    (sqrt(2)*(alphaStar.*Ivec + betaStar.*(1 - sqrt1minCt))) * D; % eq. 6.16/ 7.3

noW = (length(gammaVec)-1);
deltaCell = cell(noW,1);
sigmaYCell = cell(noW,1);
sigmaZCell = cell(noW,1);

% Calculate centerline and wake expansion for all downwind turbines
for idxT = 1: noW

    % For all turbines calculate the influence on the next turbine below
    % First calculate the 
    xDistVec = 0: xDistT(idxT+1) -  xDistT(idxT);

    x0 = x0Vec(idxT);
    for idx2 = 1: length(xDistVec) % calculate the speed over the 
        xDist =  xDistVec(idx2);
        % Wake width
        if xDist <= x0
            deltaCell{idxT}(idx2) = ThetaCo(idxT) * xDist;
            sigmaYCell{idxT}(idx2) = D * cos(gammaVec(idxT))/sqrt(8);
            sigmaZCell{idxT}(idx2) = D /sqrt(8);
        else

            % wake defelection delta
            term1 = ThetaCo(idxT) * x0;
            % 7.2 eq. for sigmaY and sigmaZ
            sigmaY = b * (xDist - x0) + D * cos(gammaVec(idxT))/sqrt(8);
            sigmaZ = b * (xDist - x0) + D/sqrt(8);

            sqrtCt = sqrt(Ct1(idxT));
            tempTerm = 1.6 * sqrt((8*sigmaY*sigmaZ)/(D^2*cos(gammaVec(idxT))));

            term2 = ThetaCo(idxT) .* D/14.7 .* sqrt(cos(gammaVec(idxT))./(b^2*Ct1(idxT))) .*...
                (2.9 + 1.3 .*sqrtCt - Ct1(idxT)) .* ...
                log( (1.6 + sqrtCt).*(tempTerm - sqrtCt)./...
                ((1.6 - sqrtCt).*(tempTerm + sqrtCt)));

            deltaCell{idxT}(idx2) = term1 + term2;
            sigmaYCell{idxT}(idx2) = sigmaY;
            sigmaZCell{idxT}(idx2) = sigmaZ;
        end
    end

    % Get wake deflection and expansion for the next downwind turbines
    d = deltaCell{idxT}(end);

    sigmaY  = sigmaYCell{idxT}(end); % wake widths at xDist;
    %R = sigmaY;% wake widths at xDist;
    sigmaZ = sigmaZCell{idxT}(end); 

    %Eq. 7.1 far wake turbine deficit: sigma
    fracCT = (Ct1(idxT)*cos(gammaVec(idxT))) / (8*(sigmaY*sigmaZ/D^2));

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

    dUvecCell{idxT} = dUVec;

    dURWake = mean(dUVec(:),'omitnan');
    dWT2 = (1 - dURWake);

    deltaV(idxT) = dWT2; % deltaV =  [0,(c12* a1 * A12ToA1Vec)];
end

U = Vinf .* cumprod([1,deltaV]); %(1 - 2 * deltaV).* %[cos(gamma), 1]; %disp(V) ToDo: This is only for two turbines
Cp1 = 1; %interp1(cpVec(1,:),cpVec(2,:),gammaDeg)/cpVec(2,1);
Cp = 4*a1 .* (sqrt(1- a1.*(2*cos(gammaVec)-a1)).^2);
Cp = [Cp(1)*Cp1, Cp(2:end)]; %*0.6

A = r^2 *pi;
rho = 1.225;
P = 1/2* rho * A* Cp .* (U .* cos(gammaVec)).^3;

PT = - sum(P); % negative total sum for mini