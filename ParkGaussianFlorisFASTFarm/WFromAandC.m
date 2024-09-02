function W = WFromAandC(c,aVec,isDiff)

if nargin < 3
    isDiff = [];
end
 
rowN = length(aVec); %which lin
c_kTo1 = c(rowN,1:rowN);% C lines sel

sqrt1 = sqrt(sum((c_kTo1.*aVec').^2));

if isempty(isDiff)
    W = (1 - 2 * sqrt1)^3; %disp(V)
elseif sqrt1 == 0
    W = 1;
else
    W = 3*((1 - 2 * sqrt1)^2)*(-1/sqrt1)*2*c_kTo1(1)^2*aVec(1);
end
