function cp = CpFromA(a,isDiff)
% calculated power coefficient cp from axial induction a
if nargin < 2
    isDiff = [];
end

if isempty(isDiff)
    cp = 4*a*(1-a)^2;
else
    cp = (1-a)*(4-12*a); % 4*(1-a)^2 - 8*a*(1-a) 
end