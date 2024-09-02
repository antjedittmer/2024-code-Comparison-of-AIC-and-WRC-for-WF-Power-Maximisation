function [dPT,dPTVec] = calculate_dPT_sum(aVec,c)

m = length(aVec);
dPTVec = nan(m);
isDiff = 1;
aVecN = [aVec; 1/3];
for idx = 1:m
    a = aVec(idx);
    dPTVec(idx) = CpFromA(a,isDiff);
    tmp = 0;
    for idxP = idx:m
        aW = aVecN(idxP+1);
        tmp = tmp + CpFromA(aW) * WFromAandC(c,aVecN(idxP:end),isDiff);
    end
    dPTVec(idx) = dPTVec(idx) + tmp;
end

dPT = sum(dPTVec);