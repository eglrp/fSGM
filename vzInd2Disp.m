function D = vzInd2Disp(w, O, vMax, n)
    vzRatio = w./n*vMax;
    vzInd = vzRatio ./ (1-vzRatio);
    D = O.*vzInd;
end


