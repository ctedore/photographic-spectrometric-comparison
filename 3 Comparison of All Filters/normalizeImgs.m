function allData = normalizeImgs(allData, indices_for_maxVals)

    % find maxes in each relevant channel
    max_vals = zeros(size(indices_for_maxVals, 2), 1);
    iter = 1;
    for k = indices_for_maxVals
        max_vals(iter) = max(max(allData{k}));
        iter = iter + 1;
    end
    max_E_for_norm = max(max_vals);
    
    % normalize relevant channels
    for k = indices_for_maxVals
        allData{k} = allData{k} / max_E_for_norm;    
    end

end