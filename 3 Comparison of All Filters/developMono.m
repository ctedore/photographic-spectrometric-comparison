function [] = developMono(filterType, data, indices)
    incr = 1;
    for k = indices
        filename = sprintf(['%d_', filterType, '.png'], incr);
        imwrite(data{k}, filename, 'png');
        incr = incr + 1;
    end
end