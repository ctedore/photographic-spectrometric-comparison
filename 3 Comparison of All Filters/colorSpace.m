function [chroma, coords] = colorSpace(q)
    sumq = sum(q);
    q = q / sumq;
    % Pike et al 2012 
    X = (1/sqrt(2)) * (q(1) - q(2));
    if numel(q) > 2
        Y = (sqrt(2) / sqrt(3)) * (q(3) - ((q(1) + q(2)) / 2));
    end
    if numel(q) > 3
        Z = (sqrt(3) / 2) * (q(4) - ((q(1) + q(2) + q(3)) / 3));
    end
    if numel(q) == 2
        coords = X;
        chroma = sqrt(X.^2); 
    elseif numel(q) == 3
        coords = [X Y];
        chroma = sqrt(X.^2 + Y.^2); 
    elseif numel(q) == 4
        coords = [X Y Z];
        chroma = sqrt(X.^2 + Y.^2 + Z.^2); 
    end
end