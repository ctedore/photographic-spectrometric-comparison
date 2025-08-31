function deltaS = euclidDist(patch1, patch2)
    Xdiff = patch1(1) - patch2(1);
    deltaS = sqrt(Xdiff.^2);
    if length(patch1) == 2
        Ydiff = patch1(2) - patch2(2);
        deltaS = sqrt(Xdiff.^2 + Ydiff.^2);
    elseif length(patch1) == 3
        Ydiff = patch1(2) - patch2(2);
        Zdiff = patch1(3) - patch2(3);
        deltaS = sqrt(Xdiff.^2 + Ydiff.^2 + Zdiff.^2);
    end
end