function [] = compareImgs(img1, img2)
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc(img1)
    colormap gray
    figure('units','normalized','outerposition',[0 0 1 1])
    imagesc(img2)
    colormap gray
    pause
end
