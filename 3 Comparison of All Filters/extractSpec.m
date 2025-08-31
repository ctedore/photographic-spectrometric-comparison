function spectrum = extractSpec(fID,w,convert)
    textscan(fID,'%f %f','HeaderLines',100);
    data = textscan(fID,'%f %f',740);
    if convert == 1
        data{1,2} = 5.05E15 * data{1,1} .* data{1,2} * 0.01; % converts from uW/cm^2 to photons/m^2 
    end
    spectrum = interp1(data{1,1}, data{1,2}, w);
end