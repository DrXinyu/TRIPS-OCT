function [colorfulstokes] = stokes2color(stokes)
    [azimuth,elevation] = stokes2polar(stokes);
    
    hsvmap(:,:,1) = (azimuth+pi)./(2*pi);
    saturation = stokes(:,:,3);
    saturation(stokes(:,:,3)<0) = 0;
    saturation = 1-saturation;
    
    value = stokes(:,:,3);
    value(elevation>0) = 0;
    value = 1+value;
    
    hsvmap(:,:,2) = saturation;
    hsvmap(:,:,3) = value;
    
    colorfulstokes = hsv2rgb(hsvmap);
    
end

