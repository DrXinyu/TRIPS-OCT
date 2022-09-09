function [azimuth,elevation] = stokes2polar(stokes)

    azimuth = atan2(stokes(:,:,2),stokes(:,:,1));
    elevation = real(asin(complex(stokes(:,:,3))));

end

