function power = objective_fun(k_fitting,dispersion_fitting,fringe,CalStru,ROI)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    CalStru = updateCalStru(k_fitting,dispersion_fitting,CalStru);

    image = fringe2image(fringe,CalStru);
    
    power = -max(image(ROI,:));
    

end

