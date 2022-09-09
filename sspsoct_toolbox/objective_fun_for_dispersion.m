function power = objective_fun_for_dispersion(k_fitting,dispersion_fitting,fringe,CalStru,ROI)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    CalStru.ShowChange = 0;
    CalStru = updateCalStru(k_fitting,dispersion_fitting,CalStru);

    image = gather(fringe2image(fringe,CalStru));
%     figure(1)
%     imagesc(log(image(ROI{1},:)));
    
    
    if isfield(CalStru,'Fullrange')
        if CalStru.Fullrange
            image = fftshift(image,1);
        else
        end
    else
    end
    power = 0;
    for ROIindex = 1:length(ROI)
        power = mean(-max(image(ROI{ROIindex},:)))+power;
    end
    

end

