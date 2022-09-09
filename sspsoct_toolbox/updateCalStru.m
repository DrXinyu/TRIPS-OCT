function CalStru = updateCalStru(k_poly,dispersion_poly,CalStru)

    if isfield(CalStru,'ShowChange')
        if CalStru.ShowChange
            mad_fit = polyval(k_poly,CalStru.k_fitting.x);           
            figure("name","k change")
            
            MAmean = cumsum([0;mad_fit]);
            MAmean = MAmean - (MAmean(1000) - CalStru.MAmean(1000));
            plot(MAmean)

            hold on
            
            old_mad_fit = polyval(CalStru.k_fitting.poly,CalStru.k_fitting.x);        
            MAmean = cumsum([0;old_mad_fit]);
            MAmean = MAmean - (MAmean(1000) - CalStru.MAmean(1000));
            plot(MAmean)

            legend("updated","old");
            figure("name","dispersion change")
            plot(polyval(dispersion_poly,CalStru.dispersion_fitting.x));
            hold on
            plot(polyval(CalStru.dispersion_fitting.poly,CalStru.dispersion_fitting.x));
            legend("updated","old");
            
            

        else
        end
    else
    end

    mad_fit = polyval(k_poly,CalStru.k_fitting.x);
    
    
    mad_fit(mad_fit<=0)= -mad_fit(mad_fit<=0);
    %fit_mask = (mad_fit>=0);

    MAmean = cumsum([0;mad_fit]);
    MAmean = MAmean - (MAmean(1000) - CalStru.MAmean(1000));
    
    xi = (linspace(MAmean(1),MAmean(end),4096))';
 
    CalStru.MAmean = MAmean;
    CalStru.xi = xi;
    
    dp_fit = polyval(dispersion_poly,CalStru.dispersion_fitting.x);
    p_dpl = polyfit((0:(length(dp_fit))-1)',dp_fit,1);
    dp_fit = dp_fit - polyval(p_dpl,(0:(length(dp_fit)-1))');
    
    CArray = (exp(-(1i.*dp_fit)));

    CalStru.CArray = CArray;
    
    CalStru.k_fitting.poly = k_poly;
    CalStru.dispersion_fitting.poly = dispersion_poly;    
    
end

