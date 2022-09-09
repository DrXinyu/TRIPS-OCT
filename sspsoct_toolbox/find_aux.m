function newCalStru = find_aux(fringeA,fringeB,CalStru)
    % locate the auxiliary signal peak
    
    oldCalStru = CalStru;
    
    if isfield(CalStru,'afrq')
        afrq = CalStru.afrq;  
    else
        afrq = 0.15;
    end    
    
    CalStru.Background = 'mean';
    CalStru.Dispersion = 0;
    CalStru.Shiftw = afrq*pi;
    CalStru.Shifting = 0;
    CalStru.GPU = 0;
    CalStru.NFFT = 4096*16;
    CalStru.Binning = 0;
    CalStru.MinDepth = 0;
    I1 = fringe2image(fringeA,CalStru);
    I2 = fringe2image(fringeB,CalStru);

    I = I1+I2;

    auxiliary_confidence = round((CalStru.NFFT/2)*(afrq-0.1)):round((CalStru.NFFT/2)*(afrq+0.1));
    [~,auxiliary_location] = max(mean(I(auxiliary_confidence,:),2));
    afrq = gather(auxiliary_location-1+(CalStru.NFFT/2)*(afrq-0.05))./(CalStru.NFFT/2);
    a_filter = setup_filter(afrq);
    
    newCalStru = oldCalStru;
    
    newCalStru.afrq = afrq;
    newCalStru.a_filter = a_filter;
    

end

