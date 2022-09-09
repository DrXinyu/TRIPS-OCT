function [stablized_fringeA,stablized_fringeB,CalStru] = k_stabilize(fringeA_o,fringeB_o,CalStru,stabilize)
%
%   This function uses an auxiliary signal to stablize the fringe.
%   frq is the normalized frequency of axuiliary signal ranging from 0-1.
%   by Xinyu Liu

%% pre-interpolation    
    
    oldCalStru = CalStru;
    if isfield(CalStru,'GPU')
        if CalStru.GPU == 1
            fringeA_o = gpuArray(single(fringeA_o)); 
            fringeB_o = gpuArray(single(fringeB_o));    
        else  
        end
    else
    end
    
    fringeA_o = fringeA_o(CalStru.Ang_Cutting_Limit(1):CalStru.Ang_Cutting_Limit(2),:);
    fringeB_o = fringeB_o(CalStru.Ang_Cutting_Limit(1):CalStru.Ang_Cutting_Limit(2),:);
    
    fringeA_o = interp1(CalStru.MAmean,fringeA_o,CalStru.xi);
    fringeB_o = interp1(CalStru.MAmean,fringeB_o,CalStru.xi); 
    
    
    fringeA = fringeA_o-smooth(mean(fringeA_o,2),50);
    fringeB = fringeB_o-smooth(mean(fringeB_o,2),50);
    

%     fringeA = double(gather(fringeA));
%     fringeB = double(gather(fringeB));
%     
    
    CalStru.Background = 'smooth';
    CalStru.Dispersion = 0;
    CalStru.K_stabilized = 1;
    CalStru.NFFT = 4096*16;
    CalStru.Complex = 1;
    CalStru.MinDepth = 0;
    CalStru.Shifting = 0;
    CalStru.Binning = 0;

    if isfield(CalStru,'GPU')
        if CalStru.GPU == 1
            fringeA_o = gpuArray(single(fringeA_o));
            fringeB_o = gpuArray(single(fringeB_o));  
        else
        end
    else
    end    

    
%% re-interpolation on scaling


    [pixnum,ascannum] = size(fringeA_o);
    fringeA_o = interp1(CalStru.MAmean_fit2,fringeA_o,CalStru.xi2,'linear','extrap');
    fringeB_o = interp1(CalStru.MAmean_fit2,fringeB_o,CalStru.xi2,'linear','extrap');
    
    stablized_fringeA = fringeA_o;
    stablized_fringeB = fringeB_o;
    
    if stabilize
    
        CalStru.Background = 'smooth';
        CalStru.Dispersion = 0;
        CalStru.K_stabilized = 1;
        CalStru.NFFT = 4096*4;
        CalStru.Complex = 1;
        CalStru.MinDepth = 0;
        CalStru.Shifting = 0;

        auxiliary_confidence = round((CalStru.NFFT/2)*(CalStru.afrq-0.02)):round((CalStru.NFFT/2)*(CalStru.afrq+0.02));
        
        CA = fringe2image(gather(fringeA_o),CalStru);
        IoA = abs(CA);

        [~,pLA] = max(IoA(auxiliary_confidence,:),[],1);
        pLA = pLA+auxiliary_confidence(1)-1;

        meanpLA = median(pLA);
        ind = round(CalStru.meanpLA);
        scaleA = pLA./meanpLA;

        phaseA = rotate_unwrap([angle(CA(ind,:)),0]);
        phaseA = unwrap(phaseA-phaseA(end));


        CB = fringe2image(gather(fringeB_o),CalStru);
        IoB = abs(CB);
        [~,pLB] = max(IoB(auxiliary_confidence,:),[],1);
        pLB = pLB+auxiliary_confidence(1)-1;

        meanpLB = median(pLB);
        ind = round(CalStru.meanpLB);
        scaleB = pLB./meanpLB;

        phaseB = rotate_unwrap([angle(CB(ind,:)),0]);
        phaseB = unwrap(phaseB-phaseB(end));

        scale = (scaleA+scaleB)/2;
        scale = scale'-smooth(scale,100);
        
        scale(scale<0.998) =0.998;
        scale(scale>1.002) = 1.002;
        
        phase = (phaseA(1:end-1)+phaseB(1:end-1))/2;
        phase = phase'-smooth(phase,100);

        oldCalStru.test_phase = phase;
        
        MAmean = CalStru.xi2;
        xi = MAmean;

        MA_M = repmat(MAmean,1,ascannum);
        MA_M = bsxfun(@plus,MA_M,phase');  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scale
        %MA_M = bsxfun(@times,MA_M,scale'); 
        %MA_B = bsxfun(@multiple,MAmean,scaleB);    

        for index = 1:ascannum
            stablized_fringeA(:,index) = interp1(MA_M(:,index),fringeA_o(:,index),xi,'linear','extrap');
            stablized_fringeB(:,index) = interp1(MA_M(:,index),fringeB_o(:,index),xi,'linear','extrap');   

        end

    end
    CalStru = oldCalStru;
    CalStru.K_stabilized = 1;
    
    %%
    

end
