function [stablized_fringeA,stablized_fringeB,CalStru] = k_stabilize_setup(fringeA_o,fringeB_o,CalStru)
%
%   This function uses an auxiliary signal to stablize the fringe.
%   frq is the normalized frequency of axuiliary signal ranging from 0-1.
%   by Xinyu Liu

    %% pre-interpolation    
    
    oldCalStru = CalStru;
    if isfield(CalStru,'GPU')
        if CalStru.GPU == 1
            %auxiliaryA = gpuArray(single(auxiliaryA));
            fringeA_o = gpuArray(single(fringeA_o)); 
            %auxiliaryB = gpuArray(single(auxiliaryB));
            fringeB_o = gpuArray(single(fringeB_o));    
        else
            
        end
    else
        
    end
    
    fringeA_o = fringeA_o(CalStru.Ang_Cutting_Limit(1):CalStru.Ang_Cutting_Limit(2),:);
    fringeB_o = fringeB_o(CalStru.Ang_Cutting_Limit(1):CalStru.Ang_Cutting_Limit(2),:);
    
    fringeA_o = interp1(CalStru.MAmean,fringeA_o,CalStru.xi);
    fringeB_o = interp1(CalStru.MAmean,fringeB_o,CalStru.xi);    

    
    padding = mean(fringeA_o(1,:));
    
    
    fringeA = fringeA_o-smooth(mean(fringeA_o,2),50);
    fringeB = fringeB_o-smooth(mean(fringeB_o,2),50);
    

    fringeA = double(gather(fringeA));
    fringeB = double(gather(fringeB));
    

    auxiliaryA = kFilter(fringeA,CalStru.a_filter);
    auxiliaryB = kFilter(fringeB,CalStru.a_filter);
    
    CalStru.Background = 'smooth';
    CalStru.Dispersion = 0;
    CalStru.K_stabilized = 1;
    CalStru.NFFT = 4096*16;
    CalStru.Complex = 1;
    CalStru.MinDepth = 0;
    CalStru.Shifting = 0;
    CalStru.Binning = 0;
%     auxiliaryA = reshape_fringe(auxiliaryA);
%     auxiliaryB = reshape_fringe(auxiliaryB);
    
    if isfield(CalStru,'GPU')
        if CalStru.GPU == 1
            auxiliaryA = gpuArray(single(auxiliaryA));
            fringeA_o = gpuArray(single(fringeA_o));
            auxiliaryB = gpuArray(single(auxiliaryB));
            fringeB_o = gpuArray(single(fringeB_o));  
        else
        end
    else
    end    

    
%% re-interpolation on scaling

    [pixnum,ascannum] = size(auxiliaryA);
    HyA = hilbert(auxiliaryA(500:end-500,:));
    HyB = hilbert(auxiliaryB(500:end-500,:));
    
    [~,T1] = (max(mean(abs(HyA)+abs(HyB),2)));
    T1 = T1+500;
    
    %T1 = round((pixnum/3)*2);

    Hy = hilbert(auxiliaryA);
    angA = angle(Hy);
    HyAng = force_increasing_unwrap(unwrap(angA)); 
    HyAng = bsxfun(@minus,HyAng,HyAng(T1,:));
    HyAngA = mean(HyAng,2);
    
    Hy = hilbert(auxiliaryB);
    angB = angle(Hy);
    HyAng = force_increasing_unwrap(unwrap(angB)); 

    HyAng = bsxfun(@minus,HyAng,HyAng(T1,:));
    HyAngB = mean(HyAng,2);
    MA_auxiliary = (HyAngA+HyAngB)/2;
    
    
%     
    fa = [5,pixnum-5];
    xo = (fa(1):fa(2)-1)';
    normxo = (xo-mean(xo))/std(xo);
    mad = gather(diff(MA_auxiliary(fa(1):fa(2))));

    k_fit_x = ((1:((length(MA_auxiliary))-1))'-mean(xo))/std(xo);
    k_poly = polyfit(normxo,mad,7);

    mad_fit = polyval(k_poly,k_fit_x);
    mad_fit(mad_fit<=0)= 0.0001;
% 
    MAmean_fit = cumsum([MA_auxiliary(1);mad_fit]);    
    
    xi = (linspace((min(MAmean_fit)),(max(MAmean_fit)),pixnum))';
    CalStru.MAmean_fit2 = MAmean_fit;
    CalStru.xi2 = xi;

    
    fringeA_o = interp1(CalStru.MAmean_fit2,fringeA_o,CalStru.xi2,'linear','extrap');
    fringeB_o = interp1(CalStru.MAmean_fit2,fringeB_o,CalStru.xi2,'linear','extrap');
%     auxiliaryA = interp1(MAmean_fit,auxiliaryA,xi,'linear','extrap');
%     auxiliaryB = interp1(MAmean_fit,auxiliaryB,xi,'linear','extrap');    
    
    
%     Hy = hilbert(auxiliaryA);
%     angA = angle(Hy);
%     HyAng = force_increasing_unwrap(unwrap(angA)); 
%     HyAng = bsxfun(@minus,HyAng,HyAng(T1,:));
%     HyAngA = mean(HyAng,2);
%     
%     Hy = hilbert(auxiliaryB);
%     angB = angle(Hy);
%     HyAng = force_increasing_unwrap(unwrap(angB)); 
% 
%     HyAng = bsxfun(@minus,HyAng,HyAng(T1,:));
%     HyAngB = mean(HyAng,2);
%     MA_auxiliary = (HyAngA+HyAngB)/2;
%     xi = linspace((min(MA_auxiliary)),(max(MA_auxiliary)),pixnum);
    
    
    
    CalStru.Background = 'smooth';
    CalStru.Dispersion = 0;
    CalStru.K_stabilized = 1;
    CalStru.NFFT = 4096*4;
    CalStru.Complex = 1;
    CalStru.MinDepth = 0;
    CalStru.Shifting = 0;
    
    auxiliary_confidence = round((CalStru.NFFT/2)*(CalStru.afrq-0.02)):round((CalStru.NFFT/2)*(CalStru.afrq+0.02));

    CA = fringe2image(fringeA_o,CalStru);
    IoA = abs(CA);

    [~,pLA] = max(IoA(auxiliary_confidence,:),[],1);
    pLA = pLA+auxiliary_confidence(1)-1;
    
    meanpLA = median(pLA);
    ind = round(meanpLA);
    scaleA = pLA./meanpLA;
    
    phaseA = rotate_unwrap([angle(CA(ind,:)),0]);
    phaseA0 = mean(phaseA);

    phaseA = unwrap(phaseA-phaseA(end));


    CB = fringe2image(fringeB_o,CalStru);
    IoB = abs(CB);
    [~,pLB] = max(IoB(auxiliary_confidence,:),[],1);
    pLB = pLB+auxiliary_confidence(1)-1;
    meanpLB = median(pLB);
    ind = round(meanpLB);
    scaleB = pLB./meanpLB;

    phaseB = rotate_unwrap([angle(CB(ind,:)),0]);    
    phaseB0 = mean(phaseB);
    phaseB = unwrap(phaseB-phaseB(end));
    

    scale = (scaleA+scaleB)/2;
    scale = scale'-smooth(scale,40);
    scale(scale<0.998) =0.998;
    scale(scale>1.002) = 1.002;
    
    phase = (phaseA(1:end-1)+phaseB(1:end-1))/2;
    phase = phase'-smooth(phase,40);
    
    MAmean = xi;
    xi = MAmean;

    MA_M = repmat(MAmean,1,ascannum);
    MA_M = bsxfun(@plus,MA_M,phase');    
    %MA_M = MA_M+phase;   
    MA_M = bsxfun(@times,MA_M,scale'); 
    %MA_B = bsxfun(@multiple,MAmean,scaleB);    

    for index = 1:ascannum
        stablized_fringeA(:,index) = interp1(MA_M(:,index),fringeA_o(:,index),xi,'linear','extrap');
        stablized_fringeB(:,index) = interp1(MA_M(:,index),fringeB_o(:,index),xi,'linear','extrap');   
        
%         stablized_auxiliaryA(:,index) = interp1(MA_M(:,index),auxiliaryA(:,index),xi,'linear');
%         stablized_auxiliaryB(:,index) = interp1(MA_M(:,index),auxiliaryB(:,index),xi,'linear');
    end

    
    CalStru = oldCalStru;
    
    CalStru.K_stabilized = 1;
    CalStru.MAmean_fit2 = MAmean_fit;
    CalStru.xi2 = xi;
    CalStru.meanpLA = meanpLA;
    CalStru.meanpLB = meanpLB;
    
    
    
end

