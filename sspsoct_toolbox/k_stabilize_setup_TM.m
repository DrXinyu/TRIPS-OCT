function [stablized_fringeA,stablized_fringeB,CalStru] = k_stabilize_setup_TM(fringeA_o,fringeB_o,CalStru)
%
%   This function uses an auxiliary signal to stablize the fringe.
%   frq is the normalized frequency of axuiliary signal ranging from 0-1.
%   by Xinyu Liu

    %% pre-interpolation    
    
    oldCalStru = CalStru;
    
    CalStru.GPU = 0;
    if isfield(CalStru,'GPU')
        if CalStru.GPU == 1
            %auxiliaryA = gpuArray(single(auxiliaryA));
            fringeA_o = gpuArray(single(fringeA_o)); 
            %auxiliaryB = gpuArray(single(auxiliaryB));
            fringeB_o = gpuArray(single(fringeB_o));    
        else
            fringeA_o = gather(fringeA_o);
            fringeB_o = gather(fringeB_o);
        end
    else
        
    end
        
    CalStru.Background = 'smooth';
    CalStru.Dispersion = 0;
    CalStru.Shifting = 0;
    CalStru.GPU = 0;
    CalStru.NFFT = 4096*16;
    CalStru.Binning = 0;
    CalStru.MinDepth = 0;
    I1 = fringe2image(fringeA_o,CalStru);
    I2 = fringe2image(fringeB_o,CalStru);

    I = I1+I2;

    auxiliary_confidence = round((CalStru.NFFT/2)*(CalStru.afrq-0.1)):round((CalStru.NFFT/2)*(CalStru.afrq+0.1));
    [~,auxiliary_location] = max(mean(I(auxiliary_confidence,:),2));
    afrq = gather(auxiliary_location-1+(CalStru.NFFT/2)*(CalStru.afrq-0.1))./(CalStru.NFFT/2);
    a_filter = setup_filter(afrq);
    
    
    fringeA_o = fringeA_o(CalStru.Ang_Cutting_Limit(1):CalStru.Ang_Cutting_Limit(2),:);
    fringeB_o = fringeB_o(CalStru.Ang_Cutting_Limit(1):CalStru.Ang_Cutting_Limit(2),:);
    
    fringeA_o = interp1(CalStru.MAmean,fringeA_o,CalStru.xi);
    fringeB_o = interp1(CalStru.MAmean,fringeB_o,CalStru.xi);    
    
    padding = mean(fringeA_o(1,:));
    
    
    fringeA = fringeA_o-smooth(mean(fringeA_o,2),50);
    fringeB = fringeB_o-smooth(mean(fringeB_o,2),50);
    

    fringeA = double(gather(fringeA));
    fringeB = double(gather(fringeB));
    
    auxiliaryA = kFilter(fringeA,a_filter);
    auxiliaryB = kFilter(fringeB,a_filter);
    
    
    
    
    
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
    MA_auxiliary = HyAngA;%(HyAngA+HyAngB)/2;
    
    
%     
    fa = [5,pixnum-5];
    xo = (fa(1):fa(2)-1)';
    normxo = (xo-mean(xo))/std(xo);
    mad = gather(diff(MA_auxiliary(fa(1):fa(2))));

    k_fit_x = ((1:((length(MA_auxiliary))-1))'-mean(xo))/std(xo);
    k_poly = polyfit(normxo,mad,15);

    mad_fit = polyval(k_poly,k_fit_x);
    mad_fit(mad_fit<=0)= 0.0001;
% 
    MAmean_fit = cumsum([MA_auxiliary(1);mad_fit]);    
    
    xi = (linspace((min(MAmean_fit)),(max(MAmean_fit)),pixnum))';
    CalStru.MAmean_fit2 = MAmean_fit;
    CalStru.xi2 = xi;

    
    stablized_fringeA = interp1(CalStru.MAmean_fit2,fringeA_o,CalStru.xi2,'linear','extrap');
    stablized_fringeB = interp1(CalStru.MAmean_fit2,fringeB_o,CalStru.xi2,'linear','extrap');

    
    CalStru = oldCalStru;
    
    CalStru.K_stabilized = 1;
    CalStru.MAmean_fit2 = MAmean_fit;
    CalStru.xi2 = xi;
    CalStru.Shiftw = afrq*pi;
       
    
end

