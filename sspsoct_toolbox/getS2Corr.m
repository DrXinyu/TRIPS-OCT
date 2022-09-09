function [Acorr,thit] = getS2Corr(fringeA,fringeB,CalStru)
%S1S2
%
    CalStru.Complex = 1;
    CalStru.Binning = 0;
    CalStru.Background = 'smooth';
    CalStru.MinDepth = 1;

    CalStru.Dispersion = 1;

    CalStru.Shifting = 0;
    J1o = fringe2image(fringeA,CalStru);
    J2o = fringe2image(fringeB,CalStru);


    %S1 = squeeze(S1(:,:,3,:));

    CalStru.Shifting = 1;
    auxiliary = (round((CalStru.NFFT/2)*(CalStru.afrq-0.01)):round((CalStru.NFFT/2)*(CalStru.afrq+0.01)))- round((CalStru.Shiftw/(pi*2))*CalStru.NFFT);
    auxiliary(auxiliary<1)=[];
    
    
    JS1o = fringe2image(fringeA,CalStru);
    JS2o = fringe2image(fringeB,CalStru);
%     JS1o(auxiliary,:) = 0;
%     JS2o(auxiliary,:) = 0;
%     
    %S1b = tom2Stokes(JS1,JS2);

    JS1 = JS1o-J1o;
    JS2 = JS2o-J2o;

    IS = abs(JS1o)+abs(JS2o);
    ISsnr = 20*log10(IS)-mean(mean(20*log10(IS(100:150,:))));
    ISsnr(auxiliary,:) = 0;
    %ISsnr(ISsnr<10) = 0;
   
    

    S1 = tom2Stokes(J1o,J2o);
    S2 = tom2Stokes(JS1,JS2);
    
    if isfield(CalStru,'GPU')
        if CalStru.GPU
            S1 = gpuArray(single(S1));
            S2 = gpuArray(single(S2));
        else
        end
    else
    end

    SNR_mask = repmat(ISsnr<20,1,1,3);
    
    S2sel = S2;
    S1sel = S1;
    S1sel(SNR_mask) = NaN;
    S2sel(SNR_mask) = NaN;
    S1sel(1:100,:,:) = NaN;
    S2sel(1:100,:,:) = NaN;
    S1sel(end-100:end,:,:) = NaN;
    S2sel(end-100:end,:,:) = NaN;    
    
    I1 = mean(sqrt(dot(S1sel,S1sel,3)),'omitnan');
    I2 = mean(sqrt(dot(S2sel,S2sel,3)),'omitnan');


    Acorr = sqrt(mean(gather(I2./I1),'omitnan'));
    
%     JS1 = JS1o-J1o.*Acorr;
%     JS2 = JS2o-J2o.*Acorr;
%     S2 = tom2Stokes(JS1,JS2);
%     S2sel = S2;
%     S2sel(SNR_mask) = NaN;
%     S2sel(1:100,:,:) = NaN;
%     S2sel(end-100:end,:,:) = NaN;    
    ms2p = squeeze(mean(S2sel,1,'omitnan'));
    ms1p = squeeze(mean(S1sel,1,'omitnan'));
    
    %Acorr = reshape(Acorr,[1,1,5]);

    n_ms2p = ms2p./repmat(max(sqrt(dot(ms2p,ms2p,2)),1e-9),[1,3]);
    n_ms1p = ms1p./repmat(max(sqrt(dot(ms1p,ms1p,2)),1e-9),[1,3]);

    tp = cross(n_ms2p,n_ms1p,2);
    Qarray = repmat([0,1,0],length(tp),1);
    tQ = cross(Qarray,n_ms1p,2);

    n_tQ = tQ./repmat(max(sqrt(dot(tQ,tQ,2)),1e-9),[1,3]);
    n_tp = tp./repmat(max(sqrt(dot(tp,tp,2)),1e-9),[1,3]);

    thit = gather(unwrap(atan2(sum(cross(n_tQ,n_tp,2).*n_ms1p,2),...
        (sum(n_tQ.*n_tp,2)-sum(n_ms1p.*n_tQ,2).^2))));
    thit(isnan(thit)) = 0;
    thit = (thit)';
    
end

