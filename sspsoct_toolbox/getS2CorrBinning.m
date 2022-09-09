function [Acorr,thit] = getS2CorrBinning(fringeA,fringeB,CalStru)
%S1S2
%
    %CalStru.NFFT = 4096*2;
    CalStru.Complex = 1;
%    CalStru.Binning = 0;
    CalStru.Background = 'mean';
    CalStru.MinDepth = 1;
    CalStru.GPU = 1;
    w = 2*CalStru.Binning-1;
    CalStru.Shifting = 0;
    J1o = fringe2image(fringeA,CalStru);
    J2o = fringe2image(fringeB,CalStru);


    %S1 = squeeze(S1(:,:,3,:));

    CalStru.Shifting = 1;
    JS1o = fringe2image(fringeA,CalStru);
    JS2o = fringe2image(fringeB,CalStru);
    %S1b = tom2Stokes(JS1,JS2);



    JS1 = JS1o-J1o;
    JS2 = JS2o-J2o;

    IS = abs(JS1o)+abs(JS2o);
    ISsnr = 20*log10(IS)-mean(mean(mean(20*log10(IS(50:75,:)))));
    ISsnr(ISsnr<9) = 0;


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
    


    SNR_mask = repmat(ISsnr<9,1,1,1,3);
    S2sel = S2;
    S1sel = S1;
    S1sel(SNR_mask) = NaN;
    S2sel(SNR_mask) = NaN;
    ms2p = squeeze(mean(S2sel,1,'omitnan'));
    ms1p = squeeze(mean(S1sel,1,'omitnan'));
    E1 = sqrt(dot(ms1p,ms1p,3));
    E2 = sqrt(dot(ms2p,ms2p,3));
    Acorr = gather(median(E2./E1,'omitnan'));
    Acorr = reshape(Acorr,[1,1,w]);
    n_ms2p = ms2p./repmat(max(sqrt(dot(ms2p,ms2p,3)),1e-9),[1,1,3]);
    n_ms1p = ms1p./repmat(max(sqrt(dot(ms1p,ms1p,3)),1e-9),[1,1,3]);


    tp = cross(n_ms2p,n_ms1p,3);
    Q(1,1,1) = 0; Q(1,1,2) = 1; Q(1,1,3) = 0;
    Qarray = repmat(Q,size(tp(:,:,1)));
    tQ = cross(Qarray,n_ms1p,3);

    n_tQ = tQ./repmat(max(sqrt(dot(tQ,tQ,3)),1e-9),[1,1,3]);
    n_tp = tp./repmat(max(sqrt(dot(tp,tp,3)),1e-9),[1,1,3]);
    
    n_ms1p = gather(n_ms1p);
    n_tQ = gather(n_tQ);
    n_tp = gather(n_tp);
    thit = gather(unwrap(atan2(sum(cross(n_tQ,n_tp,3).*n_ms1p,3),(sum(n_tQ.*n_tp,3)-sum(n_ms1p.*n_tQ,3).^2))));


    thit(isnan(thit)) = 0;
    thit = reshape(thit,size(IS(1,:,:)));

end

