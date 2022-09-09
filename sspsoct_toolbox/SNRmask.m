function [snr_mask,I1,I2] = SNRmask(I1,I2,lowlevel,highlevel)
    dim = size(I1);
    kernel = (gausswin(15,1)./sum(gausswin(15,1)))';
    I1 = imfilter(I1,kernel,'replicate');
    I2 = imfilter(I2,kernel,'replicate');
    kernel = kernel';
    I1 = imfilter(I1,kernel,'replicate');
    I2 = imfilter(I2,kernel,'replicate');

    I = (I1+I2)/2;            
    I = reshape(I,dim(1)*dim(2),[]);

    ground_noise = quantile(I,0.1);
    SNR = I./ground_noise;
    SNR = reshape(SNR,dim);
    snr_mask1 = (20.*log10(SNR)) > lowlevel;
    snr_mask2 = (20.*log10(SNR))< highlevel;

    snr_mask = snr_mask1 & snr_mask2;
end

