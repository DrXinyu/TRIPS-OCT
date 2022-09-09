function [shiftw,CalStru] = alignShift(fringeA,fringeB,CalStru)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    old_CalStru= CalStru;   
    
    CalStru.NFFT = 4096*16;
    NFFT = CalStru.NFFT;
    if CalStru.Shiftw/(2*pi) > 0.08
        CalStru.Shiftw = 0.08*(2*pi);
    else
        
    end
    
    
    encode_depth_n = round((CalStru.Shiftw/(pi*2))*NFFT);
    wscale = linspace(0,2*pi,NFFT);
    CalStru.K_stablize = 0;
    
    CalStru.Complex = 0;
    CalStru.Binning = 0;
    CalStru.Background = 'mean';
    CalStru.Shifting = 0;
    CalStru.MinDepth = 0;
    CalStru.Window = 1;
    CalStru.K_stabilized = 1;
    
    Alinex = fringe2image(fringeA,CalStru);
    Aliney = fringe2image(fringeB,CalStru);
    
    image = abs(Alinex)+abs(Aliney);
    
    image = 20*log10(image) - mean(mean(20*log10(image(200:300,:))));
    image = image(1:encode_depth_n*2,:);
    
    shallow_image = image;
    shallow_image(encode_depth_n-1000:end,:) = 0;
    [SAmp,shallow_position] = max(shallow_image);
    
    SAmp_mask = SAmp>5;

    deep_image = image;
    deep_image([1:encode_depth_n+1000,encode_depth_n*2-1000:end],:) = 0;
    [DAmp,deep_position] = max(deep_image);
    
    norm_shallow_image = (shallow_image-mean(shallow_image(:)))./var(shallow_image(:));
    norm_deep_image = (deep_image-mean(deep_image(:)))./var(deep_image(:));
    
    for xindex = 4000:encode_depth_n+1000
         s_i1 = circshift(norm_shallow_image,xindex,1);
         xce(xindex) = sum(sum((s_i1+norm_deep_image).^2));
    end
    
    [~,mn] = max(xce);
    shiftw = (mn/(NFFT-1))*2*pi;

    CalStru = old_CalStru; 
    CalStru.Shiftw = shiftw;
    CalStru.Deepamp = 1;
end

