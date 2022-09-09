function flatten_curve = setup_flattening(IBscan,CalStru)
% locate the max line in the image
    [D,L] = size(IBscan); 
    if isfield(CalStru,'GPU')
        if CalStru.GPU
            IBscan = gpuArray(single(IBscan));
        else
        end
    else
    end
    
    If = imgaussfilt(IBscan,8);
    If = 20*log10(If)-mean(20*log10(If(150,:)));
    BW = If>2;
    se = strel('disk',80);
    BW = imopen(BW,se);
    se = strel('disk',80);
    BW = imclose(BW,se);
    se = strel('line',180,90);
    BW = imerode(BW,se);
    If = If.*BW;

    [ss,peak_depth] = max(If(350:end-150,:));
    peak_depth = gather(peak_depth+350);
    threshold = gather(quantile(ss,0.1));
    fitmask = gather(ss>threshold);
    if sum(fitmask(:)) < 10
        fitmask = ss>min(ss(:));
    end
    fitx = 1:L;
    flatten_curve = gather(polyfit(fitx(fitmask),peak_depth(fitmask),6));

end

