function image = fringe2image(fringe,CalStru)
    
    if isfield(CalStru,'GPU')
        if CalStru.GPU == 1
            fringe = gpuArray(single(fringe));
        else
        end
    else
    end
    
    if isfield(CalStru,'Background')
        if strcmp(CalStru.Background,'smooth')
            fringe = fringe-smooth(mean(fringe,2),20);
        elseif strcmp(CalStru.Background,'mean')
            fringe_amp = max(fringe)-min(fringe);
            fringe = fringe-mean(fringe(:,fringe_amp<quantile(fringe_amp,0.8)),2);
        elseif strcmp(CalStru.Background,'none')            
        end
    else
    end
    
    if isfield(CalStru,'K_stabilized')
        if CalStru.K_stabilized == 1
            yC = fringe;
        else
            yC = k_sample_fringe(fringe,CalStru);
        end
    else
        yC = k_sample_fringe(fringe,CalStru);
    end  
    [pn,L] = size(yC);
    

    %plot(yC)
    
    
    % windowing    
    
    if isfield(CalStru,'Window') 
        if isfield(CalStru,'Binning')
            if CalStru.Binning ~= 0
                CalStru.Window = 0;
            end
        end
        if CalStru.Window == 1
             win = hanning(pn);
             win = win/(sqrt(sum(win.^2))/sqrt(pn));
             yC = bsxfun(@times,yC,win);
        else
        end
    else
    end
    
    % dispersion correction
    
    if isfield(CalStru,'Dispersion')
        if CalStru.Dispersion
            yCdp = bsxfun(@times,yC,CalStru.CArray);
        else
            yCdp = yC;
        end
    else
        yCdp = bsxfun(@times,yC,CalStru.CArray);
    end
    
    % shifting
    
    if isfield(CalStru,'Shifting')
        
        if CalStru.Shifting
            yCdp = yCdp.*exp(-1i*CalStru.Shiftw.*(0:(length(CalStru.xi)-1))');    
        end
        
    else  
    end      
    
    
    % binning
    
    if isfield(CalStru,'Binning')
        
        if CalStru.Binning > 1
            fract = CalStru.Binning;
            wnum = 2*fract - 1;
            N = pn;
            W = round(N/fract);
            scale = sqrt(sum(hanning(N).^2)/sum(hanning(W).^2));
            for wind = 1:wnum
                window(:,wind) = circshift(cat(1,hanning(W),zeros(N-W,1)),round((wind-1)*N/(wnum+1)))*scale;
            end
            window = permute(shiftdim(window,-1),[2,1,3]);
            yCdp = bsxfun(@times,yCdp,window); 
        end
        
    else  
    end
    
  
    
    
    
    % FFT
    
    if isfield(CalStru,'NFFT')
        NFFT = CalStru.NFFT;
    else
        NFFT = 4096*4;
    end

    
    
    % complex image
    if isfield(CalStru,'Complex')
        if CalStru.Complex
            imagef = fft(yCdp,NFFT);           
        else
            imagef = abs(fft(yCdp,NFFT));
        end
    else
        imagef = abs(fft(yCdp,NFFT));
    end

    
    
    % full range

    if isfield(CalStru,'Fullrange')
        if CalStru.Fullrange
            image = (imagef);
        else
            image = ((imagef(1:NFFT/2,:)));
        end
    else
        image = ((imagef(1:NFFT/2,:)));
    end
    
    
    
    if isfield(CalStru,'MinDepth')
        if CalStru.MinDepth > 0
            encode_depth_n = round(round((CalStru.Shiftw/(2*pi))*NFFT).*CalStru.MinDepth);
            image = image(1:encode_depth_n,:);
        else
        end
    else
    end    
    % binning output
    
    if isfield(CalStru,'Binning')
        
        if CalStru.Binning > 1
            [D,W] = size(image);
            image = reshape(image,D,[],wnum);
        end
        
    else  
    end
    
    
    
    
end

function yC = k_sample_fringe(fringe,CalStru)
    [spec,line] = size(fringe);
    
    if spec == length(CalStru.xi)
        yC = fringe;    
    elseif spec ~= length(CalStru.MAmean)
        fringe = fringe(CalStru.Ang_Cutting_Limit(1):CalStru.Ang_Cutting_Limit(2),:);
        yC = interp1(CalStru.MAmean,fringe,CalStru.xi);   
%         if isfield(CalStru,'MAmean2')
%             yC = interp1(CalStru.MAmean2,yC,CalStru.xi2); 
%         end
        
    else
        yC = interp1(CalStru.MAmean,fringe,CalStru.xi);   
%         if isfield(CalStru,'MAmean2')
%             yC = interp1(CalStru.MAmean2,yC,CalStru.xi2); 
%         end        
    end
   
    if isfield(CalStru,'Search')
        if CalStru.Search == 1
            yC = yC./(sqrt(mean((yC).^2)));
        else
        end
    else
    end
    


end

