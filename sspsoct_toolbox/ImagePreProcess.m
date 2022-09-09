classdef ImagePreProcess  < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        GPU
        NFFT
    end
    
    methods
        function self = ImagePreProcess(CalStru)
            self.GPU = 0;
            if isfield(CalStru,'GPU')
                if CalStru.GPU
                    self.GPU = 1;
                else
                end
            else
            end
            
            self.NFFT = CalStru.NFFT;
        end

        function flatten_curve = setup_flattening(self,IBscan)
            
            [D,L] = size(IBscan); 

            if self.GPU
                IBscan = gpuArray(single(IBscan));
            end


            If = imgaussfilt(IBscan,8);
            If = 20*log10(self.SNRmap(If));
            
            BW = If>4;            
            BW(1:round(self.NFFT.*0.0073),:) = 0;
            BW(round(self.NFFT.*0.07):end,:) = 0;
            

            se = strel('disk',round(self.NFFT.*0.0015));
            BW = imopen(BW,se);
            se = strel('disk',round(self.NFFT.*0.0015));
            BW = imclose(BW,se);
            se = strel('line',round(self.NFFT.*0.01),90);
            BW = imerode(BW,se);
            If = If.*BW;

            [ss,peak_depth] = max(If);
            peak_depth = gather(peak_depth);
            threshold = gather(quantile(ss,0.1));
            fitmask = gather(ss>threshold);
            if sum(fitmask(:)) < 10
                fitmask = ss>min(ss(:));
            end
            fitx = 1:L;
            flatten_curve = gather(polyfit(fitx(fitmask),peak_depth(fitmask),6));

        end
        function showit(self,image,flatten_curve)
            image = gather(image);
            IS = size(image);
            D = IS(1);            
            L = IS(2);
            flaten_s = round(polyval(flatten_curve,1:L));
            
            flaten_s(flaten_s<1) = 10;
            flaten_s(flaten_s>(D-10)) = D-10;
            fp = sub2ind(IS,flaten_s,1:L);
            image(fp) = max(image(:));
            imagesc(image);
 
        end
        
        
        
        
        
        function new_im_seq = flatten_image(self,image,flatten_curve)
        %flatten image
            image = gather(image);
            IS = size(image);
            L = IS(2);
            D = IS(1);

            flaten_s = round(polyval(flatten_curve,1:L)+D);
            flaten_s(flaten_s<(D+1)) = 10+D;
            flaten_s(flaten_s>(2*D-10)) = 2*D-10;

            flaten_mask_row = repmat(flaten_s,2*D,1);
            flaten_mask_row = flaten_mask_row - repmat((D:-1:-D+1)',1,L);
            flaten_mask_col = repmat((1:L),2*D,1);

            im_seq = reshape(image,D,L,[]);
            [~,~,n] = size(im_seq);
            new_im_seq = ones(D*2,L,n);

            for index = 1:n
                newimage = ones(D*3,L).*1e-9;
                newimage(D+1:2*D,:) = im_seq(:,:,index);
                indxp = sub2ind(size(newimage),flaten_mask_row,flaten_mask_col);
                aaa = newimage(indxp);
                new_im_seq(:,:,index) = aaa;  

            end
            IS(1) = 2*D;
            new_im_seq = reshape(new_im_seq,IS);


        end
        
        
        function SNR = SNRmap(self,I)
            
            ground_noise = quantile(I(:),0.1);
            SNR = I./ground_noise;
  
        end
        
        
        
  
        
    end
end

