classdef Dop_contrast < handle

    properties
        ground_noise
        dop_SNR_curve
        
    end
    
    methods
        function self = Dop_contrast()
            if ~exist('DOPCache', 'dir')
               mkdir('DOPCache')
            end
        end
        
        function load_calibration_frame(self,dop,intensity,fileID)
            cali_frame.dop = dop;
            cali_frame.intensity = intensity;
            save(fullfile('DOPCache',strcat(fileID,'.mat')),"cali_frame"); 
        end
        
        function calibrate(self)
            
            filenameN = dir('DOPCache');
            dop = [];
            intensity = [];
            for file_index = 1:length(filenameN)
                if filenameN(file_index).isdir
                    continue
                end
                
                data_file = fullfile(filenameN(file_index).folder,filenameN(file_index).name);
                frame = importdata(data_file);
                dop = cat(1,frame.dop(:),dop);
                intensity = cat(1,frame.intensity(:),intensity);
            end


            
            SNRdata = 20*log10(self.SNRmap(intensity));
            large_end = max(SNRdata);
            figure
            plot(SNRdata,dop)
            hold on
            binsetp = 1;
            for SNRbindex = 1:binsetp:large_end
                SNRbin_low = SNRbindex-binsetp;
                SNRbin_high = SNRbindex;
                select = (SNRdata>=SNRbin_low)&(SNRdata<SNRbin_high);
                Edop(SNRbindex) = gather(mean(dop(select),'omitnan'));
            end
             plot(Edop,'r')
            Edop(isnan(Edop)) = Edop(circshift(isnan(Edop),-1));
            self.dop_SNR_curve = fit(gather((1:binsetp:large_end)'),Edop','smoothingspline');   
            
        end
        
        
        function [ddop,nanddop] = get_contrast(self,dop,intenisty)
            if length(size(intenisty)) ~= length(size(dop))
                intenisty = mean(intenisty,3);
            end
            
            SNRm = 20*log10(self.SNRmap(intenisty));
            %SNRm = 20*log10(intenisty./self.ground_noise);
            tabledop = reshape(feval(self.dop_SNR_curve,SNRm),size(intenisty));
            ddop = (tabledop-dop);
            nanddop = ddop;
            mask = (SNRm)<0;
            ddop(mask) = 0;
            nanddop(mask) = nan;
            %ddop(ddop<0) = 0; 
        end
        
        
        function dop = averagedStokes2Dop(self,avgS,avgI)
            
            dim = size(avgS);
            SI = dot(avgS,avgS,length(dim));
            dop = SI./avgI.^2;
            
            dop = mean(dop,length(dim)-1);

            
            
            
        end
        
        
        
        
        function SNR = SNRmap(self,I)
            
            %I = reshape(I,self.dim(1)*self.dim(2),[]);
            self.ground_noise = quantile(I(:),0.1);
            SNR = I./self.ground_noise;
            %SNR = reshape(SNR,self.dim(1:end-1));
  
        end
        
        
        
        
        
        
    end
end

