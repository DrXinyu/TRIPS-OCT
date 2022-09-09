classdef Masks_and_Filters < handle
   
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function self = Masks_and_Filters()
        end
        
        function map = snr_mask(self,intensity,level)
            ground_noise = quantile(intensity,0.1);
            SNR = intensity./ground_noise;
            SNRlog = 20*log10(SNR);
            map = (SNRlog>level(1))&(SNRlog<level(2));
        end
        
        function snr = int2snr(self,intensity)
            ground_noise = quantile(intensity,0.1);
            SNR = intensity./ground_noise;
            snr = 20*log10(SNR);
        end        
        
        
        
        
        
        function map = orthogonality_mask(self)
            
        end
        
        function map = dop_mask(self,S,fw,level)
            
            I = sqrt(dot(S,S,3));
            kernel = (gausswin(fw,1)./sum(gausswin(fw,1)))';

            Sf = imfilter(S,kernel,'replicate');
            If = imfilter(I,kernel,'replicate');

            kernel = gausswin(fw,1)./sum(gausswin(fw,1));

            Sf = imfilter(Sf,kernel,'replicate');
            If = imfilter(If,kernel,'replicate');

            % normalization
            QUVf = sqrt(dot(Sf,Sf,3));
            dop = QUVf./If;

            map = dop>level;
        end        
        
        
        
        
    end
end

