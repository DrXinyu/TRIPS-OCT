classdef EnfaceStokesMap < handle

    
    properties
       maf
       enface_Stokes
       enface_intensity
       ONH_location
    end
    
    methods
        
        function self = EnfaceStokesMap()
            self.maf = Masks_and_Filters();
            self.ONH_location = [0 0];
        end
        
        function add_frame(self,Stokes_frame,intensity,num)
            snr = self.maf.snr_mask(intensity,[2 30]);
            dop = self.maf.dop_mask(Stokes_frame,10,0.6);
            
            nanmask = double(snr&dop);
            nanmask((snr&dop)==0) = nan;
            
            Stokes_frame = nanmask.*Stokes_frame;
            intensity_frame = nanmask.*intensity;
            
            line_intensity = mean(intensity_frame,1,'omitnan');
            line_Stokes = median(Stokes_frame,1,'omitnan');
            
            self.enface_Stokes(num,:,:) = gather(line_Stokes);
            self.enface_intensity(num,:) = gather(line_intensity);
            
        end
        
        
        function setONH(self,xy)
            self.ONH_location = xy;
        end
        
        
        function [S1n, S2n, S3n] = tripod_Stokes(self,S1,S2)
            
            S1 = S1./repmat(max(sqrt(dot(S1,S1,3)),1e-9),[1 1 3]);
            S2 = S2./repmat(max(sqrt(dot(S2,S2,3)),1e-9),[1 1 3]);
            
            na = S1 + S2;
            nb = S1 - S2;
            na = na./repmat(max(sqrt(dot(na,na,3)),1e-9),[1 1 3]);
            nb = nb./repmat(max(sqrt(dot(nb,nb,3)),1e-9),[1 1 3]);
%             S1n = (na + nb)/sqrt(2);
%             S2n = (na - nb)/sqrt(2);

            S1n = na;
            S2n = nb;
            
            S3n = cross(S1n,S2n,3);

        end
        
        
        function [ret, PA] = dual_input_mapping(self,em2)
            
            S1plus = self.enface_Stokes;
            S1plus(isnan(S1plus)) = 1e-20;
            S1plus = imgaussfilt(S1plus,2);
            
            S1minus = S1plus(self.ONH_location(1),self.ONH_location(2),:).*ones(size(S1plus));
            
            S2plus = em2.enface_Stokes;
            S2plus(isnan(S2plus)) = 1e-20;
            S2plus = imgaussfilt(S2plus,2);
            
            S2minus = S2plus(em2.ONH_location(1),em2.ONH_location(2),:).*ones(size(S2plus));            
            
            [S1n_minus, S2n_minus, S3n_minus] = self.tripod_Stokes(S1minus,S2minus);
            [S1n_plus, S2n_plus, S3n_plus] = self.tripod_Stokes(S1plus,S2plus);
            
            PA1 = cross(S1n_minus-S1n_plus,S2n_minus-S2n_plus,3);
            PA2 = cross(S2n_minus-S2n_plus,S3n_minus-S3n_plus,3); 
            PA3 = cross(S1n_minus-S1n_plus,S3n_minus-S3n_plus,3); 
            
            [~,N]=max(abs(PA1)+abs(PA2)+abs(PA3),[],3);
            
            PA1 = sign(PA1(N)+1e-20).*PA1;
            PA2 = sign(PA2(N)+1e-20).*PA2;
            PA3 = sign(PA3(N)+1e-20).*PA3;
            PA = PA1+PA2+PA3;
            PA  = PA./repmat(max(sqrt(dot(PA,PA,3)),1e-9),[1 1 3]);
            
            norm_pm = sign(dot(cross(S1n_plus,S1n_minus,3),PA,3));            
            PA = PA.*norm_pm;  
            temp = dot(S1n_plus,PA,3).^2;

%            pm = sign((1-dot(S1plus,S1minus)).*(dot(S1plus-S1minus,S2plus+S2minus)));
            ret= real(acos(complex(dot(S1n_plus,S1n_minus,3)-temp)./(1-temp)));
         
        end
        
        
        function orth = orthogonolity(self,em2)
            
            S1= self.enface_Stokes;
            S1 = imgaussfilt(S1,2);
            S1  = S1./repmat(max(sqrt(dot(S1,S1,3)),1e-9),[1 1 3]);
            S2 = em2.enface_Stokes;
            S2 = imgaussfilt(S2,2);
            S2  = S2./repmat(max(sqrt(dot(S2,S2,3)),1e-9),[1 1 3]);
            orth = acosd(gather(dot(S1,S2,3)));
            
        end
        
        
        
        
    end
end

