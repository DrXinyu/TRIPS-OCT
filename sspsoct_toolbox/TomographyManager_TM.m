classdef TomographyManager_TM < handle

    properties
        CalStru
        fringeA
        fringeB
        Amp_factor =1
        I
        pn

        J1
        J2
        JS1
        JS2
        snr_map
        
        frame_index
        shifts
        last_frame
        
        x_filter_width = 6
        z_filter_width = 6
        
        
        
        DepthROI = [];
        
    end
    
    methods
        
        function self = TomographyManager_TM(CalStru)
            self.CalStru = CalStru;
            self.CalStru.whole_depth = 1.5493e+04;
            %self.CalStru.NFFT = 4096;        
            self.frame_index = 0;
            self.CalStru.MinDepth = 2;
            
        end
        
        function setDepthROI(self,roi)
            self.DepthROI = roi;
        end

        function load_fringes(self,fringeA_o,fringeB_o)
            
            [self.fringeA,self.fringeB,~] = k_stabilize(fringeA_o,fringeB_o,self.CalStru,0);

            self.CalStru.Complex = 1;
            self.CalStru.Binning = 0;
            self.CalStru.Shifting = 0;
            
            self.J1 = fringe2image(self.fringeA,self.CalStru);
            self.J2 = fringe2image(self.fringeB,self.CalStru)*1.4;%1.25
            self.I = sqrt(abs(self.J1).^2+abs(self.J2).^2);   
            
            [self.pn,~] =size(self.fringeA);
              
        end  
        
        function load_simulation_fringes(self,fringeA_o,fringeB_o)
            self.fringeA = fringeA_o;
            self.fringeB = fringeB_o;
            

            %self.CalStru.MinDepth = 1;
            self.CalStru.Complex = 1;
            self.CalStru.Binning = 0;
            self.CalStru.Shifting = 0;
            
            self.J1 = fringe2image(self.fringeA,self.CalStru);
            self.J2 = fringe2image(self.fringeB,self.CalStru);
            self.I = sqrt(abs(self.J1).^2+abs(self.J2).^2);   

        end 
    
        
        function I = intensity(self)   
            I = self.I;
            if ~isempty(self.DepthROI)
                I = I(self.DepthROI(1):self.DepthROI(2),:);
            end
        end       

        
        function S1 = get_Stokes(self,cfactor)
            if nargin<2
                cfactor = 1.25;
            end
            S1 = tom2Stokes(self.J1,self.J2*cfactor);      
            
            
            if ~isempty(self.DepthROI)
                S1 = S1(self.DepthROI(1):self.DepthROI(2),:,:);
                
            end

        end
        
        function S = get_Stokes_binning(self,binning,cfactor)

            if nargin<3
                cfactor = 1.25;
            end
            
            
            %self.CalStru.MinDepth = 1;
            self.CalStru.Complex = 1;
            self.CalStru.Binning = binning;
            
            self.CalStru.Shifting = 0;
            
            J1 = fringe2image(self.fringeA,self.CalStru);
            J2 = fringe2image(self.fringeB,self.CalStru)*cfactor; %1.25
            
            S = tom2Stokes(J1,J2);      
            
            if ~isempty(self.DepthROI)
                S = S(self.DepthROI(1):self.DepthROI(2),:,:,:);
            end            

        end        
        
        
        
        function [S,I] = filter_by_bin(self,S)
            
            fract = self.CalStru.Binning;
            wnum = 2*fract - 1;
            N = self.pn;
            W = round(N/fract);
            window = cat(1,hanning(W),zeros(N-W,1));
            filt = (abs(fft(window,self.CalStru.NFFT)));
            filt = cat(1,flipud(filt(2:51)),filt(1:50));
            filt = filt.*(1./(sum(filt)));
            
            I = sqrt(dot(S,S,4));
            S = imfilter(S,filt,'same');
            I = imfilter(I,filt,'same');

        end
        
        
         function [S,I] = lateral_filtering(self,S,I,fwx)
             
            filter = (gausswin(fwx*2,2.3)./sum(gausswin(fwx*2,2.3)))';    
            
            S = imfilter(S,filter,'same');
            I = imfilter(I,filter,'same');

        end
               
        
        
        
        
        
        
        
        function fk = get_bin_kernal(self,binning)
            self.CalStru.Binning = binning;
            fract = self.CalStru.Binning;
            wnum = 2*fract - 1;
            N = self.pn;
            W = round(N/fract);
            window = cat(1,hanning(W),zeros(N-W,1));
            filt = (abs(fft(window,self.CalStru.NFFT)));
            filt = cat(1,flipud(filt(2:51)),filt(1:50));
            fk = filt.*(1./(sum(filt)));
        end
            
        
        
        
        function JM = get_Jones_binning(self,binning)
            
            %self.CalStru.MinDepth = 1;
            self.CalStru.Complex = 1;
            self.CalStru.Binning = binning;
            
            self.CalStru.Shifting = 0;
            
            J1 = fringe2image(self.fringeA,self.CalStru);
            J2 = fringe2image(self.fringeB,self.CalStru);

            JM = cat(4,J1,J2);
            
            if ~isempty(self.DepthROI)
                JM = JM(self.DepthROI(1):self.DepthROI(2),:,:,:);
            end   
        end           
        
        
        
        
        function [snr_mask,I1,I2] = SNRmask(self,I1,I2,lowlevel,highlevel)
            
            
            SNR = self.SNRmap();
            snr_mask1 = (20.*log10(SNR)) > lowlevel;
            snr_mask2 = (20.*log10(SNR))< highlevel;

            snr_mask = snr_mask1 & snr_mask2;
            
        end
       
        
        function snr_map = SNRmap(self,I)
                   
            ground_noise = quantile(I(:),0.1);
            snr_map = I./ground_noise;
   
        end
        

        
        function nim = avg_imresize(self,im,newW)
            dim = size(im);
            kernal_W = round(dim(2)./newW);
            if kernal_W >= 1
                k = ones(1,kernal_W);
                im = imfilter(im,k);
            end
            if kernal_W == dim(2)/newW
                nim = im(:,1:kernal_W:end,:,:);
            else
                im = permute(im,[2 1 3 4]);
                im = reshape(im,dim(2),[]);
                nx = linspace(1,dim(2),newW);
                im = interp1((1:dim(2))',im,nx);
                im = reshape(im,newW,dim(1),[]);
                im = permute(im,[2 1 3 4]);
                dim(2) = newW;
                nim = reshape(im,dim);

            end
        end
        
        
        
        


        
        function single_frame_stabilize_setup(self,signature_frame,this_frame,frameindex,alineROI)
            self.frame_index = self.frame_index+1;
            if frameindex ~= self.frame_index
                disp("frame stabilize error: not the right index");
                return
            end
            
            if self.frame_index == 1
                self.shifts(self.frame_index) = 0;
            else
                signature_frame = imgaussfilt(signature_frame,20);
   
                if isempty(alineROI)      
                    [value,shiftframe]= max(self.frame_correlation(this_frame,signature_frame),[],2);
                else
                    [value,shiftframe]= max(self.frame_correlation(this_frame(:,alineROI(1):alineROI(2)),signature_frame(:,alineROI(1):alineROI(2))),[],2);
                end
                level = quantile(value,0.75);
                shiftframe = gather(shiftframe);
                self.shifts(self.frame_index) = (mean(shiftframe(value>level)))-21;
            end            

        end
        

        function stabilize_setup(self,signature_frame,frameindex,alineROI)
            
            self.frame_index = self.frame_index+1;
            if frameindex ~= self.frame_index
                disp("frame stabilize error: not the right index");
                return
            end
            
            if self.frame_index == 1
                self.shifts(self.frame_index) = 0;
            else
                signature_frame = imgaussfilt(signature_frame,3);
   
                if isempty(alineROI)      
                    [value,shiftframe]= max(self.frame_correlation(signature_frame,self.last_frame),[],2);
                else
                    [value,shiftframe]= max(self.frame_correlation(signature_frame(:,alineROI(1):alineROI(2)),self.last_frame(:,alineROI(1):alineROI(2))),[],2);
                end
                level = quantile(value,0.75);
                shiftframe = gather(shiftframe);
                self.shifts(self.frame_index) = (mean(shiftframe(value>level)))-21;
            end            
            self.last_frame = signature_frame;
            
        end
        
        
        
        function xc = frame_correlation(self,frame1,frame2)
            frame1 = frame1 - mean(frame1);
            frame2 = frame2 - mean(frame2);
            frame1array = repmat(frame1,1,1,41);
            for findex = 1:41
                frame1array(:,:,findex) = circshift(frame1array(:,:,findex),findex-21,1);
            end
            xc = squeeze(sum(abs(frame1array+frame2),1));
        end 
        
        function frames = stabilize_shift(self,frames,smoothsize,frame_index)
            if isempty(smoothsize)
                motion = self.CalStru.FrameShift' - mean(self.CalStru.FrameShift);            
            else
                motion = self.CalStru.FrameShift' - smooth(self.CalStru.FrameShift,smoothsize);   
            end

            frames = circshift(frames,round(motion(frame_index)),1);  
        end
        
        
        function S = Stokes_medfilt2(self,S,kernal)
            if nargin<3
                kernal = [3 3];
            end
            [H,W,A,B] = size(S);
            
            S = reshape(S,H,[]);
            S = medfilt2(S,kernal);
            
            S = reshape(S,H,W,A,B);
            
        
        end
        
        
        
    
    end
end

