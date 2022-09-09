classdef FlattenManager < handle

    properties
        flat_surface
        line = ones(1,40)/40;
        flat_surface_array = [];
        RPE = [];
        sclera=[];
        ex_vivo_surface = [];
        
    end
    
    methods
        function self = FlattenManager()
        end
        
        function PSdual_Manager = flattenMulluerMatrix(self,mask,PSdual_Manager)
            MMc = permute(PSdual_Manager.MMc,[2,3,1]);
            [~,surface] = max(mask>0);
            surface = round(smooth(surface,100))';
            
            [H, W] = size(mask);
            D = max(surface);
            surface_sub_X = repmat(surface,H,1)+repmat((1:H)',1,W);
            surface_sub_Y = repmat(1:W,H,1);            
            indxp = sub2ind([H+D,W],surface_sub_X,surface_sub_Y);
            
            [H, W, B] = size(MMc);

            for index =1:B
                MMcex = cat(1,MMc(:,:,index),zeros(D,W));
                MMc(:,:,index) = MMcex(indxp);
            end
            PSdual_Manager.MMc = permute(MMc,[3,1,2]);
        
        end
        
        
        function [I,S] = flatten2edge(self,mask,I,S)

            [~,surface] = max(mask>0);
            surface = round(smooth(surface,100))';
            surface = surface+5;
            
            [H, W] = size(mask);
            D = max(surface);
            surface_sub_X = repmat(surface,H,1)+repmat((1:H)',1,W);
            surface_sub_Y = repmat(1:W,H,1);            
            indxp = sub2ind([H+D,W],surface_sub_X,surface_sub_Y);
            if ~isempty(I)
                
                Iex = cat(1,I,zeros(D,W));
                I = Iex(indxp);
                
            end

            if ~isempty(S)
            
                [H, W, B, St] = size(S);
                Sres = reshape(S,H,W,[]);

                for index =1:B*St
                    Sex = cat(1,Sres(:,:,index),zeros(D,W));
                    Sres(:,:,index) = Sex(indxp);
                end

                S = reshape(Sres,H,W,B,St); 
            end
        end

        function intelligent_flatten_setup_using_RPE(self,dop,snr)
            [H, W] = size(snr);

            mask = (snr);
            mask = diff(mask);
            [~,sleek_surface] = max(mask);
            sleek_surface = gather(sleek_surface);
% %             
%             [p,~,mu] = polyfit(1:length(sleek_surface),sleek_surface,7);
%             sleek_surface = polyval(p,1:length(sleek_surface),[],mu);
%             

            [flat_curve,is_signal] = self.extract_curvature(snr);
            
            edge_sig = diff(is_signal>0.01);
            edge_sig(edge_sig<0) = 0;
            level_sig = cumsum([0 edge_sig]);
            l2 = quantile(level_sig,0.75);
            if (l2 == 1) || (l2 == 0)
                shift = median(flat_curve(is_signal>0)-sleek_surface(is_signal>0),2);
                self.flat_surface = round(flat_curve-shift);
            elseif l2 > 1
            %l2 = quantile(level_sig,0.75);
                self.flat_surface = flat_curve;
                shift = median(flat_curve((is_signal>0)&(level_sig<=1))-sleek_surface((is_signal>0)&(level_sig<=1)),2);
                self.flat_surface(level_sig<=1) = round(flat_curve(level_sig<=1)-shift);
                for Lindex = 2:max(level_sig)  
                    shift = median(flat_curve((is_signal>0)&(level_sig>=Lindex))-sleek_surface((is_signal>0)&(level_sig>=Lindex)),2);
                    self.flat_surface(level_sig>=Lindex) = round(flat_curve(level_sig>=Lindex)-shift);        
                end
            end
                %             surface(abs(surface-self.flat_surface)>40) = self.flat_surface(abs(surface-self.flat_surface)>40);
            %self.flat_surface = (surface+self.flat_surface)/2;
            
%             b_surface = [ones(1,50)*self.flat_surface(1) self.flat_surface ones(1,50)*self.flat_surface(end)];
%             self.flat_surface = round(conv(b_surface,self.line,'same'));
%             self.flat_surface = self.flat_surface(51:end-50);
            %self.flat_surface = round(sleek_surface_fit);
            self.flat_surface(self.flat_surface<1)=1;
            self.flat_surface(isnan(self.flat_surface))=1;
            
            self.flat_surface = self.inter_frame_smooth(3);
            
        end        
        
        function surface_flatten_setup(self,dop,snr)
            [H, W] = size(snr);
            
            snr = imgaussfilt(snr,3);
            mask = (snr>1.5)&(dop>0.3);
            
            [~,rough_surface] = max(mask>0);
            
            rough_surface = round(medfilt1(gather(rough_surface),10));
            L = 5;
            
            rough_surface(rough_surface>(H-L-1))=1;
            
            surface_sub_X = repmat(rough_surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            upper_retina = snr(indxp);
            
            
%             % generate horizontal edge emphasis kernel
%             h = fspecial('sobel');
%             % invert kernel to detect vertical edges
%             h = h';
%             upper_retina = imfilter(upper_retina,h);
%             
            
            [~,sleek_surface] = max(upper_retina);
            sleek_surface = sleek_surface +rough_surface;
            sleek_surface = round(medfilt1(gather(sleek_surface),10));
            
            dataf = (smooth(sleek_surface',50))';
            
            self.flat_surface = dataf;
            self.flat_surface(self.flat_surface<1)=1;
            self.flat_surface(isnan(self.flat_surface))=1;
            
            self.flat_surface = self.inter_frame_smooth(5);
       
        end
        
        
        function find_RPE(self,snr,rough_estimator)

            snr = imgaussfilt(log(snr),3);
            difsnr = diff(snr);
            [H, W] = size(difsnr);            
            surface = self.flat_surface;

            L = 40; 
            rough_estimator = rough_estimator-5;
            surface(surface>(H-L-1))=1;
            
            surface_sub_X = repmat(surface+rough_estimator,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            around_RPE = difsnr(indxp);
            [~,self.RPE] = max(around_RPE);
            
            self.RPE = surface +self.RPE+rough_estimator;
            self.RPE = round(medfilt1(gather(self.RPE),10))+4;
            
        end
        
        function find_sclera(self,snr,ret)
            
            [H, W] = size(snr);      
            snr = imgaussfilt(log(snr),2);
            ret = imgaussfilt(ret.*(snr>1.5),5);
            retmask = ret>0.06;
            
            surface = self.RPE+10;
            L = 40;

            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            belowRPE = ret(indxp);
            [~,self.sclera] = max(belowRPE);
            self.sclera = self.RPE+self.sclera;
            self.sclera = medfilt1(gather(self.sclera),30);
            self.sclera = (smooth(self.sclera,30))';
            self.sclera = round(self.sclera);
            
            L = 20;
            surface = round(self.sclera)-L;
            surface(surface<self.RPE) = self.RPE(surface<self.RPE);

            difsnr = diff(snr);
            [H, W] = size(difsnr);  
            
            surface(surface>(H-L-1))=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,self.sclera] = max(above_sclera);
            
            self.sclera = surface+self.sclera;
            self.sclera = medfilt1(gather(self.sclera),10);
            self.sclera = (smooth(self.sclera,10))';
            self.sclera = round(self.sclera);            
             
        end       

        function find_pigmented_sclera(self,snr,ret)
            
            [H, W] = size(snr);      
            snr = imgaussfilt(log(snr),2);
            ret = imgaussfilt(ret.*(snr>1.5),5);
            retmask = ret>0.06;
            
            surface = self.RPE+10;
            L = 40;

            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            belowRPE = ret(indxp);
            [~,self.sclera] = max(belowRPE);
            self.sclera = self.RPE+self.sclera;
            self.sclera = medfilt1(gather(self.sclera),30);
            self.sclera = (smooth(self.sclera,30))';
            self.sclera = round(self.sclera);
            
            L = 20;
            surface = round(self.sclera)-L;
            surface(surface<self.RPE) = self.RPE(surface<self.RPE);

            difsnr = diff(snr);
            [H, W] = size(difsnr);  
            
            surface(surface>(H-L-1))=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,sclera_upper] = max(above_sclera);
            
            sclera_upper = surface+sclera_upper;
            
            
            L = 10;
            surface = round(sclera_upper);
            surface((surface+L)<self.sclera) = self.sclera((surface+L)<self.sclera)-L;
            
            difsnr = diff(snr);
            [H, W] = size(difsnr);  
            
            surface(surface>(H-L-1))=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,self.sclera] = min(above_sclera);            

            self.sclera = surface+self.sclera;
            self.sclera = medfilt1(gather(self.sclera),10);
            self.sclera = (smooth(self.sclera,10))';
            self.sclera = round(self.sclera);            
             
        end   
        
        
        
        
        
        
        
        
        function ex_vivo_sclera_surface_setup(self,snr)
            
            [H, W] = size(snr);    
            snrdB = 20*log10(snr);
            snrdB(snrdB>30)=30;
            snrdB(snrdB<20)=0;
            blur = imgaussfilt(snrdB,10);
            [~,line] = max(blur);
            L = 40;   
            surface = line - round(L);
            snrdB = imgaussfilt(snrdB,2);
            difsnr = diff(snrdB);
            [H, W] = size(difsnr);      
          
            surface(surface>(H-L-1))=1;
            surface(surface<1) = 1;
            
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,ss] = max(above_sclera);
            ss = medfilt1(gather(ss),10);
            ss = (smooth(ss,5))';
            
            
            self.flat_surface = round(surface+ss-1);
            
        end
        

        
        function dataf = fitting_filter(self,data)
            
            [p,~,mu] = polyfit(1:length(data),data,4);
            dataf = polyval(p,1:length(data),[],mu);
            
        end
       
        function intelligent_flatten_setup(self,dop,snr)
            [H, W] = size(snr);

            mask = (snr>1)&(dop>0.5);
            [~,rough_surface] = max(mask>0);
            rough_surface = round(medfilt1(gather(rough_surface),10));
            rough_surface(rough_surface>(H-16))=1;
            
            surface_sub_X = repmat(rough_surface,15,1)+repmat((1:15)',1,W);
            surface_sub_Y = repmat(1:W,15,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            upper_retina = snr(indxp);
            [~,sleek_surface] = max(upper_retina);
            sleek_surface = sleek_surface +rough_surface;

            [flat_curve,is_signal] = self.extract_curvature(snr);
            
            edge_sig = diff(is_signal>0.01);
            edge_sig(edge_sig<0) = 0;
            level_sig = cumsum([0 edge_sig]);
            l2 = quantile(level_sig,0.75);
            if (l2 == 1) || (l2 == 0)
                shift = median(flat_curve(is_signal>0)-sleek_surface(is_signal>0),2);
                self.flat_surface = round(flat_curve-shift);
            elseif l2 > 1
            %l2 = quantile(level_sig,0.75);
                self.flat_surface = flat_curve;
                shift = median(flat_curve((is_signal>0)&(level_sig<=1))-sleek_surface((is_signal>0)&(level_sig<=1)),2);
                self.flat_surface(level_sig<=1) = round(flat_curve(level_sig<=1)-shift);
                for Lindex = 2:max(level_sig)  
                    shift = median(flat_curve((is_signal>0)&(level_sig>=Lindex))-sleek_surface((is_signal>0)&(level_sig>=Lindex)),2);
                    self.flat_surface(level_sig>=Lindex) = round(flat_curve(level_sig>=Lindex)-shift);        
                end
            end
            
%             surface(abs(surface-self.flat_surface)>40) = self.flat_surface(abs(surface-self.flat_surface)>40);
%             self.flat_surface = (surface+self.flat_surface)/2;
            
%             b_surface = [ones(1,50)*self.flat_surface(1) self.flat_surface ones(1,50)*self.flat_surface(end)];
%             self.flat_surface = round(conv(b_surface,self.line,'same'));
%             self.flat_surface = self.flat_surface(51:end-50);
            
            self.flat_surface(self.flat_surface<1)=1;
            self.flat_surface(isnan(self.flat_surface))=1;
            self.flat_surface = self.inter_frame_smooth(3);
            
        end
        
        function smf = inter_frame_smooth(self,N)
            dim = size(self.flat_surface_array);
            if dim(1) > N
                self.flat_surface_array(1,:) = [];
            end
            self.flat_surface_array = cat(1,self.flat_surface_array,self.flat_surface);
            smf = round(mean(self.flat_surface_array,1));
        end
        
        
        function image = test_flatten(self,image)
            [H, W] = size(image);

            if W ~= length(self.flat_surface)
                surface =  round(interp1(1:length(self.flat_surface),self.flat_surface,linspace(1,length(self.flat_surface),W)));
            else
                surface = self.flat_surface;
            end
            surface(isnan(surface))= round(mean(surface,'omitnan'));
            surface(surface<1) = 1;
            surface(surface>H) = H;
            indxp = sub2ind([H,W],surface,1:W);   
            image(indxp) = 0;     
            
            if W ~= length(self.RPE)
                surface =  round(interp1(1:length(self.RPE),self.RPE,linspace(1,length(self.RPE),W)));
            else
                surface = self.RPE;
            end
            surface(isnan(surface))= round(mean(surface,'omitnan'));
            surface(surface<1) = 1;
            surface(surface>H) = H;
            indxp = sub2ind([H,W],surface,1:W);   
            image(indxp) = 0;       
            
            
            if W ~= length(self.sclera)
                surface =  round(interp1(1:length(self.sclera),self.sclera,linspace(1,length(self.sclera),W)));
            else
                surface = self.sclera;
            end
            surface(isnan(surface))= round(mean(surface,'omitnan'));
            surface(surface<1) = 1;
            surface(surface>H) = H;
            indxp = sub2ind([H,W],surface,1:W);   
            image(indxp) = 0;             
            
%             if W ~= length(self.ex_vivo_surface)
%                 surface =  round(interp1(1:length(self.ex_vivo_surface),self.ex_vivo_surface,linspace(1,length(self.ex_vivo_surface),W)));
%             else
%                 surface = self.ex_vivo_surface;
%             end
%             surface(isnan(surface))= round(mean(surface,'omitnan'));
%             surface(surface<1) = 1;
%             surface(surface>H) = H;
%             indxp = sub2ind([H,W],surface,1:W);   
%             image(indxp) = 0;              
%             

        end
        
        
        
        
        
        
        function  imagef = flatten_image(self,image)
            [H, W] = size(image);
            if W ~= length(self.flat_surface)
                surface =  round(interp1(1:length(self.flat_surface),self.flat_surface,linspace(1,length(self.flat_surface),W)));
            else
                surface = self.flat_surface;
            end
            surface(isnan(surface))= round(mean(surface,'omitnan'));
            surface(surface<1) = 1;
            surface(surface>H) = H;
            D = max(surface);
            surface_sub_X = repmat(surface,H,1)+repmat((1:H)',1,W);
            surface_sub_Y = repmat(1:W,H,1);            
            indxp = sub2ind([H+D,W],surface_sub_X,surface_sub_Y);    
            Imex = cat(1,image,zeros(D,W));
            imagef = Imex(indxp); 
        end          
                  
        
        
        
        function l = get_layers(self)
            l = cat(1,self.flat_surface, self.RPE, self.sclera);
        end
        
                     

        
        
        function PSTriple_Manager = flatten_Mueller(self,PSTriple_Manager)

            [St,H, W] = size(PSTriple_Manager.Mmueller);  
            Sres = PSTriple_Manager.Mmueller;
            if W ~= length(self.flat_surface)
                surface =  round(interp1(1:length(self.flat_surface),self.flat_surface,linspace(1,length(self.flat_surface),W)));
            else
                surface = self.flat_surface;
            end
            surface(isnan(surface))= round(mean(surface,'omitnan'));
            surface(surface<1) = 1;
            surface(surface>H) = H;

            D = max(surface);
            surface_sub_X = repmat(surface,H,1)+repmat((1:H)',1,W);
            surface_sub_Y = repmat(1:W,H,1);            
            indxp = sub2ind([H+D,W],surface_sub_X,surface_sub_Y);     
            
            for index =1:St
                Sex = cat(1,squeeze(Sres(index,:,:)),zeros(D,W));
                Sres(index,:,:) = Sex(indxp);
            end
            
            PSTriple_Manager.Mmueller = Sres;
            
        end        
        
        
        function imagef = reverse_flatten_image(self,image)
            [H, W] = size(image);
            if W ~= length(self.flat_surface)
                surface =  round(interp1(1:length(self.flat_surface),self.flat_surface,linspace(1,length(self.flat_surface),W)));
            else
                surface = self.flat_surface;
            end
            surface(isnan(surface))= round(mean(surface,'omitnan'));
            surface(surface<1) = 1;
            surface(surface>H) = H;
            D = max(surface);
            surface_sub_X = D-repmat(surface,H,1)+repmat((1:H)',1,W);
            surface_sub_Y = repmat(1:W,H,1);            
            indxp = sub2ind([H+D,W],surface_sub_X,surface_sub_Y);    
            Imex = cat(1,zeros(D,W),image);
            imagef = Imex(indxp); 
        end
        
        
        
        
   
        

        
        function [flat,is_signal] =  extract_curvature(self,I)
            dim = size(I);
            line = [1 1]./2;
            frame = imfilter(I,line);     
            frame = imresize(frame,[dim(1),dim(2)/2]);

            frame = frame - mean(frame);
            frame2 = circshift(frame,1,2);
            framearray = repmat(frame,1,1,41);
            for findex = 1:41
                framearray(:,:,findex) = circshift(framearray(:,:,findex),findex-21,1);
            end
            xc = squeeze(sum(abs(framearray+frame2),1));
            [value,shi] = max(xc,[],2);
            flat = cumsum(smooth(shi-21,2));
            is_signal = value>750;
            flat = mean(flat(is_signal))-flat;
            flat =  interp1(1:2:dim(2),flat,1:dim(2),'linear','extrap');   
            is_signal =  interp1(1:2:dim(2),double(is_signal),1:dim(2),'linear','extrap');   
            
        end 
        
        
        
        function nim = avg_imresize(self,im,newW)
            dim = size(im);
            kernal_W = round(dim(2)./newW);
            k = ones(1,kernal_W)./kernal_W;
            im = imfilter(im,k);
            nim = imresize(im,[dim(1) newW]);
        end

        
    end
end

