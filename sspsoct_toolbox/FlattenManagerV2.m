classdef FlattenManagerV2 < handle

    properties
        flat_surface
        line = ones(1,40)/40;
        retina_surface = [];
        RPE = [];
        sclera =[];
        sclera_upper_surface=[];
        sclera_lower_surface=[];        
        ex_vivo_surface = [];
        
        width
        height
        
        low_choroid
        
    end
    
    methods
        function self = FlattenManagerV2()
        
        end
         
        function surface_flatten_setup(self,dop,snr)
            
            [H, W] = size(snr);
            
            self.width = W;
            self.height = H;
            
            snr = (imgaussfilt(snr,2));
            difsnr = (diff((snr)));
            difsnr = (imgaussfilt(difsnr,1.5));
            %%mask = (snr>1.5)&(dop>0.5);

            mask = (snr>8)&(dop>0.9); % this is for making supplementary figure

            
            [~,rough_surface] = max(mask>0);
            
            %rough_surface = round(medfilt1(gather(rough_surface),10));
            rough_surface = rough_surface-10;
            L = 20;
            
            rough_surface(rough_surface>(H-L-1))=1;
            rough_surface(rough_surface<1)=1;
            surface_sub_X = repmat(rough_surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H-1,W],surface_sub_X,surface_sub_Y);  
            upper_retina = difsnr(indxp);

%             % generate horizontal edge emphasis kernel
%             h = fspecial('sobel');
%             % invert kernel to detect vertical edges
%             h = h';
%             upper_retina = imfilter(upper_retina,h);
       
            [~,sleek_surface] = max(upper_retina);
            sleek_surface = sleek_surface +rough_surface+2;
            %sleek_surface = round(medfilt1(gather(sleek_surface),10));
            
            dataf = (smooth(sleek_surface',5))';
            
            self.flat_surface = dataf;
            self.flat_surface(self.flat_surface<1)=1;
            self.flat_surface(isnan(self.flat_surface))=1;
            self.flat_surface = round(self.flat_surface);
        end
        
        
        function GP_find_RPE(self,snr,rough_estimator)

            snr = imgaussfilt(log(snr),5);
            difsnr = diff(snr);
            [H, W] = size(difsnr);            
            surface = (self.flat_surface);

            L = 20; 
            rough_estimator = rough_estimator-7;
            surface = surface + rough_estimator;
            surface(surface>(H-L-1))=1;
            surface(surface<1)=1;
            
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            around_RPE = difsnr(indxp);
            [~,self.RPE] = max(around_RPE);
            
            self.RPE = surface +self.RPE;
            self.RPE = round(medfilt1(gather(self.RPE),10))+4;
            
        end
        
        function AGP_find_sclera(self,snr,ret)
            
            [H, W] = size(snr);      
            snr = imgaussfilt(log(snr),2);
            ret = imgaussfilt(ret.*(snr>1.5),5);
            retmask = ret>0.06;
            
            surface = self.RPE+10;
            L = 40;

            surface(surface>(H-L-1))=1;
            surface(surface<1)=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            belowRPE = ret(indxp);
            [~,self.sclera] = max(belowRPE);
            self.sclera = surface+self.sclera;
            self.sclera = medfilt1(gather(self.sclera),10);
            self.sclera = (smooth(self.sclera,10))';
            self.sclera = round(self.sclera);
            
            L = 30;
            difsnr = diff(snr);
            [H, W] = size(difsnr);  
            surface = round(self.RPE);
            
            surface(surface>(H-L-1))=1;
            
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            
            [~,self.low_choroid] = min(above_sclera);
            self.low_choroid = surface+self.low_choroid;
            self.low_choroid(self.low_choroid>self.sclera) = self.RPE(self.low_choroid>self.sclera)+3;
            
            self.low_choroid = medfilt1(gather(self.low_choroid),10);
            self.low_choroid = (smooth(self.low_choroid,10))';
            self.low_choroid  = round(self.low_choroid);  
            
            
            
            
            L = 40;
            surface = round(self.sclera)-40;
            surface(surface<self.low_choroid) = self.low_choroid(surface<self.low_choroid);
%             surface(surface<low_choroid) = low_choroid(surface<low_choroid);
            
            difsnr = diff(snr);
            [H, W] = size(difsnr);  
            
            surface(surface>(H-L-1))=1;
            surface(surface<1)=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,self.sclera_upper_surface] = max(above_sclera);

            self.sclera_upper_surface = surface+self.sclera_upper_surface;
            self.sclera_upper_surface(self.sclera_upper_surface>self.sclera) = self.sclera(self.sclera_upper_surface>self.sclera);
            
            self.sclera_upper_surface = medfilt1(gather(self.sclera_upper_surface),10);
            self.sclera_upper_surface = (smooth(self.sclera_upper_surface,10))';
            self.sclera_upper_surface = round(self.sclera_upper_surface);  
            
            
            
            
            
            L = 20;
            surface = round(self.sclera);
            
%             difsnr = diff(snr);
%             [H, W] = size(difsnr);  
            
            surface(surface>(H-L-1))=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,self.sclera_lower_surface] = min(above_sclera);
            
            self.sclera_lower_surface = surface+self.sclera_lower_surface;
            self.sclera_lower_surface = medfilt1(gather(self.sclera_lower_surface),10);
            self.sclera_lower_surface = (smooth(self.sclera_lower_surface,10))';
            self.sclera_lower_surface = round(self.sclera_lower_surface);  

        end       
          

        function PGP_find_sclera(self,snr,ret)
            
            [H, W] = size(snr);      
            snr = log(imgaussfilt((snr),2));
            ret = imgaussfilt(ret.*(snr>1.0),5);
            %snr = log(snr);
            retmask = ret>0.06;
            
            surface = self.RPE+5;
            L = 60;

            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            belowRPE = ret(indxp);
            [~,self.sclera] = max(belowRPE);
            self.sclera = surface+self.sclera;
            self.sclera = medfilt1(gather(self.sclera),30);
            self.sclera = (smooth(self.sclera,30))';
            self.sclera = round(self.sclera);


            L = 35;
            surface = round(self.sclera)-L;
            
            C =5;
            surface(surface<(self.RPE+C)) = self.RPE(surface<(self.RPE+C))+C;

            
            difsnr = diff(snr);
            [H, W] = size(difsnr);  
            
            surface(surface>(H-L-1))=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,low_choroid] = max(above_sclera);
            
            low_choroid = surface+low_choroid;
            low_choroid = round((smooth(low_choroid,10))');
            self.low_choroid = round(low_choroid);
            
            
            L = 15;
            surface = low_choroid;
            aaa = self.sclera;
            surface((surface+L)>aaa) = aaa((surface+L)>aaa)-L;
            surface(surface<low_choroid) = low_choroid(surface<low_choroid);
            
            
            
            surface(surface>(H-L-1))=1;
            surface(surface<1)=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,self.sclera_upper_surface] = min(above_sclera);     
            self.sclera_upper_surface = surface+self.sclera_upper_surface;
            self.sclera_upper_surface = medfilt1(gather(self.sclera_upper_surface),10);
            self.sclera_upper_surface = (smooth(self.sclera_upper_surface,10))';
            self.sclera_upper_surface = round(self.sclera_upper_surface);                
            

            L =50;
            surface = round(self.sclera);
            surface((surface+L)<self.sclera) = self.sclera((surface+L)<self.sclera)-L;
            
            difsnr = diff(snr);
            [H, W] = size(difsnr);  
            
            surface(surface>(H-L-1))=1;
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            above_sclera = difsnr(indxp);
            [~,self.sclera_lower_surface] = min(above_sclera);            

            self.sclera_lower_surface = surface+self.sclera_lower_surface;
            self.sclera_lower_surface = medfilt1(gather(self.sclera_lower_surface),10);
            self.sclera_lower_surface = (smooth(self.sclera_lower_surface,10))';
            self.sclera_lower_surface = round(self.sclera_lower_surface);            
             
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
%             
            if W ~= length(self.low_choroid)
                surface =  round(interp1(1:length(self.low_choroid),self.low_choroid,linspace(1,length(self.low_choroid),W)));
            else
                surface = self.low_choroid;
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
            
            if W ~= length(self.sclera_upper_surface)
                surface =  round(interp1(1:length(self.sclera_upper_surface),self.sclera,linspace(1,length(self.sclera_upper_surface),W)));
            else
                surface = self.sclera_upper_surface;
            end
            surface(isnan(surface))= round(mean(surface,'omitnan'));
            surface(surface<1) = 1;
            surface(surface>H) = H;
            indxp = sub2ind([H,W],surface,1:W);   
            image(indxp) = 0;      
            
            if W ~= length(self.sclera_lower_surface)
                surface =  round(interp1(1:length(self.sclera_lower_surface),self.sclera,linspace(1,length(self.sclera_lower_surface),W)));
            else
                surface = self.sclera_lower_surface;
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
        
        function PSTriple_Manager = reverse_flatten_Mueller(self,PSTriple_Manager)
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
            surface_sub_X = D-repmat(surface,H,1)+repmat((1:H)',1,W);
            surface_sub_Y = repmat(1:W,H,1);            
            indxp = sub2ind([H+D,W],surface_sub_X,surface_sub_Y);   


            for index =1:St
                Sex = cat(1,ones(D,W),squeeze(Sres(index,:,:)));
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
        

        function thickness = get_choroid_thickness(self)
            thickness = self.sclera_upper_surface-self.RPE;
        end
        function thickness = get_sclera_thickness(self)
            thickness = self.sclera_lower_surface-self.sclera_upper_surface; 
        end  
        
        
        function m = get_mask(self,layer)            
            H = self.height;
            W = self.width;
            image = ones(H,W);

            if layer == 1
                b = self.RPE-self.flat_surface;
                b = repmat(b,H,1)+repmat((1:H)',1,W);
                b(b<1) = 1;
                b(b>H) = H;               
                indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
                image(indxp) = 0;    
            end
            if layer == 2
                b = self.RPE-self.flat_surface;
                b = repmat(b,H,1)+repmat((1:H)',1,W);
                b(b<1) = 1;
                b(b>H) = H;               
                indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
                image(indxp) = 0;    
                image = 1-image;
                
                b = self.sclera_upper_surface-self.flat_surface;
                b = repmat(b,H,1)+repmat((1:H)',1,W);
                b(b<1) = 1;
                b(b>H) = H;               
                indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
                image(indxp) = 0;                    
            end
            if layer == 3
                b = self.sclera_upper_surface-self.flat_surface;
                b = repmat(b,H,1)+repmat((1:H)',1,W);
                b(b<1) = 1;
                b(b>H) = H;               
                indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
                image(indxp) = 0;    
                image = 1-image;
                
                b = self.sclera_lower_surface-self.flat_surface; 
                b = repmat(b,H,1)+repmat((1:H)',1,W);
                b(b<1) = 1;
                b(b>H) = H;               
                indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
                image(indxp) = 0;                    
            end
            
            m = image;
        
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

