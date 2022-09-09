classdef FlattenManagerForHuman < handle

    properties
        flat_surface
        line = ones(1,40)/40;
        flat_surface_array = [];
        
        flat_surface_array2 = [];
        flat_surface_array3 = [];
        
        RPE = [];
        sclera=[];
        RNFL=[];
        ex_vivo_surface = [];
        height;
        width;
        inner_sclera = [];
        
        dop;
        CalStru
        frameindex
    end
    
    methods
        function self = FlattenManagerForHuman()
        end
        
        
        function useCalStru(self,CalStru,frameindex)
            self.CalStru = CalStru;
            self.frameindex = frameindex;
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

 
        
        function surface_flatten_setup(self,dop,snr)
             
            [H, W] = size(snr);
             
            self.dop = dop;
            
            if isfield(self.CalStru,'RNFLface')
                sm = SurfaceManager(self.CalStru.RNFLface);
                self.flat_surface = sm.get_line(self.frameindex);
                
               if W ~= length(self.flat_surface )
                    self.flat_surface =  round(interp1(1:length(self.flat_surface),self.flat_surface ,linspace(1,length(self.flat_surface),W)));
               end
               return
            end

            dop(:,end-30:end) =  dop(:,end-30:end)+0.2;
            
            dop(:,1:30) =  dop(:,1:30)+0.2;
            dop(dop>1)=1;
            snr = imgaussfilt((snr),3);
            %mask = (snr>2)&(dop>0.4);
            %mask = (snr>1.7)&(dop>0.4);
            
            
            mask = (snr>1.7)&(dop>0.3); % this is used for most data
%             %mask = ((20*log10(snr))>1.4)&(dop>0.6);  % I cannot remember
%             %what is this used
%             SE = strel('disk',15);
%             
%             mask = imopen(mask,SE);            
%             mask = imclose(mask,SE); 
%            mask(1:30,:)=0;
            
            
            [~,rough_surface] = max(mask>0);
%             rough_surface(rough_surface==1) = mean(rough_surface(~(rough_surface==1)));
%             rough_surface(isnan(rough_surface)) = 100;
            
            rough_surface = round(medfilt1(gather(rough_surface),30,'omitnan','truncate'));
            
            
            L = 10;
            difsnr = diff(log(snr));
            rough_surface(rough_surface>(H-L-1))=1;
            
            
            surface_sub_X = repmat(rough_surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);
            aaa = log(snr);
            upper_retina = aaa(indxp);
            
            
%             % generate horizontal edge emphasis kernel
%             h = fspecial('sobel');
%             % invert kernel to detect vertical edges
%             h = h';
%             upper_retina = imfilter(upper_retina,h);
%             
            
            [~,sleek_surface] = max(upper_retina);
            sleek_surface = sleek_surface +rough_surface-5;
            sleek_surface = round(medfilt1(gather(sleek_surface),10,'omitnan','truncate'));
            
            dataf = (smooth(sleek_surface',10))';
            
            self.flat_surface = dataf;
            self.flat_surface(self.flat_surface<1)=1;
            self.flat_surface(isnan(self.flat_surface))=1;
            
            self.flat_surface = round(self.flat_surface);
            %self.flat_surface = self.inter_frame_smooth(5);
       
        end
        
     
        
        
        
        function find_RPE(self,snr,rough_estimator)
            
            
            [H, W] = size(snr);
            self.height = H;
            self.width = W;
            
            if isfield(self.CalStru,'RPEface')
                sm = SurfaceManager(self.CalStru.RPEface);
                self.RPE = sm.get_line(self.frameindex);
                if W ~= length(self.RPE )
                    self.RPE  =  round(interp1(1:length(self.RPE ),self.RPE ,linspace(1,length(self.RPE ),W)));
                end  
                return
            end

            snr = imgaussfilt(log(snr),5);
            difsnr = diff(snr);
            [H, W] = size(difsnr);            
            surface = self.flat_surface;
            
            rough_estimator = ones(size(surface)).*rough_estimator;
            
            rough_estimator(surface<50)=5;
            
            surface = surface+rough_estimator-5;
            
            
            

            L = 60;
            surface(surface>(H-L-1))=1;
            surface(surface<1)=1;
            
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  

            around_RPE = difsnr(indxp);
            [~,self.RPE] = max(around_RPE);
            
            self.RPE = surface +self.RPE;
            
            self.RPE = medfilt1(gather(self.RPE),30,'omitnan','truncate');
            self.RPE = round((smooth((self.RPE)',10))');
        end
        
        function m = RNFL_mask_after_flatten(self,snr)
            
            [H, W] = size(snr);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            m = self.get_RNFL_mask(); 
            snr = imgaussfilt(log(snr),2);
            snr = snr.*m;
            
            maxI = max(snr(:));
            level = maxI/5;
            snr(snr<level) = 0;
            snr = snr-level;
            snr(snr<level) = 0;
            
            difsnr = diff(snr);
            [~,self.RNFL] = max(difsnr);
            
            self.RNFL = round(medfilt1(gather(self.RNFL),10,'omitnan','truncate'));
            
            image = ones(H,W);

            b = self.RNFL;
            b = repmat(b,H,1)+repmat((1:H)',1,W);
            b(b<1) = 1;
            b(b>H) = H;               
            indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
            image(indxp) = 0;    
           
            m = image; 

        end   
        
        function find_sclera(self,snr,ret,dop)
            
            if isfield(self.CalStru,'Scleraface')
                
                sm = SurfaceManager(self.CalStru.Scleraface);
                self.sclera = sm.get_line(self.frameindex,2.5);
                [H, W] = size(snr);                 
                
                if ~isempty(self.sclera) 

                    if W ~= length(self.sclera)
                        self.sclera =  round(interp1(1:length(self.sclera),self.sclera,linspace(1,length(self.sclera),W)));
                    end
                    return
                else
                end
            end

            
            [H, W] = size(snr);      
            snr = imgaussfilt((snr),3);
            
            dop(:,end-30:end) = dop(:,end-30:end)+0.1;
            ret = imgaussfilt(ret.*(snr>1.2).*(dop>0.6),5);
            
            snr = diff(log(snr));
            H = H-1;
            
            %ret = diff(ret);
            
            ret = ret(2:end,:);
            ret(1:50,:)=0;
            ret(end-50:end,:)=0;
            retf = imgaussfilt((ret),30);
            
            
            surface = self.RPE;
            L = 200;
                    
            surface(surface>(H-L-1))=H-L-1;
            surface(surface<1) = 1;
            
            %figure;imagesc(ret);hold on;plot(surface);   
            

            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            ret_belowRPE = retf(indxp);  
            %figure;imagesc(ret_belowRPE);   
            
            
            [~,middle_sclera] = max(ret_belowRPE); %
 
            
            middle_sclera = medfilt1(gather(middle_sclera),30,'omitnan','truncate');
            middle_sclera = round((smooth(middle_sclera,30))')+(self.RPE);
            
            %figure;imagesc(ret);hold on;plot(middle_sclera);   

            
            sc_mask = self.get_mask(middle_sclera,snr,120);
            
            snr(~sc_mask)=0; 
            
            middle_sclera(middle_sclera<1) = 1;
            middle_sclera(middle_sclera>H) = H;
            indxp = sub2ind([H,W],middle_sclera,1:W);   
            snr(indxp) = 0;     
             
            snr_bw = snr>quantile(snr(:),0.7);
            
            SE = strel('rectangle',[5 2]);
            
            snr_bw = imopen(snr_bw,SE);
            snr_bw = imclose(snr_bw,SE);

            
            retmask = 0.08;
            ret_bw = ret>retmask;
            
            %figure;imshow(ret_bw)
                     %ret_bw = imopen(ret>retmask,SE);
            %SE = strel('rectangle',[8 3]);            
            %figure;imshow(ret_bw)
            
            surface = self.RPE;
            L = 130;
            
            surface(surface>(H-L-1))=H-L-1;
            surface(surface<1) = 1;
            
            
            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            belowRPE_ret = ret_bw(indxp);
            belowRPE_int = snr_bw(indxp);

            %sclera_int = sclera_mask.*belowRPE_int;
            belowRPE_int = flipud(belowRPE_int);
            %belowRPE_ret = flipud(belowRPE_ret);
            
            [~,aaa] = max(belowRPE_int);
            [~,bbb] = max(belowRPE_ret);
            aaa = medfilt1(gather(aaa),31,'omitnan','truncate');
            bbb = L-bbb;
            
            %figure;imagesc(snr);hold on;plot(aaa);   
            %figure;imagesc(ret);hold on;plot(bbb);   
%             aaa = medfilt1(gather(aaa),31);
%             bbb = medfilt1(gather(bbb),31);
            
           %figure;plot(aaa);hold on;plot(bbb);
            
            ccc = (aaa+bbb)/2;
            
            self.sclera = ccc;%min(aaa,bbb);
            self.sclera(1:30) = aaa(1:30);
            self.sclera(end-30:end) = aaa(end-30:end);
            
            
            self.sclera = (self.RPE)+ (L-self.sclera);
            %figure;plot(self.sclera);hold on;
            onh = self.find_onh_region(self.sclera);
            
            %figure;plot(self.sclera);hold on;plot(onh*400)

            onh(1:5) = 1;
            onh(end-5:end) = 1;
            
            %self.sclera = medfilt1(gather(self.sclera),30);
            self.sclera = self.inter_frame_smooth(self.sclera,3);
            self.sclera = self.fitting_filter(self.sclera,1,onh);
            
            self.sclera = self.inter_frame_smooth3(self.sclera,3);
            %plot(self.sclera)
            %[~,self.sclera] = max(aaa,bbb);
            %self.sclera = (self.RPE)+ (L-self.sclera);
%             self.sclera = medfilt1(gather(self.sclera),10);
%             self.sclera = (smooth(self.sclera,10))';

            self.sclera = round(self.sclera);

        end       
        
        
        
        

        function inner_mean = find_inner_sclera(self,snr,ret,phi)
            
            [H, W] = size(snr);  
            surface = self.sclera;
            L = 100;
            
            surface(surface>(H-L-1))=H-L-1;
            surface(surface<1) = 1;

            surface_sub_X = repmat(surface,L,1)+repmat((1:L)',1,W);
            surface_sub_Y = repmat(1:W,L,1);        
            
            indxp = sub2ind([H,W],surface_sub_X,surface_sub_Y);  
            
            belowSclera_ret = gather(ret(indxp));
            belowSclera_phi = gather(phi(indxp));

            cpl = belowSclera_ret.*exp(1i*belowSclera_phi);
            inner_mean = (mean(cpl(1:20,:)))./abs(mean(cpl(1:20,:)));
            cpl = cpl./(inner_mean);
            
            %cpl_plus = circshift(cpl,1,1);
            diffphi1 = abs(angle(cpl));
            
            diffphi1(1:10,:)=0;
            
            diffphi1 = imgaussfilt(diffphi1,5);
            
            dp_bw = diffphi1>1;
            
            [~,aaa] = max(dp_bw);

            self.inner_sclera = (self.sclera)+aaa;  
            self.inner_sclera = medfilt1(gather(self.inner_sclera),50,'omitnan','truncate');
            self.inner_sclera = (smooth(self.inner_sclera,50))';            
            %figure;plot(self.inner_sclera);hold on;plot(self.sclera);
            self.inner_sclera = self.inter_frame_smooth2(self.inner_sclera,5);
            
            self.inner_sclera(self.inner_sclera<gather(self.sclera))=gather(self.sclera(self.inner_sclera<gather(self.sclera)));
            
            self.inner_sclera = round(self.inner_sclera);

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
            ss = medfilt1(gather(ss),10,'omitnan','truncate');
            ss = (smooth(ss,5))';
            
            
            self.flat_surface = round(surface+ss-1);
            
        end
        

        
        function dataf = fitting_filter(self,data,cut,exclude)
            

            X = 1:length(data);
            Y = data;
            
            X(exclude) = [];
            Y(exclude) = [];
            
            
            X = X(cut:end-cut);
            Y = Y(cut:end-cut);
            [p,~,mu] = polyfit(X,Y,5);
            dataf = polyval(p,1:length(data),[],mu);
            
        end
       
        function onh = find_onh_region(self,data)
            
            
            %data = self.RPE;
            datamf = medfilt1(gather(data),11,'omitnan','truncate');
            
            x = 10:length(data)-10;
            y = datamf(10:end-10);

            
            [p,~,mu] = polyfit(x,y,3);
            dataf = polyval(p,1:length(data),[],mu);            
            
            error = (data-dataf);
            
            onh = error(20:end-20) > 15;
            
            x(onh) = [];
            y(onh) = [];
            
            
            [p,~,mu] = polyfit(x,y,3);
            dataf = polyval(p,1:length(data),[],mu);                   
            
            %figure;plot(dataf);hold on;plot(data);
            
            
            error = (data-dataf);
            
            onh = error > 25;            
%             onh(1:30)=0;
%             onh(end-30:end)=0;            
            
            SE = [1 1 1 1 1 1 1 1 1 1];
            onh = imdilate(onh,SE);
            


        end
        
        
        function smf = inter_frame_smooth(self,surface,N)
            dim = size(self.flat_surface_array);
            if dim(1) > N
                self.flat_surface_array(1,:) = [];
            end
            self.flat_surface_array = cat(1,self.flat_surface_array,surface);
            smf = (mean(self.flat_surface_array,1));
        end
        
        
        function smf = inter_frame_smooth2(self,surface,N)
            dim = size(self.flat_surface_array2);
            if dim(1) > N
                self.flat_surface_array2(1,:) = [];
            end
            self.flat_surface_array2 = cat(1,self.flat_surface_array2,surface);
            smf = (mean(self.flat_surface_array2,1));
        end
                
        function smf = inter_frame_smooth3(self,surface,N)
            dim = size(self.flat_surface_array3);
            if dim(1) > N
                self.flat_surface_array3(1,:) = [];
            end
            self.flat_surface_array3 = cat(1,self.flat_surface_array3,surface);
            smf = (mean(self.flat_surface_array3,1));
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
            
%             if W ~= length(self.inner_sclera)
%                 surface =  round(interp1(1:length(self.inner_sclera),self.inner_sclera,linspace(1,length(self.inner_sclera),W)));
%             else
%                 surface = self.inner_sclera;
%             end
%             surface(isnan(surface))= round(mean(surface,'omitnan'));
%             surface(surface<1) = 1;
%             surface(surface>H) = H;
%             indxp = sub2ind([H,W],surface,1:W);   
%             image(indxp) = 0;                 
            
            
            
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
                  
        function  imagef = flatten_image_at_RPE(self,image)
            [H, W] = size(image);
            if W ~= length(self.RPE)
                surface =  round(interp1(1:length(self.RPE),self.RPE,linspace(1,length(self.RPE),W)));
            else
                surface = self.RPE;
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
        
        
        function  imagef = flatten_image_at_sclera(self,image)
            [H, W] = size(image);
            if W ~= length(self.sclera)
                surface =  round(interp1(1:length(self.sclera),self.sclera,linspace(1,length(self.sclera),W)));
            else
                surface = self.sclera;
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
            l = cat(1,self.flat_surface, self.RPE, self.sclera,self.inner_sclera);
        end
        
        
        function t = get_choroid_thickness(self)
            t =  self.sclera-self.RPE;
        end
        
        function t = get_inner_retina_thickness(self)
            t = self.RPE - self.flat_surface;
        end
        function t = get_inner_sclera_thickness(self)
            t =  self.inner_sclera-self.sclera;
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
        
        
        
        function m = get_RNFL_mask(self)            
            H = self.height;
            W = self.width;
            image = ones(H,W);

            b = self.RPE-15-self.flat_surface;
            b = repmat(b,H,1)+repmat((1:H)',1,W);
            b(b<1) = 1;
            b(b>H) = H;               
            indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
            image(indxp) = 0;    
           
            m = image;
        
        end
        
        
        function m = get_inner_sclera_mask(self)            
            H = self.height;
            W = self.width;
            image = ones(H,W);

            b = self.sclera - 30;
            b = repmat(b,H,1)+repmat((1:H)',1,W);
            b(b<1) = 1;
            b(b>H) = H;               
            indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
            image(indxp) = 0;    
            image = 1-image;

            b = self.sclera + 30;
            b = repmat(b,H,1)+repmat((1:H)',1,W);
            b(b<1) = 1;
            b(b>H) = H;               
            indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
            image(indxp) = 0;   

            m = image;
        
        end
        
        function m = get_mask(self,surface,image,L)            
            [H,W] = size(image);
            image = ones(H,W);

            b = surface-L;
            b = repmat(b,H,1)+repmat((1:H)',1,W);
            b(b<1) = 1;
            b(b>H) = H;               
            indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
            image(indxp) = 0;    
            image = 1-image;

            b = surface;
            b = repmat(b,H,1)+repmat((1:H)',1,W);
            b(b<1) = 1;
            b(b>H) = H;               
            indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
            image(indxp) = 0;   

            m = image;
        
        end
        
        
        
        
        function m = get_true_inner_sclera_mask(self)            
            H = self.height;
            W = self.width;
            image = ones(H,W);

            b = self.sclera;
            b = repmat(b,H,1)+repmat((1:H)',1,W);
            b(b<1) = 1;
            b(b>H) = H;               
            indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
            image(indxp) = 0;    
            image = 1-image;

            b = self.inner_sclera
            b = repmat(b,H,1)+repmat((1:H)',1,W);
            b(b<1) = 1;
            b(b>H) = H;               
            indxp = sub2ind([H,W],b,repmat(1:W,H,1));   
            image(indxp) = 0;   

            m = image;
        
        end        
        
        
        
        
        
        
        
        
                
        

        
    end
end

