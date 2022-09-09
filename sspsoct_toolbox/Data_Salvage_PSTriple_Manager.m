classdef Data_Salvage_PSTriple_Manager < handle
    % this class use two input system to solve retardation and optic axis
    
    properties

        dim % size of the Stokes frame
        sdim % dimension of Stokes
        bins % bin number
        pixel_depth % the real length of one pixel
        zfilt_kernal
        
        Mmueller
        I1
        I2
        I3
        
        CalStru
        
        Mcorr % symetrix correction
        PMDCorr % PMD correction
        DPMDCorrTD 
        
        
        retardation
        
        nan_mask
        
        MMc
        mask
        retina_mask
        last_surface_WW = 0;
        
        
        ret_level = 4;
        
        
        %% 
        
        order
        
        %%
        cumulative_r
        
        
        %% salvage
        s
        cs12
        cs13
        cs23
        invs
        isis
    end
    
    methods
        
        function self = Data_Salvage_PSTriple_Manager(CalStru)
            %% in this function, do orthonolize
            self.pixel_depth = CalStru.whole_depth/CalStru.NFFT;
            self.CalStru = CalStru;
            %self.zfilt_kernal = zfilt_kernal;
            if isfield(CalStru,'binning_rc')
                try
                    self.Mcorr = self.CalStru.binning_rc.Mcorr;
                    self.PMDCorr = self.CalStru.binning_rc.PMDCorr;
                    self.DPMDCorrTD = self.CalStru.binning_rc.DPMDCorrTD;   
                    self.order = self.CalStru.binning_rc.order;  
                catch
                end
            end
            
            
            
            %% salvage
            self.s = importdata('InputS.mat');
             
            s1 = self.s(1:3,:);
            s2 = self.s(4:6,:);
            s3 = self.s(7:9,:);
            cs12 = sqrt(1 - sum(s1.*s2,1));
            cs13 = sqrt(1 - sum(s1.*s3,1));
            cs23 = sqrt(1 - sum(s2.*s3,1));            
            self.cs12 = permute(cs12,[1 3 4 2]);
            self.cs13 = permute(cs13,[1 3 4 2]);
            self.cs23 = permute(cs23,[1 3 4 2]);
            
            self.invs = matrix3inv(self.s);
            self.isis = MatrixMultiply(self.invs,transpose3x3(self.invs));
            
            
        end
%% 
        function obtainSortCorrection(self,S1,S2,S3)
            
            % compute corresponding intensity signal
            I1 = sqrt(sum(S1.^2,4));
            I2 = sqrt(sum(S2.^2,4));
            I3 = sqrt(sum(S3.^2,4));

            S1n = S1./I1;
            S2n = S2./I2;
            S3n = S3./I3;

            dotProd12 = (sum(S1n.*S2n,4));
            dotProd13 = (sum(S1n.*S3n,4));
            dotProd23 = (sum(S2n.*S3n,4));

            mask1 = 10*log10(I1)>quantile(10*log10(I1(:)),0.6);
            mask2 = 10*log10(I2)>quantile(10*log10(I2(:)),0.6);
            mask3 = 10*log10(I3)>quantile(10*log10(I3(:)),0.6);
            
            mask = mask1&mask2&mask3;

            Sto = {S1,S2,S3};
            ind = 3;
            
            angle3 = acosd(dotProd12(:,:,ind));
            angle2 = acosd(dotProd13(:,:,ind));
            angle1 = acosd(dotProd23(:,:,ind));
            a1 = median(angle1(mask(:,:,ind)),'all'); 
            a2 = median(angle2(mask(:,:,ind)),'all'); 
            a3 = median(angle3(mask(:,:,ind)),'all');             

            [~,i1] = min([a1,a2,a3]);
            self.order(1) = gather(i1);
            
            ind = 5;
            
            angle3 = acosd(dotProd12(:,:,ind));
            angle2 = acosd(dotProd13(:,:,ind));
            angle1 = acosd(dotProd23(:,:,ind));
            a1 = median(angle1(mask(:,:,ind)),'all'); 
            a2 = median(angle2(mask(:,:,ind)),'all'); 
            a3 = median(angle3(mask(:,:,ind)),'all'); 
            
            [~,i2] = min([a1,a2,a3]);
            self.order(2) = gather(i2);


        end
%%
        function [S1,S2,S3] = sortS123(self,S1,S2,S3)
            Sto = {S1,S2,S3};
            S1 = Sto{self.order(1)};
            S2 = Sto{self.order(2)};
            id = [1 2 3];
            id(self.order) = [];
            S3 = Sto{id}; 
            
            
        end
        
        
        
        %% load_Stokes>>tri_input
        function load_Stokes(self,S1,S2,S3)
            
            self.dim = size(S1);
            self.sdim = length(self.dim);
            self.bins = self.dim(self.sdim-1);

            %self.checkStokes(S1,S2,S3);
            
            [S1,S2,S3] = self.sortS123(S1,S2,S3);
            %self.checkStokes(S1,S2,S3);
            
            if isfield(self.CalStru,'GPU')
                if self.CalStru.GPU == 1
                    S1 = gpuArray(single(S1));
                    S2 = gpuArray(single(S2));
                    S3 = gpuArray(single(S3));
                    self.nan_mask = gpuArray(single(self.nan_mask));
                else
                end
            else
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            S1 = permute(S1,[self.sdim, 1:self.sdim-1]);
            S2 = permute(S2,[self.sdim, 1:self.sdim-1]);
            S3 = permute(S3,[self.sdim, 1:self.sdim-1]);
            
            self.I1 = sqrt(sum(S1.^2,1));
            self.I2 = sqrt(sum(S2.^2,1));
            self.I3 = sqrt(sum(S3.^2,1));

            % compute initial c-values
            c12 = (sqrt(abs(self.I1.*self.I2 - sum(S1.*S2,1))))./self.cs12;
            c13 = (sqrt(abs(self.I1.*self.I3 - sum(S1.*S3,1))))./self.cs13;
            c23 = (sqrt(abs(self.I2.*self.I3 - sum(S2.*S3,1))))./self.cs23;

            % (arbitrarily) define their mean
            cEst = (c12 + c13 + c23)/3;

            % find how S1, S2, and S3 have to be scaled to achieve cEst
            alpha = cEst.*c23./c12./c13;
            beta = cEst.*c13./c12./c23;
            gamma = cEst.*c12./c13./c23;
            
            Msub = cat(1,S1,S2,S3);

            S1 = alpha.*S1;
            S2 = beta.*S2;
            S3 = gamma.*S3;

            I1 = sqrt(sum(S1.^2,1));
            I2 = sqrt(sum(S2.^2,1));
            I3 = sqrt(sum(S3.^2,1));

            m1 = I1./cEst;
            m2 = I2./cEst;
            m3 = I3./cEst;
            x = cat(1,m1,m2,m3);
            
    %         Msub = cat(1,S1,S2,S3);   
    %         aaa = det3x3(Msub);
            % normalized a vector
            [~,H,W,~] = size(S1);
            invsT_pm = permute(transpose3x3(self.invs),[1 3 4 2]);
            invsT_rp = repmat(invsT_pm,[1,H,W,1]);            
            s_pm = permute(self.s,[1 3 4 2]);
            s_rp = repmat(s_pm,[1,H,W,1]);       
            %ddet = S1(1,:).*S2(2,:).*S3(3,:) + S1(2,:).*S2(3,:).*S3(1,:) + S1(3,:).*S2(1,:).*S3(2,:) - S1(3,:).*S2(2,:).*S3(1,:) - S1(1,:).*S2(3,:).*S3(2,:) - S1(2,:).*S2(1,:).*S3(3,:);
            %ddet = det3x3((cat(1,S1,S2,S3)./cEst));
            ddet = det3x3(MatrixMultiply(cat(1,S1,S2,S3)./cEst,invsT_rp));
            
            [cha,cha2] = self.solve_cha(m1,m2,m3);
            cha(ddet<0) = cha2(ddet<0);
            
            cha(cha<1) = NaN;

            sha = sinh(acosh(cha));
            
            an1 = dot(invsT_rp([1 4 7],:,:,:),((x-cha)./sha),1);
            an2 = dot(invsT_rp([2 5 8],:,:,:),((x-cha)./sha),1);
            an3 = dot(invsT_rp([3 6 9],:,:,:),((x-cha)./sha),1);
            
            an = cat(1,an1,an2,an3);
%         	mD = (cha-1).*self.kron_vector_expand(an);
% 
%             mD(1,:) = mD(1,:)+1;
%             mD(5,:) = mD(5,:)+1;
%             mD(9,:) = mD(9,:)+1;
            mD = reshape(eye(3),9,1) + (cat(1,an(:,:,:,:).*an(1,:,:,:),an(:,:,:,:).*an(2,:,:,:),an(:,:,:,:).*an(3,:,:,:)).*(cha(:,:,:,:)-1));
            
            mE = MatrixMultiply(mD,s_rp);
            mE = mE + sha.*cat(1,an,an,an);
            
            mEinv = matrix3inv(mE);
            
            mR = MatrixMultiply(cat(1,S1,S2,S3)./cEst,mEinv);
            
%             mR =  MatrixMultiply(Msub./cEst,mEinv);
%             mR = round(mR.*1e5)./1e5;
%             mR = euclideanRotation(mR);
%             
%             out.D = an.*tanh(acosh(cha));
%             out.R = decomposeRot(mR);


            self.Mmueller = self.constructMueller2(cha,sha.*an,mD,mR).*cEst;
            
            self.nan_mask = ~isnan(mean(self.Mmueller,1));

        end
        
        function Mmueller = constructMueller(self,m00,D,mD,mR)
            dimM = size(D);
            dimM(1) = 16;
            Mmueller = ones(dimM);
            m00 = reshape(m00,1,[]);
            D = reshape(D,3,[]);
            mD = reshape(mD,9,[]);
            mR = reshape(mR,9,[]);
            Mmueller = reshape(Mmueller,16,[]);
            if isfield(self.CalStru,'GPU')
                if self.CalStru.GPU == 1
                    Mmueller = gpuArray(single(Mmueller));
                else
                end
            else
            end
            msub = MatrixMultiply(mR,mD);
            
            Mmueller(1,:) = m00(1,:);
            Mmueller(2,:) = dot(mR([1 4 7],:),D,1);
            Mmueller(3,:) = dot(mR([2 5 8],:),D,1);
            Mmueller(4,:) = dot(mR([3 6 9],:),D,1);
            Mmueller(5,:) = D(1,:).*m00(1,:);
            Mmueller(9,:) = D(2,:).*m00(1,:);
            Mmueller(13,:) = D(3,:).*m00(1,:);
            Mmueller(6:8,:) = msub(1:3,:);
            Mmueller(10:12,:) = msub(4:6,:);
            Mmueller(14:16,:) = msub(7:9,:);
            Mmueller = reshape(Mmueller,dimM);
            
        end
        
        function look_nan(self,m,clim)
            figure
            for bindex = 1:self.bins
                subplot(1,self.bins,bindex)
                imagesc(gather(squeeze(m(1,:,:,bindex))),clim)
            end
           
        end
        
        
        function dxdT = kron_vector_expand(self,D)
            dimD = size(D);
            D = reshape(D,3,[]);
            reshapedim = size(D);
            reshapedim(1) = 9;
            dxdT = ones(reshapedim);
            if isfield(self.CalStru,'GPU')
                if self.CalStru.GPU == 1
                    dxdT = gpuArray(single(dxdT));
                else
                end
            else
            end
            dxdT(1,:) = D(1,:).^2;
            dxdT(2,:) = D(1,:).*D(2,:);
            dxdT(3,:) = D(1,:).*D(3,:);
            dxdT(4,:) = dxdT(2,:);
            dxdT(5,:) = D(2,:).^2;
            dxdT(6,:) = D(2,:).*D(3,:);
            dxdT(7,:) = dxdT(3,:);
            dxdT(8,:) = dxdT(6,:);
            dxdT(9,:) = D(3,:).^2;
            dimD(1) = 9;
            dxdT = reshape(dxdT,dimD);
            
        end
        
        function testCorrection(self)
            mask = 1;
            self.MMc = checkSymmetric(self.MM,mask,self.CalStru.symmetric_wcorr);                 
            OA = decomposeRot(self.MMc);       
            
            dimmm = size(self.MMc);
            bsize = round(dimmm(4)/2);

            figure

            for bindex = 1:2:dimmm(4)
                for quvindex = 1:3
                    subplot(bsize,3,quvindex+((bindex-1)/2)*3)
                    imagesc(squeeze(OA(quvindex,:,:,bindex)),[-2 2])
                end
            end
            %% 
            figure
            
            for bindex = 1:2:dimmm(4)
                S_plot = self.MMc(1:3,:,:,bindex);
                S_plot = permute(S_plot,[2,3,1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,2,(bindex))
                imshow(CS1n) 
            end              
            hold on   
            
            self.MMc = alignToCentralBin(self.MMc,mask,self.CalStru.binning_rc);

            for bindex = 1:2:dimmm(4)
                S_plot = self.MMc(1:3,:,:,bindex);
                S_plot = permute(S_plot,[2,3,1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,2,(bindex+1))
                imshow(CS1n) 
            end    
            hold off
            
        end

        function obtainSymmetrizationCorrection(self,imask)
            
            %imask = permute(imask,[16 1 2 4]);
            mask = repmat(self.nan_mask,[16,1,1,1]);
            mask = mask.*imask;
            
            for bindex = 1:self.bins
                MM = self.Mmueller(:,:,:,bindex);
                maskbin = mask(:,:,:,bindex);
                selectedMueller = reshape(MM(maskbin>0),16,[]);
                selectedMueller = selectedMueller./(selectedMueller(1,:));
                MMproj = mean(selectedMueller,2,'omitnan');
                MMprojsave(:,bindex) = MMproj;
                A = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0]/sqrt(2);
                F = conj(A'*reshape(MMproj,[4,4])*A);
                H = reshape(F([1,9,3,11,5,13,7,15,2,10,4,12,6,14,8,16]'),[4,4]); 
                [a,b] = eig(diag([-1,1,-1,1])*H*diag([-1,1,-1,1]));
                [~,pos] = min(abs(diag(b)));
                cc = a(:,pos);
                Jcorr = [cc(2),cc(4);cc(1),cc(3)];
                Jcorr = Jcorr/sqrt(det(Jcorr));
                [rr,dd] = JonesDecomp(Jcorr);% convert back to Stokes, ignoring diattanuation component

                % collect correction rotation vector for all spectral bins
                wcorr(:,bindex) = real(rr);
                ddcorr(:,bindex) = real(dd);

                % error computation
                jvec = [Jcorr(2,1),Jcorr(1,1),Jcorr(2,2),Jcorr(1,2)].';
                errFinal= real(jvec'*diag([-1,1,-1,1])*H*diag([-1,1,-1,1])*jvec/2); 
                % factor of 2 is as reference for entire energy; each pixel within ww is a normalized Jones matrix with sum(abs(J(:)).^2) = 2
                jvec = [0, 1, 1, 0]';
                errInit = real(jvec'*diag([-1,1,-1,1])*H*diag([-1,1,-1,1])*jvec/2);
                errFinal
                errInit
                self.Mcorr(:,bindex) = gather(self.Jones2Mueller(Jcorr(:)));
                self.Mcorr
                
            end

            %unwrap correction across spectral bins
            ss = cat(2,0,mod(cumsum(sqrt(sum(diff(wcorr,[],2).^2,1))>pi),2));
            if sum(ss)>self.bins/2
                ss = 1-ss;
            end
            retcorr = sqrt(sum(wcorr.^2,1));
            wcorr = wcorr - bsxfun(@times,bsxfun(@rdivide,wcorr,retcorr),ss*2*pi);

            % maintain a retardation <pi
            if mean(sqrt(sum(wcorr.^2,1)))>pi
                wcorr = wcorr - bsxfun(@rdivide,wcorr,sqrt(sum(wcorr.^2,1)))*2*pi;
            end
            
            figure(1)
            plot(wcorr')
            figure(2)
            plot(ddcorr')

        end
        
        function checkSymmetrizationCorrection(self)
            v = [6 10 14]; %[5 9 13];%
            vT = [6 7 8];%[2 3 4];%;
            
            bsize = round(self.bins/2);
            figure
            d_part = self.Mmueller(v,:,:,:);
            for bindex = 1:bsize
                S_plot = squeeze(d_part(:,:,:,2*bindex-1));
                S_plot = permute(S_plot,[2 3 1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,4,bindex*4-3)
                imshow(CS1n) 
            end   
            d_part = self.Mmueller(vT,:,:,:);
            d_part(3,:,:,:) = -d_part(3,:,:,:);
            for bindex = 1:bsize
                S_plot = squeeze(d_part(:,:,:,2*bindex-1));
                S_plot = permute(S_plot,[2 3 1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,4,bindex*4-2)
                imshow(CS1n) 
            end               
            
            for bindex = 1:self.bins
                corrMmueller(:,:,:,bindex) = MatrixMultiply(self.Mcorr(:,bindex),self.Mmueller(:,:,:,bindex));
            end
            
            d_part = corrMmueller(v,:,:,:);
            for bindex = 1:bsize
                S_plot = squeeze(d_part(:,:,:,2*bindex-1));
                S_plot = permute(S_plot,[2 3 1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,4,bindex*4-1)
                imshow(CS1n) 
            end
            d_part = corrMmueller(vT,:,:,:);
            d_part(3,:,:,:) = -d_part(3,:,:,:);
            for bindex = 1:bsize
                S_plot = squeeze(d_part(:,:,:,2*bindex-1));
                S_plot = permute(S_plot,[2 3 1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,4,bindex*4)
                imshow(CS1n) 
            end               
        end
        
        
        function applySymmetrizationCorrection(self)
            for bindex = 1:self.bins
                self.Mmueller(:,:,:,bindex) = MatrixMultiply(self.Mcorr(:,bindex),self.Mmueller(:,:,:,bindex));
            end
            MmuellerT = transpose4x4(self.Mmueller);
            MmuellerT([4 8 12 13 14 15],:,:,:) = - MmuellerT([4 8 12 13 14 15],:,:,:);
            self.Mmueller = (MmuellerT+self.Mmueller)/2;

        end
        
        
        function obtainPMDCorrection(self,imask)
            %imask = permute(imask,[3 1 2 4]);
            mask = repmat(self.nan_mask,[16,1,1,1]);
            mask = mask.*imask;
            

            v0 = [0 0 0 0 0 0];            

            cbi = round(self.bins/2);

            v(:,cbi) = [0 0 0 0 0 0];
            
            for indw = cbi+1:self.bins
                indw
                v(:,indw) = self.search_for_align_vectors(v0,self.Mmueller(:,:,:,indw)./self.Mmueller(1,:,:,indw),self.Mmueller(:,:,:,cbi)./self.Mmueller(1,:,:,cbi),imask);
                v0 = v(:,indw);
            end
            
            v0 = [0 0 0 0 0 0]; 
            
            for indw = cbi-1:-1:1 
                v(:,indw) = self.search_for_align_vectors(v0,self.Mmueller(:,:,:,indw)./self.Mmueller(1,:,:,indw),self.Mmueller(:,:,:,cbi)./self.Mmueller(1,:,:,cbi),imask);
                v0 = v(:,indw);
            end

            
            for indw =1:self.bins
                omeg = [0 v(4,indw) v(5,indw) v(6,indw)
                        v(4,indw) 0 -v(3,indw) v(2,indw)
                       v(5,indw)  v(3,indw) 0 -v(1,indw)
                        v(6,indw) -v(2,indw) v(1,indw) 0];
                rotation = expm(omeg);
                D = [1 0 0 0
                     0 1 0 0
                     0 0 1 0
                     0 0 0 -1];  
                DrotationTD =D*rotation.'*D;

                self.PMDCorr(:,indw) = rotation(:);
                self.DPMDCorrTD(:,indw) = DrotationTD(:);
                
            end
                
            figure
            plot(v')

        
        end
        
        
        function v = search_for_align_vectors(self,v0,M1,M2,mask)
            
            f = @(x)objective_fun_for_align_vector(x,M1,M2,mask,self.CalStru.GPU);
            %options = optimset('Display','iter','MaxFunEvals',3000,'PlotFcns',@optimplotfval);
            %v = fminsearch(f,v0,options);
            v = fminsearch(f,v0);

        end        
        
        
        function flip_Mueller(self)
            self.Mmueller = flip(self.Mmueller,2);
        end
        

        function checkPMDCorrection(self)
            v = [5 9 13];%[6 10 14]; %
            vT = [10 11 12];%;[2 3 4];%
            
            bsize = round(self.bins/2);
            figure
            d_part = self.Mmueller(v,:,:,:);
            for bindex = 1:bsize
                S_plot = squeeze(d_part(:,:,:,2*bindex-1));
                S_plot = permute(S_plot,[2 3 1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,4,bindex*4-3)
                imshow(CS1n) 
            end   
            d_part = self.Mmueller(vT,:,:,:);
            d_part(3,:,:,:) = -d_part(3,:,:,:);
            for bindex = 1:bsize
                S_plot = squeeze(d_part(:,:,:,2*bindex-1));
                S_plot = permute(S_plot,[2 3 1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,4,bindex*4-2)
                imshow(CS1n) 
            end               
            
            for bindex = 1:self.bins
                corrMmueller(:,:,:,bindex) = MatrixMultiply(MatrixMultiply(self.PMDCorr(:,bindex),self.Mmueller(:,:,:,bindex)),self.DPMDCorrTD(:,bindex));
            end
            
            d_part = corrMmueller(v,:,:,:);
            for bindex = 1:bsize
                S_plot = squeeze(d_part(:,:,:,2*bindex-1));
                S_plot = permute(S_plot,[2 3 1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,4,bindex*4-1)
                imshow(CS1n) 
            end
            d_part = corrMmueller(vT,:,:,:);
            d_part(3,:,:,:) = -d_part(3,:,:,:);
            for bindex = 1:bsize
                S_plot = squeeze(d_part(:,:,:,2*bindex-1));
                S_plot = permute(S_plot,[2 3 1]);
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,4,bindex*4)
                imshow(CS1n) 
            end   
            
%             MVar = var(corrMmueller,0,4,'omitnan');
%             figure
%             ivar = squeeze(mean(MVar,1,'omitnan'));
%             imagesc(ivar);
            
            
        end
        
        function applyPMDCorrection(self)
            for bindex = 1:self.bins
                self.Mmueller(:,:,:,bindex) = MatrixMultiply(MatrixMultiply(self.PMDCorr(:,bindex),self.Mmueller(:,:,:,bindex)),self.DPMDCorrTD(:,bindex));
            end  
            self.Mmueller = mean(self.Mmueller,4,'omitnan');
            %self.Mmueller = self.Mmueller(:,:,:,5);
            
        end

         
        
        function ret = get_local_retardance(self)
            
            mR = self.polarDecomposition(self.Mmueller);
            mRdz = circshift(mR,-1,2);
            mLocal = MatrixMultiply(mR,matrix3inv(mRdz));     
            r = decomposeRot(mLocal);   
            self.cumulative_r = r;
            ret= squeeze(sqrt(sum(r.^2,1))/2/self.pixel_depth/pi*180);
            
        end
        
        function [mR, mD] = polarDecomposition(self,M)
            dim = size(M);
            if dim(1) ~= 16
                error('fisrt dimnetion must be 16');
            end
            M = reshape(M,16,[]); 
            m00 = M(1,:);
            M = M./m00;
            D = M([5 6 13],:);
            cos2k = dot(D,D,1);
            %self.nan_mask(cos2k>1) = 0;
            cos2k(cos2k>1) = 1-10e-5;
            sink = sqrt(1-cos2k);         
            mD = ((1-sink)./cos2k).*self.kron_vector_expand(D);
            mD(1,:) = mD(1,:)+sink(1,:);mD(5,:) = mD(5,:)+sink(1,:);mD(9,:) = mD(9,:)+sink(1,:);
            
            
            [iMmDI,~] = matrix3inv(mD);
            mR = MatrixMultiply(M([6 7 8 10 11 12 14 15 16],:),iMmDI);
            mR = euclideanRotation(mR);
            dim(1) = 9;
            mR = reshape(mR,dim); 
            mD = reshape(mD,dim); 
        end
        
        function  remove_linear_birefringence_from_surface(self)
            
            dim = size(self.Mmueller);
            
            % normalize the Mueller matrix
            surface_Mueller = squeeze(self.Mmueller(:,1,:)./self.Mmueller(1,1,:));
            surface_Mueller(isnan(surface_Mueller)) = 1;
            
            % dop masking
            dopmask = self.get_depolarization_index()>0.90;
            dopmask = dopmask(1,:);  

            surface_Mueller = self.Mulluer_fit_filter(surface_Mueller,dopmask);
                            %surface_Mueller = medfilt1(gather(surface_Mueller),30,[],2);

                            %plot(surface_Mueller')

            
            % simple logm on all the surface Mueller matrix to get
            % retardance and diattenuation vectors
            [ret,dia] = self.slow_Mueller2vectors(surface_Mueller);
            
            
            % I tried to unwrap, but not working
                                %%%%%%%%%%%%%%%%%
                    %             ret = ret;
                    %             thit = sqrt(dot(ret,ret,1));
                    %             thit_f = medfilt1(thit,5,[],2);
                    %             ret = (thit_f./thit).*ret;
                    %            ret = self.layer_rotV_align(ret,dopmask);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % reconstruct the Mueller matrix
            
            inv_half_surface_Mueller = self.slow_vectors2Mueller(-ret./2,-dia./2);
            
            if isfield(self.CalStru,'GPU')
                if self.CalStru.GPU == 1
                    inv_half_surface_Mueller= gpuArray(single(inv_half_surface_Mueller));
                else
                end
            else
            end
            
            % apply to image Mueller matrix
            inv_half_surface_Mueller = repmat(permute(inv_half_surface_Mueller,[1 3 2]),1,dim(2),1);
            self.Mmueller = MatrixMultiply(inv_half_surface_Mueller,MatrixMultiply(self.Mmueller,inv_half_surface_Mueller));
            
        end
        
        function remove_circular_diatenuation(self)
            
            dim = size(self.Mmueller);
            mask = self.get_depolarization_index()>0.6;
            
            normalize_Mueller = (self.Mmueller(:,:,:)./self.Mmueller(1,:,:)).*permute(mask,[3 1 2]);%
            
            normalize_Mueller(isnan(normalize_Mueller)) = 0;
            frame_M = (sum(normalize_Mueller,2));%normalize_Mueller;%
            mask = (sum(mask,2));

            tanh2z = (-(2.*frame_M(13,:,:))./(frame_M(1,:,:)-frame_M(16,:,:)));
            
            tanh2z(tanh2z>1) = 0.99;
            tanh2z(tanh2z<-1) = -0.99;
            
            z = atanh(tanh2z)/2;
            z = medfilt1(gather(z),30,[],3);
            z(isnan(z)) = 0;

            [MDC,invMDC] = self.circularDiattenuator(z);

            MDC = repmat(MDC,1,dim(2),1);
            invMDC = repmat(invMDC,1,dim(2),1);
            
            self.Mmueller = MatrixMultiply(invMDC,MatrixMultiply(self.Mmueller,MDC));
            
        end
        

        function [MDC,invMDC] = circularDiattenuator(self,z)
            
            dim = size(z);
            if dim(1) ~= 1
                error('should be 1 in first dimension')
            end
            dim(1) = 16;
            MDC = zeros(dim);
            if isfield(self.CalStru,'GPU')
                if self.CalStru.GPU == 1
                    MDC = gpuArray(single(MDC));
                else
                end
            else
            end
            MDC(1,:) = cosh(z(1,:));
            MDC(4,:) = sinh(z(1,:));
            MDC(6,:) = 1;
            MDC(11,:)= 1;
            MDC(13,:)=sinh(z(1,:));
            MDC(16,:)=cosh(z(1,:));
            
            invMDC = MDC;

            invMDC(4,:)  = -MDC(4,:);
            invMDC(13,:) = -MDC(13,:);            
        end
        
        
        function SNR = SNRmap(self,I)
            
            %I = reshape(I,self.dim(1)*self.dim(2),[]);
            ground_noise = quantile(I(:),0.15);
            SNR = I./ground_noise;
            %SNR = reshape(SNR,self.dim(1:end-1));
  
        end
        
        function ma = getMuellerMask(self)
            ma = gather(squeeze(sum(self.nan_mask,4)));
        end
        
        function D = getCumDiattenuation(self)
            M = self.Mmueller./self.Mmueller(1,:,:);
            D = squeeze(M(5,:,:).^2+M(9,:,:).^2+M(13,:,:).^2);
            D = D.*(self.get_depolarization_index()>0.6);
            
        end
            
         
        
        function [m1,m2] = get_retina_n_sclera_mask(self,I,ret)
            
            if self.CalStru.GPU == 1
                I = gpuArray(single(I));
                ret = gpuArray(single(ret));
            else
            end 

            Im = (imgaussfilt(I,3));
            %retm = (imgaussfilt(ret,3));
            
            retm = ret < self.ret_level;

            snrI = self.SNRmap(Im);
%             
%             %%%%%%%%%%
            Im = snrI>1;
            dop = self.dopmap()>0.5;
            
            mr = dop&retm&Im;
            %%%%%%%%%%%%%%%%%%%%%%%%
            se1 = strel('disk',1);
            mr = imopen(mr,se1);
            self.retina_mask = mr;
            
            m1 = self.retina_mask;
            se1 = strel('disk',10);
            mr = imopen(mr,se1);
            se2 = strel('rectangle',[15,15]);
            m2 = imclose(mr,se2); 
            
%             scleramask = (ret > 0.25)&(m1);
%             scleramask = imopen(scleramask,se1);
%             scleramask = imclose(scleramask,se2);   
%             m2 = scleramask;
            %m2 =retm&m1;
            
        end
        
        function [CalStru] = getCalStru(self)
            self.CalStru.binning_rc.Mcorr = self.Mcorr;
            self.CalStru.binning_rc.PMDCorr = self.PMDCorr;
            self.CalStru.binning_rc.DPMDCorrTD = self.DPMDCorrTD;
            self.CalStru.binning_rc.order = self.order;
            
            CalStru = self.CalStru;
        end
        
       
        
        function mask = snrmask(self,level)
            SNR = mean((self.SNRmap(self.I1)+self.SNRmap(self.I2)+self.SNRmap(self.I3))/3,1);
            SNR = 10*log10(SNR);
            mask = SNR>level;
        end       


        function line = all_depth_ret(self,fret,oa)
            
            fret(fret>self.ret_level) = 0;
            rec = fret.*exp(1i.*oa); %why we have NANs
            rec(isnan(rec)) = 0;
            line = abs(mean(rec));
            
        end

        
        
        function [retFinal,phi] = local_optic_axis(self)
            %%
            % generate mask to exclude points with low DOP from correcting the local
            % retardation
            %mmask = self.get_depolarization_index()>0.6;

            
            dimm = size(self.Mmueller);     
            
            
            mR = self.polarDecomposition(self.Mmueller);
            mRdz = circshift(mR,-1,2);
            mLocal = MatrixMultiply(mR,matrix3inv(mRdz));     
            self.cumulative_r = decomposeRot(mR);   
            r = decomposeRot(mLocal);   
            
            
            MMM = mR;%self.Mmueller([6,7,8,10,11,12,14,15,16],:,:);
            
            %MMM = euclideanRotation(MMM);
            
            
            nanmask = squeeze(sum(isnan(MMM),1));
            
            if isfield(self.CalStru,'GPU')              
                if self.CalStru.GPU == 1                
                    w = gpuArray(single(zeros([3,dimm(2),dimm(3)])));
                else
                    w = zeros([3,dimm(2),dimm(3)]);
                end
            else

            end


            mmask = ones(dimm(2),dimm(3));

            N = repmat([1;0;0;0;1;0;0;0;1],[1,1,dimm(3)]);
            if isfield(self.CalStru,'GPU')              
                if self.CalStru.GPU == 1                
                    N = gpuArray(single(N));
                else
                end
            else

            end            
            
%             N = permute(makeRot(cat(1,zeros(2,numel(WW)),-WW)),[1,3,2]);
      
            
            for indz = 2:150 %self.dim(1)%refInds(3):roiz(end)

                % D*N(indz)*D*Mtot(indz+1)*N'(indz)
                nloc = MatrixMultiply(bsxfun(@times,N,[1;1;-1;1;1;-1;-1;-1;1]),MatrixMultiply(MMM(:,indz,:,:),N([1,4,7,2,5,8,3,6,9],:,:,:)));
                Wloc = decomposeRot(nloc);   
                w(:,indz,:) = squeeze(Wloc);  
                nlocroot = makeRot(Wloc/2);

                Nnext = MatrixMultiply(nlocroot,N);
                %N(:,:,:) = Nnext(:,:,:);
                N(:,:,mmask(indz,:)>0) = Nnext(:,:,mmask(indz,:)>0);

            end
            
            %Omega = bsxfun(@times,permute(pa,[2,3,1]),ret);
            %Omegaf = imfilter(permute(w,[2,3,1]),ones(self.fwz,1)/self.fwz);
            %Omegaf = imfilter(permute(w,[2,3,1]),self.zfilt_kernal);
            Omegaf = permute(w,[2,3,1]);
            retFinal= squeeze(sqrt(sum(r.^2,1))/2/self.pixel_depth/pi*180);
            %retFinal = sqrt(sum(Omegaf.^2,3))/2/self.pixel_depth/pi*180;
            phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
            %phi = bsxfun(@rdivide,Omegaf,sqrt(sum(Omegaf.^2,3)));
            
        end
        
        
        function [retarder,diattenuator] = slow_Mueller2vectors(self,M)
             warning ('off','all');
             dm = real(self.my_logm(M));
             warning ('on','all');
             diattenuator = (cat(1,dm(5,:)+dm(2,:),dm(9,:)+dm(3,:),dm(13,:)+dm(4,:)))/2;
             retarder = -(cat(1,dm(12,:)-dm(15,:),dm(14,:)-dm(8,:),dm(7,:)-dm(10,:)))/2;
        end
        
        function M = slow_vectors2Mueller(self,retarder,diattenuator)
            dim = size(retarder);
            dim(1) = 16;
            M = zeros(dim);
            M(2,:,:,:) = diattenuator(1,:,:,:);
            M(3,:,:,:) = diattenuator(2,:,:,:);
            M(4,:,:,:) = diattenuator(3,:,:,:);
            M(5,:,:,:) = diattenuator(1,:,:,:);
            M(9,:,:,:) = diattenuator(2,:,:,:);
            M(13,:,:,:) = diattenuator(3,:,:,:);       
            M(7,:,:,:) = -retarder(3,:,:,:);  
            M(8,:,:,:) = retarder(2,:,:,:); 
            M(10,:,:,:) = retarder(3,:,:,:); 
            M(12,:,:,:) = -retarder(1,:,:,:); 
            M(14,:,:,:) = -retarder(2,:,:,:); 
            M(15,:,:,:) = retarder(1,:,:,:);
            M = real(self.my_expm(M));
        end        
        
        
        function Mm = Mulluer_fit_filter(self,Mm,surf)
            if sum(surf(:))<5
                return;
            end
            [Sc, L] = size(Mm);
            xgrid = (1:L);
            for sindex = 1:Sc
                Ma = Mm(sindex,:);
                [p,~,mu] = polyfit(xgrid(surf),Ma(surf),5);
                Mm(sindex,:) = polyval(p,1:L,[],mu);
            end
        end
        

        function WW = layer_rotV_align(self,WFit,surf)

            ss = repmat(cat(2,0,mod(cumsum(sqrt(sum(diff(WFit,[],2).^2,1))>pi,2),2)),[3,1])>0;
            delta = 2*pi*bsxfun(@rdivide,WFit,sqrt(sum(WFit.^2,1)));
            WFitc = WFit;
            WFitc(ss) = WFit(ss) - delta(ss);
         
            n = -1:1;
            for ind = 1:3
                WW = WFitc + 2*pi*n(ind)*bsxfun(@rdivide,WFitc,sqrt(sum(WFitc.^2,1)));

                %WW(~surf) = 0;
                
                meanRet(ind) = mean(mean(sqrt(sum((WW-self.last_surface_WW).^2,1))));
            end
            [~,mp] = min(meanRet);
            WW = WFitc + 2*pi*n(mp)*bsxfun(@rdivide,WFitc,sqrt(sum(WFitc.^2,1)));               
            
            self.last_surface_WW = WW;

        end
        
        
        function test_polarization(self,S)
            figure
            
            dim = size(S);
            bsize = round(dim(3)/2);
            for bindex = 1:bsize
                S_plot = squeeze(S(:,:,2*bindex-1,:));
                Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
                CS1n = stokes2color(Sj1n);
                subplot(bsize,1,bindex)
                imshow(CS1n) 
            end        
        end
        
        function dopi = get_depolarization_index(self)
            dimm = size(self.Mmueller);
            M = self.Mmueller;
            
            if dimm(1) ~= 16
                error('fisrt dimnetion must be 16');
            end
            M = reshape(M,16,[]);
            M = M./M(1,:);
            T1 = M(1,:).^2+M(2,:).^2+M(3,:).^2+M(4,:).^2;
            T2 = M(5,:).^2+M(6,:).^2+M(7,:).^2+M(8,:).^2;
            T3 = M(9,:).^2+M(10,:).^2+M(11,:).^2+M(12,:).^2;
            T4 = M(13,:).^2+M(14,:).^2+M(15,:).^2+M(16,:).^2;
            dopi = sqrt((T1+T2+T3+T4-1)/3);
            %dopi(dopi>1) = 0;
            dimm(1) = 1;
            dopi(isnan(dopi)) = 0;
            dopi = squeeze(reshape(dopi,dimm));
            
        end
        
        
        function retina_mask = get_retina_mask(self,oaxu,dop,Isnr,ret)
            
            %dop = self.get_depolarization_index();
            retina_mask1 = (dop>0.2)& oaxu;% 
            retina_mask = (retina_mask1&(Isnr>1.5)&(ret<6));
            SE = strel('disk',3);
            mask = imopen(retina_mask,SE);
            
            retina_mask(150:end,:)=0;
        end

        function array = zfiltering(self,array)
            array = permute(array,[2 3 4 1]);
            dim = size(array);
            array = reshape(array,dim(1),dim(2),[]);
            array = imfilter(array,self.zfilt_kernal,'same');
            array = reshape(array,dim);
            array = permute(array,[4 1 2 3]);

        end

        
        function lm = my_logm(self,M)
            M = gather(M);
            %dimj = size(M);
            cellM = num2cell(M,1);
            
            m = cellfun(@InCell_logm,cellM,'UniformOutput',false);
            lm = cell2mat(m);
        end
        
        function em = my_expm(self,M)
            M = gather(M);
            %dimj = size(M);
            cellM = num2cell(M,1);
            
            m = cellfun(@InCell_expm,cellM,'UniformOutput',false);
            em = cell2mat(m);
        end        
        

        function M = Jones2Mueller(self,J)
            
            dimm = size(J);
            J = reshape(J,4,[]);
            
            dimm(1) = 16;
            M = zeros(dimm);
            
            if isfield(self.CalStru,'GPU')
                if self.CalStru.GPU == 1
                    M = gpuArray(M);
                else
                end
            else
            end

            M(1,:) = conj(J(1,:)).*J(1,:)+conj(J(3,:)).*J(3,:)+conj(J(2,:)).*J(2,:)+conj(J(4,:)).*J(4,:);
            M(5,:) = conj(J(1,:)).*J(1,:)+conj(J(2,:)).*J(2,:)-conj(J(3,:)).*J(3,:)-conj(J(4,:)).*J(4,:);
            M(9,:) = conj(J(1,:)).*J(3,:)+conj(J(2,:)).*J(4,:)+conj(J(3,:)).*J(1,:)+conj(J(4,:)).*J(2,:);
            M(13,:) = 1i.*(conj(J(1,:)).*J(3,:)+conj(J(2,:)).*J(4,:)-conj(J(3,:)).*J(1,:)-conj(J(4,:)).*J(2,:));
            M(2,:) = conj(J(1,:)).*J(1,:)+conj(J(3,:)).*J(3,:)-conj(J(2,:)).*J(2,:)-conj(J(4,:)).*J(4,:);
            M(6,:) = conj(J(1,:)).*J(1,:)+conj(J(4,:)).*J(4,:)-conj(J(2,:)).*J(2,:)-conj(J(3,:)).*J(3,:);
            M(10,:) = conj(J(3,:)).*J(1,:)+conj(J(1,:)).*J(3,:)-conj(J(4,:)).*J(2,:)-conj(J(2,:)).*J(4,:);
            M(14,:) = 1i.*(conj(J(1,:)).*J(3,:)+conj(J(4,:)).*J(2,:)-conj(J(2,:)).*J(4,:)-conj(J(3,:)).*J(1,:));
            M(3,:) = conj(J(1,:)).*J(2,:)+conj(J(2,:)).*J(1,:)+conj(J(3,:)).*J(4,:)+conj(J(4,:)).*J(3,:);
            M(7,:) = conj(J(1,:)).*J(2,:)+conj(J(2,:)).*J(1,:)-conj(J(3,:)).*J(4,:)-conj(J(4,:)).*J(3,:);
            M(11,:) = conj(J(1,:)).*J(4,:)+conj(J(2,:)).*J(3,:)+conj(J(3,:)).*J(2,:)+conj(J(4,:)).*J(1,:);
            M(15,:) = 1i.*(conj(J(1,:)).*J(4,:)+conj(J(2,:)).*J(3,:)-conj(J(3,:)).*J(2,:)-conj(J(4,:)).*J(1,:));
            M(4,:) = 1i.*(conj(J(2,:)).*J(1,:)+conj(J(4,:)).*J(3,:)-conj(J(1,:)).*J(2,:)-conj(J(3,:)).*J(4,:));
            M(8,:) = 1i.*(conj(J(2,:)).*J(1,:)+conj(J(3,:)).*J(4,:)-conj(J(1,:)).*J(2,:)-conj(J(4,:)).*J(3,:));
            M(12,:) = 1i.*(conj(J(2,:)).*J(3,:)+conj(J(4,:)).*J(1,:)-conj(J(1,:)).*J(4,:)-conj(J(3,:)).*J(2,:));   
            M(16,:) = conj(J(4,:)).*J(1,:)+conj(J(1,:)).*J(4,:)-conj(J(3,:)).*J(2,:)-conj(J(2,:)).*J(3,:);
            M = real(M./2);
            M = reshape(M,dimm);
            
        end
        
        
        function checkStokes(self,S1,S2,S3)

            % compute corresponding intensity signal
            I1 = sqrt(sum(S1.^2,4));
            I2 = sqrt(sum(S2.^2,4));
            I3 = sqrt(sum(S3.^2,4));

            S1n = S1./I1;
            S2n = S2./I2;
            S3n = S3./I3;


            dotProd12 = (sum(S1n.*S2n,4));
            dotProd13 = (sum(S1n.*S3n,4));
            dotProd23 = (sum(S2n.*S3n,4));

            % simple intensity mask to identify areas to analyze the dot product
            mask1 = 10*log10(I1)>quantile(10*log10(I1(:)),0.6);
            mask2 = 10*log10(I2)>quantile(10*log10(I2(:)),0.6);
            mask3 = 10*log10(I3)>quantile(10*log10(I3(:)),0.6);
            
            mask =(mask1&mask2&mask3);
            
            
            for ind = 1:9
                angle3 = acosd(dotProd12(:,:,ind));
                angle2 = acosd(dotProd13(:,:,ind));
                angle1 = acosd(dotProd23(:,:,ind));
                a1(ind) = median(angle1(mask(:,:,ind)),'all'); 
                a2(ind) = median(angle2(mask(:,:,ind)),'all'); 
                a3(ind) = median(angle3(mask(:,:,ind)),'all'); 
            end
            figure
            plot(a1)
            hold on 
            plot(a2)
            plot(a3)
            hold off
            
            
            
            
            bins = linspace(-1,1,100);
            for ind = 1:9
                loc = dotProd12(:,:,ind);
                [hh,xh] = ksdensity(loc(mask(:,:,ind)),bins);
                HH12(:,ind) = hh;

                loc = dotProd13(:,:,ind);
                [hh,xh] = ksdensity(loc(mask(:,:,ind)),bins);
                HH13(:,ind) = hh;

                loc = dotProd23(:,:,ind);
                [hh,xh] = ksdensity(loc(mask(:,:,ind)),bins);
                HH23(:,ind) = hh;
            end



            %%
            % visualize a given spectral bing
            indw = 5;

            figure(1)
            clf
            subplot(2,4,1)
            %imagesc(10*log10(I))
            title('Intensity')

            subplot(2,4,5)
            imagesc(mask(:,:,indw))
            title(sprintf('Intensity bin %d',indw))

            subplot(2,4,[2,3,4])
            imagesc(cat(2,S1n(:,:,indw,1),S1n(:,:,indw,2),S1n(:,:,indw,3)))
            title(sprintf('QUV bin %d',indw))

            subplot(2,4,6)
            imagesc(dotProd12(:,:,indw))
            title(sprintf('Dot product 12, bin %d',indw))

            subplot(2,4,7)
            imagesc(dotProd13(:,:,indw))
            title(sprintf('Dot product 13, bin %d',indw))

            subplot(2,4,8)
            imagesc(dotProd23(:,:,indw))
            title(sprintf('Dot product 23, bin %d',indw))


            figure(2)
            clf
            subplot(1,3,1)
            imagesc(HH12)
            title('Histograms dotproduct 12')

            subplot(1,3,2)
            imagesc(HH13)
            title('Histograms dotproduct 13')

            subplot(1,3,3)
            imagesc(HH23)
            title('Histograms dotproduct 23')


            
            
        end
        
        function [cha1,cha2] = solve_cha(self,m1,m2,m3)
            
            
            is11 = permute(self.isis(1,:),[1 3 4 2]);
            is21 = permute(self.isis(2,:),[1 3 4 2]);
            is31 = permute(self.isis(3,:),[1 3 4 2]);
            is12 = permute(self.isis(4,:),[1 3 4 2]);
            is22 = permute(self.isis(5,:),[1 3 4 2]);
            is32 = permute(self.isis(6,:),[1 3 4 2]);
            is13 = permute(self.isis(7,:),[1 3 4 2]);
            is23 = permute(self.isis(8,:),[1 3 4 2]);
            is33 = permute(self.isis(9,:),[1 3 4 2]);
            
            cha2 = (2.*is11.*m1 + is12.*m1 + is12.*m2 + is13.*m1 + is13.*m3 + is21.*m1 + is21.*m2 + 2.*is22.*m2 + is23.*m2 + is23.*m3 + is31.*m1 + is31.*m3 + is32.*m2 + is32.*m3 + 2.*is33.*m3 + abs(4.*is11.*m1.^2 - 4.*is12 - 4.*is13 - 4.*is21 - 4.*is22 - 4.*is23 - 4.*is31 - 4.*is32 - 4.*is33 - 4.*is11 + 4.*is22.*m2.^2 + 4.*is33.*m3.^2 + is12.^2.*m1.^2 + is12.^2.*m2.^2 + is13.^2.*m1.^2 + is13.^2.*m3.^2 + is21.^2.*m1.^2 + is21.^2.*m2.^2 + is23.^2.*m2.^2 + is23.^2.*m3.^2 + is31.^2.*m1.^2 + is31.^2.*m3.^2 + is32.^2.*m2.^2 + is32.^2.*m3.^2 - 2.*is12.^2.*m1.*m2 - 2.*is13.^2.*m1.*m3 - 2.*is21.^2.*m1.*m2 - 2.*is23.^2.*m2.*m3 - 2.*is31.^2.*m1.*m3 - 2.*is32.^2.*m2.*m3 + 4.*is12.*m1.*m2 + 4.*is13.*m1.*m3 + 4.*is21.*m1.*m2 + 4.*is23.*m2.*m3 + 4.*is31.*m1.*m3 + 4.*is32.*m2.*m3 + 2.*is12.*is13.*m1.^2 - 4.*is11.*is22.*m1.^2 + 2.*is12.*is21.*m1.^2 - 4.*is11.*is22.*m2.^2 - 4.*is11.*is23.*m1.^2 + 2.*is12.*is21.*m2.^2 + 2.*is13.*is21.*m1.^2 + 2.*is12.*is23.*m2.^2 - 4.*is13.*is22.*m2.^2 + 2.*is13.*is23.*m3.^2 - 4.*is11.*is32.*m1.^2 + 2.*is12.*is31.*m1.^2 - 4.*is11.*is33.*m1.^2 + 2.*is13.*is31.*m1.^2 + 2.*is12.*is32.*m2.^2 + 2.*is21.*is23.*m2.^2 - 4.*is11.*is33.*m3.^2 + 2.*is13.*is31.*m3.^2 - 4.*is12.*is33.*m3.^2 + 2.*is13.*is32.*m3.^2 + 2.*is21.*is31.*m1.^2 + 2.*is21.*is32.*m2.^2 - 4.*is22.*is31.*m2.^2 - 4.*is21.*is33.*m3.^2 - 4.*is22.*is33.*m2.^2 + 2.*is23.*is31.*m3.^2 + 2.*is23.*is32.*m2.^2 - 4.*is22.*is33.*m3.^2 + 2.*is23.*is32.*m3.^2 + 2.*is31.*is32.*m3.^2 - 2.*is12.*is13.*m1.*m2 - 2.*is12.*is13.*m1.*m3 + 2.*is12.*is13.*m2.*m3 + 8.*is11.*is22.*m1.*m2 - 4.*is12.*is21.*m1.*m2 + 4.*is11.*is23.*m1.*m2 - 2.*is13.*is21.*m1.*m2 + 4.*is11.*is23.*m1.*m3 - 2.*is12.*is23.*m1.*m2 - 2.*is13.*is21.*m1.*m3 + 4.*is13.*is22.*m1.*m2 - 4.*is11.*is23.*m2.*m3 + 2.*is12.*is23.*m1.*m3 + 2.*is13.*is21.*m2.*m3 - 4.*is13.*is22.*m1.*m3 + 2.*is13.*is23.*m1.*m2 - 2.*is12.*is23.*m2.*m3 + 4.*is13.*is22.*m2.*m3 - 2.*is13.*is23.*m1.*m3 - 2.*is13.*is23.*m2.*m3 + 4.*is11.*is32.*m1.*m2 - 2.*is12.*is31.*m1.*m2 + 4.*is11.*is32.*m1.*m3 - 2.*is12.*is31.*m1.*m3 - 2.*is12.*is32.*m1.*m2 - 2.*is21.*is23.*m1.*m2 - 4.*is11.*is32.*m2.*m3 + 8.*is11.*is33.*m1.*m3 + 2.*is12.*is31.*m2.*m3 + 2.*is12.*is32.*m1.*m3 - 4.*is12.*is33.*m1.*m2 - 4.*is13.*is31.*m1.*m3 + 2.*is13.*is32.*m1.*m2 + 2.*is21.*is23.*m1.*m3 - 2.*is12.*is32.*m2.*m3 + 4.*is12.*is33.*m1.*m3 - 2.*is13.*is32.*m1.*m3 - 2.*is21.*is23.*m2.*m3 + 4.*is12.*is33.*m2.*m3 - 2.*is13.*is32.*m2.*m3 - 2.*is21.*is31.*m1.*m2 - 2.*is21.*is31.*m1.*m3 - 2.*is21.*is32.*m1.*m2 + 4.*is22.*is31.*m1.*m2 + 2.*is21.*is31.*m2.*m3 + 2.*is21.*is32.*m1.*m3 - 4.*is21.*is33.*m1.*m2 - 4.*is22.*is31.*m1.*m3 + 2.*is23.*is31.*m1.*m2 - 2.*is21.*is32.*m2.*m3 + 4.*is21.*is33.*m1.*m3 + 4.*is22.*is31.*m2.*m3 - 2.*is23.*is31.*m1.*m3 + 4.*is21.*is33.*m2.*m3 - 2.*is23.*is31.*m2.*m3 + 8.*is22.*is33.*m2.*m3 - 4.*is23.*is32.*m2.*m3 + 2.*is31.*is32.*m1.*m2 - 2.*is31.*is32.*m1.*m3 - 2.*is31.*is32.*m2.*m3 + 4).^(1/2))./(2.*(is11 + is12 + is13 + is21 + is22 + is23 + is31 + is32 + is33 - 1));
            cha1 = (2.*is11.*m1 + is12.*m1 + is12.*m2 + is13.*m1 + is13.*m3 + is21.*m1 + is21.*m2 + 2.*is22.*m2 + is23.*m2 + is23.*m3 + is31.*m1 + is31.*m3 + is32.*m2 + is32.*m3 + 2.*is33.*m3 - abs(4.*is11.*m1.^2 - 4.*is12 - 4.*is13 - 4.*is21 - 4.*is22 - 4.*is23 - 4.*is31 - 4.*is32 - 4.*is33 - 4.*is11 + 4.*is22.*m2.^2 + 4.*is33.*m3.^2 + is12.^2.*m1.^2 + is12.^2.*m2.^2 + is13.^2.*m1.^2 + is13.^2.*m3.^2 + is21.^2.*m1.^2 + is21.^2.*m2.^2 + is23.^2.*m2.^2 + is23.^2.*m3.^2 + is31.^2.*m1.^2 + is31.^2.*m3.^2 + is32.^2.*m2.^2 + is32.^2.*m3.^2 - 2.*is12.^2.*m1.*m2 - 2.*is13.^2.*m1.*m3 - 2.*is21.^2.*m1.*m2 - 2.*is23.^2.*m2.*m3 - 2.*is31.^2.*m1.*m3 - 2.*is32.^2.*m2.*m3 + 4.*is12.*m1.*m2 + 4.*is13.*m1.*m3 + 4.*is21.*m1.*m2 + 4.*is23.*m2.*m3 + 4.*is31.*m1.*m3 + 4.*is32.*m2.*m3 + 2.*is12.*is13.*m1.^2 - 4.*is11.*is22.*m1.^2 + 2.*is12.*is21.*m1.^2 - 4.*is11.*is22.*m2.^2 - 4.*is11.*is23.*m1.^2 + 2.*is12.*is21.*m2.^2 + 2.*is13.*is21.*m1.^2 + 2.*is12.*is23.*m2.^2 - 4.*is13.*is22.*m2.^2 + 2.*is13.*is23.*m3.^2 - 4.*is11.*is32.*m1.^2 + 2.*is12.*is31.*m1.^2 - 4.*is11.*is33.*m1.^2 + 2.*is13.*is31.*m1.^2 + 2.*is12.*is32.*m2.^2 + 2.*is21.*is23.*m2.^2 - 4.*is11.*is33.*m3.^2 + 2.*is13.*is31.*m3.^2 - 4.*is12.*is33.*m3.^2 + 2.*is13.*is32.*m3.^2 + 2.*is21.*is31.*m1.^2 + 2.*is21.*is32.*m2.^2 - 4.*is22.*is31.*m2.^2 - 4.*is21.*is33.*m3.^2 - 4.*is22.*is33.*m2.^2 + 2.*is23.*is31.*m3.^2 + 2.*is23.*is32.*m2.^2 - 4.*is22.*is33.*m3.^2 + 2.*is23.*is32.*m3.^2 + 2.*is31.*is32.*m3.^2 - 2.*is12.*is13.*m1.*m2 - 2.*is12.*is13.*m1.*m3 + 2.*is12.*is13.*m2.*m3 + 8.*is11.*is22.*m1.*m2 - 4.*is12.*is21.*m1.*m2 + 4.*is11.*is23.*m1.*m2 - 2.*is13.*is21.*m1.*m2 + 4.*is11.*is23.*m1.*m3 - 2.*is12.*is23.*m1.*m2 - 2.*is13.*is21.*m1.*m3 + 4.*is13.*is22.*m1.*m2 - 4.*is11.*is23.*m2.*m3 + 2.*is12.*is23.*m1.*m3 + 2.*is13.*is21.*m2.*m3 - 4.*is13.*is22.*m1.*m3 + 2.*is13.*is23.*m1.*m2 - 2.*is12.*is23.*m2.*m3 + 4.*is13.*is22.*m2.*m3 - 2.*is13.*is23.*m1.*m3 - 2.*is13.*is23.*m2.*m3 + 4.*is11.*is32.*m1.*m2 - 2.*is12.*is31.*m1.*m2 + 4.*is11.*is32.*m1.*m3 - 2.*is12.*is31.*m1.*m3 - 2.*is12.*is32.*m1.*m2 - 2.*is21.*is23.*m1.*m2 - 4.*is11.*is32.*m2.*m3 + 8.*is11.*is33.*m1.*m3 + 2.*is12.*is31.*m2.*m3 + 2.*is12.*is32.*m1.*m3 - 4.*is12.*is33.*m1.*m2 - 4.*is13.*is31.*m1.*m3 + 2.*is13.*is32.*m1.*m2 + 2.*is21.*is23.*m1.*m3 - 2.*is12.*is32.*m2.*m3 + 4.*is12.*is33.*m1.*m3 - 2.*is13.*is32.*m1.*m3 - 2.*is21.*is23.*m2.*m3 + 4.*is12.*is33.*m2.*m3 - 2.*is13.*is32.*m2.*m3 - 2.*is21.*is31.*m1.*m2 - 2.*is21.*is31.*m1.*m3 - 2.*is21.*is32.*m1.*m2 + 4.*is22.*is31.*m1.*m2 + 2.*is21.*is31.*m2.*m3 + 2.*is21.*is32.*m1.*m3 - 4.*is21.*is33.*m1.*m2 - 4.*is22.*is31.*m1.*m3 + 2.*is23.*is31.*m1.*m2 - 2.*is21.*is32.*m2.*m3 + 4.*is21.*is33.*m1.*m3 + 4.*is22.*is31.*m2.*m3 - 2.*is23.*is31.*m1.*m3 + 4.*is21.*is33.*m2.*m3 - 2.*is23.*is31.*m2.*m3 + 8.*is22.*is33.*m2.*m3 - 4.*is23.*is32.*m2.*m3 + 2.*is31.*is32.*m1.*m2 - 2.*is31.*is32.*m1.*m3 - 2.*is31.*is32.*m2.*m3 + 4).^(1/2))./(2.*(is11 + is12 + is13 + is21 + is22 + is23 + is31 + is32 + is33 - 1));

            
        end
        
        function Mmueller = constructMueller2(self,cc,D,mD,mR)
            dimM = size(D);
            dimM(1) = 16;
            
            
            Mmueller = ones(dimM);
%             
            if isfield(self.CalStru,'GPU')
                if self.CalStru.GPU == 1
                    Mmueller = gpuArray(single(Mmueller));
                else
                end
            else
            end
            cc = reshape(cc,1,[]);
            D = reshape(D,3,[]);
            mD = reshape(mD,9,[]);
            mR = reshape(mR,9,[]);
            Mmueller = reshape(Mmueller,16,[]);

            msub = MatrixMultiply(mR,mD);

            Mmueller(1,:) = cc;
            Mmueller(2,:) = sum(mR([1 4 7],:).*D,1);
            Mmueller(3,:) = sum(mR([2 5 8],:).*D,1);
            Mmueller(4,:) = sum(mR([3 6 9],:).*D,1);
            Mmueller(5,:) = D(1,:);
            Mmueller(9,:) = D(2,:);
            Mmueller(13,:) = D(3,:);
            Mmueller(6:8,:) = msub(1:3,:);
            Mmueller(10:12,:) = msub(4:6,:);
            Mmueller(14:16,:) = msub(7:9,:);
            Mmueller = reshape(Mmueller,dimM);

        end
        
        function mask = get_OAxU(self,level)
            r = self.cumulative_r;
            r = permute(r,[2 3 1]);
            rn = r./repmat(max(sqrt(dot(r,r,3)),1e-9),[1,1,3]);    
            square = ones(4);
            square = square./sum(square(:));
            rnf = imfilter(rn,square,'same'); 
            map = sqrt(dot(rnf,rnf,3));
            %mask=map>level;
            SE = strel('disk',3);
            mask = imclose(map>level,SE);
            mask = imopen(mask,SE);
            
            
        end
        
        function imagesc_oa(self,S_plot)
            S_plot = permute(S_plot,[2 3 1]);
            Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);         
            
            %%Sfil = imgaussfilt(Sj1n,2);
           
            
            
            
            CS1n = stokes2color(Sj1n);
            imshow(CS1n) 
        end
        
         
  
    end
end

function m = InCell_logm(M)   
    M = reshape(M,4,4);
    m = logm(M);
    m = m(:);
end   
function m = InCell_expm(M)   
    M = reshape(M,4,4);
    m = expm(M);
    m = m(:);
end  




function error = objective_fun_for_align_vector(v,M1,M2,mask,GPU)

    omeg = [0 v(4) v(5) v(6)
            v(4) 0 -v(3) v(2)
           v(5)  v(3) 0 -v(1)
            v(6) -v(2) v(1) 0];
    rotation = expm(omeg);
    D = [1 0 0 0
         0 1 0 0
         0 0 1 0
         0 0 0 -1];
     
  
    DrotationTD =D*rotation.'*D;

    rotation = rotation(:);
    DrotationTD = DrotationTD(:);
    
    if GPU == 1
        rotation = gpuArray(rotation);
        DrotationTD = gpuArray(DrotationTD);
        M1 = gpuArray(M1);
        M2 = gpuArray(M2);
        mask = gpuArray(mask);        
    end


    Mc = MatrixMultiply(MatrixMultiply(rotation,M1),DrotationTD);
    eM = Mc-M2;
    eM(isnan(eM)) = 0;
    error = double(gather(sum(abs(squeeze(eM).*mask),'all')));

end



