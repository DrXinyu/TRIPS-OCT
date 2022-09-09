classdef PSTriple_Manager_v2 < handle
    % this class use triple input system to solve retardation and optic axis
    
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
        cumulative_r
        
        nan_mask
        
        MMc
        mask
        retina_mask
        last_surface_WW = 0;
        
        
        ret_level = 4;
        
    end
    
    methods
        
        function self = PSTriple_Manager_v2(CalStru)
            %% in this function, do orthonolize
            self.pixel_depth = CalStru.whole_depth/CalStru.NFFT;
            self.CalStru = CalStru;
            %self.zfilt_kernal = zfilt_kernal;
            if isfield(CalStru,'binning_rc')
                self.Mcorr = self.CalStru.binning_rc.Mcorr;
                self.PMDCorr = self.CalStru.binning_rc.PMDCorr;
                self.DPMDCorrTD = self.CalStru.binning_rc.DPMDCorrTD;                
            end

            
        end
 
        %% load_Stokes>>tri_input
        function load_Stokes(self,S1,S2,S3)
            
            % I copied Martin's Mueller matrix reconstruction code here

            self.dim = size(S1);
            self.sdim = length(self.dim);
            self.bins = self.dim(self.sdim-1);
                       
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

            S1 = permute(S1,[self.sdim, 1:self.sdim-1]);
            S2 = permute(S2,[self.sdim, 1:self.sdim-1]);
            S3 = permute(S3,[self.sdim, 1:self.sdim-1]);
            
            self.I1 = sqrt(sum(S1.^2,1));
            self.I2 = sqrt(sum(S2.^2,1));
            self.I3 = sqrt(sum(S3.^2,1));

%            compute initial c-values
            c12 = sqrt(abs(self.I1.*self.I2 - sum(S1.*S2,1)));
            c13 = sqrt(abs(self.I1.*self.I3 - sum(S1.*S3,1)));
            c23 = sqrt(abs(self.I2.*self.I3 - sum(S2.*S3,1)));

            % (arbitrarily) define their mean
            cEst = (c12 + c13 + c23)/3;

            % find how S1, S2, and S3 have to be scaled to achieve cEst
            alpha = cEst.*c23./c12./c13;
            beta = cEst.*c13./c12./c23;
            gamma = cEst.*c12./c13./c23;

            S1 = alpha.*S1;
            S2 = beta.*S2;
            S3 = gamma.*S3;

            I1 = sqrt(sum(S1.^2,1));
            I2 = sqrt(sum(S2.^2,1));
            I3 = sqrt(sum(S3.^2,1));

            x = cat(1,I1,I2,I3)./cEst;
            
%             x1 = I1./cEst;
%             x2 = I2./cEst;
%             x3 = I3./cEst;

%             Msub = cat(1,S1,S2,S3);   
%             aaa = det3x3(Msub);

            ddet = S1(1,:).*S2(2,:).*S3(3,:) + S1(2,:).*S2(3,:).*S3(1,:) + S1(3,:).*S2(1,:).*S3(2,:) - S1(3,:).*S2(2,:).*S3(1,:) - S1(1,:).*S2(3,:).*S3(2,:) - S1(2,:).*S2(1,:).*S3(3,:);
            cha = sum(x,1)/2 - sqrt(abs(sum(x,1).^2-2*sum(x.^2,1)-2))/2;
            cha2 = sum(x,1)/2 + sqrt(abs(sum(x,1).^2-2*sum(x.^2,1)-2))/2;
            cha(ddet<0) = cha2(ddet<0);
            
            cha(cha<=1) = 1.0001;
            % normalized a vector
            an = (x-cha)./sinh(acosh(cha));

            dim = size(an);
            mDinv = reshape(eye(3),9,1) - (cat(1,an,an,an).*tanh(acosh(cha)) + reshape(cat(1,an(:,:).*an(1,:),an(:,:).*an(2,:),an(:,:).*an(3,:)).*(1-1./cha(:,:)),[9,dim(2:end)]))./(1+sum(an,1).*tanh(acosh(cha)));    

            mR = MatrixMultiply(cat(1,S1,S2,S3)./cEst,mDinv);
            
%             out.D = an.*tanh(acosh(cha));
%             out.R = decomposeRot(mR);

            mD33 = reshape(eye(3),9,1) + (cat(1,an(:,:).*an(1,:),an(:,:).*an(2,:),an(:,:).*an(3,:)).*(cha(:,:)-1));
            self.Mmueller = self.constructMueller(cha,sinh(acosh(cha)).*an,mD33,mR).*cEst;

            self.nan_mask = ~isnan(mean(self.Mmueller,1));

        end
        

        
        
        function Mmueller = constructMueller(self,cc,D,mD,mR)
            dimM = size(D);
            dimM(1) = 16;
            Mmueller = ones(dimM);
            
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
        
        
        
        
        
        
        function look_nan(self,m)
            figure
            for bindex = 1:self.bins
                subplot(1,self.bins,bindex)
                imagesc((squeeze(m(1,:,:,bindex))))
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

        function errFinal = obtainSymmetrizationCorrection(self,imask)

            %imask = permute(imask,[3 1 2 4]);
            
            mask = repmat(self.nan_mask,[16,1,1,1]);
            mask = mask.*imask;
            
            for bindex = 1:self.bins
                MM = self.Mmueller(:,:,:,bindex);
                maskbin = mask(:,:,:,bindex);
                selectedMueller = reshape(MM(maskbin>0),16,[]);
                selectedMueller = selectedMueller./(selectedMueller(1,:));
                MMproj = mean(selectedMueller,2,'omitnan');
                MMprojsave(:,bindex) = MMproj;
                
                MMproj(isnan(MMproj))=1;
                
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
                errFinal = self.Mcorr(1,:)
                
            end

%             %unwrap correction across spectral bins
%             ss = cat(2,0,mod(cumsum(sqrt(sum(diff(wcorr,[],2).^2,1))>pi),2));
%             if sum(ss)>self.bins/2
%                 ss = 1-ss;
%             end
%             retcorr = sqrt(sum(wcorr.^2,1));
%             wcorr = wcorr - bsxfun(@times,bsxfun(@rdivide,wcorr,retcorr),ss*2*pi);
% 
%             % maintain a retardation <pi
%             if mean(sqrt(sum(wcorr.^2,1)))>pi
%                 wcorr = wcorr - bsxfun(@rdivide,wcorr,sqrt(sum(wcorr.^2,1)))*2*pi;
%             end
%             
%             figure
%             plot(wcorr')
%             figure
%             plot(ddcorr')

        end
        
        
        
        
        function s = searchSymmetrizationCorrection(self,imask)
            
            
            
            mask = repmat(self.nan_mask,[16,1,1,1]);
            mask = mask.*imask;
            
            v0 = [0 0 0 0 0 0];            
% 
            cbi = round(self.bins/2);
% 
%             v(:,cbi) = [0 0 0 0 0 0];
            
            
            for bindex = 1:self.bins            

                MM = self.Mmueller(:,:,:,bindex);
                maskbin = mask(:,:,:,bindex);
                selectedMueller = reshape(MM(maskbin>0),16,[]);
                selectedMueller = selectedMueller./(selectedMueller(1,:));
                MMproj = mean(selectedMueller,2,'omitnan');
%                 MMprojsave(:,bindex) = MMproj;
                
                MMproj(isnan(MMproj))=0;
                selectedMueller(isnan(selectedMueller))=0;
                
                
                A = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0]/sqrt(2);
                F = conj(A'*reshape(MMproj,[4,4])*A);
                H = reshape(F([1,9,3,11,5,13,7,15,2,10,4,12,6,14,8,16]'),[4,4]); 
                [a,b] = eig(diag([-1,1,-1,1])*H*diag([-1,1,-1,1]));
                [~,pos] = min(abs(diag(b)));
                cc = a(:,pos);
                Jcorr = [cc(2),cc(4);cc(1),cc(3)];
                Jcorr = Jcorr/sqrt(det(Jcorr));
                [rr,dd] = JonesDecomp(Jcorr);% convert back to Stokes, ignoring diattanuation component
                
                v0 = double(gather([(rr(:))',(dd(:))']));

                bindex
                %MMproj = gather(MMproj);
                v(:,bindex) = self.search_for_sym_vectors(v0,selectedMueller);
                
            end
            
%             v0 = v(:,cbi); 
%             for bindex = cbi-1:-1:1 
%                 MM = self.Mmueller(:,:,:,bindex);
%                 maskbin = mask(:,:,:,bindex);
%                 selectedMueller = reshape(MM(maskbin>0),16,[]);
%                 selectedMueller = selectedMueller./(selectedMueller(1,:));
%                 MMproj = mean(selectedMueller,2,'omitnan');
% %                 MMprojsave(:,bindex) = MMproj;
%                 
%                 MMproj(isnan(MMproj))=0;
%                 selectedMueller(isnan(selectedMueller))=0;
%                 
%                 bindex
%                 %MMproj = gather(MMproj);
%                 v(:,bindex) = self.search_for_sym_vectors(v0,selectedMueller);
%                 v0 = v(:,bindex);
%             end

            for bindex = 1:self.bins   
                omeg = [0 v(4,bindex) v(5,bindex) v(6,bindex)
                        v(4,bindex) 0 -v(3,bindex) v(2,bindex)
                       v(5,bindex)  v(3,bindex) 0 -v(1,bindex)
                        v(6,bindex) -v(2,bindex) v(1,bindex) 0];                
%                 omeg = [0 0 0 0
%                         0 0 -v(3,bindex) v(2,bindex)
%                        0  v(3,bindex) 0 -v(1,bindex)
%                         0 -v(2,bindex) v(1,bindex) 0];
                rotation = expm(omeg);
                s(bindex) = rotation(1);
                self.Mcorr(:,bindex) = rotation(:);
            end


        end
        
        function v = search_for_sym_vectors(self,v0,M)
            
            f = @(x)objective_fun_for_symmetrization(x,M,self.CalStru.GPU);
            options = optimset('MaxFunEvals',10000);
%             options = optimset('Display','iter','MaxFunEvals',3000,'PlotFcns',@optimplotfval);
%             v = fminsearch(f,v0,options);
            v = fminsearch(f,v0,options);

        end     
    
        function checkSymmetrizationCorrection(self)
            v = [6 7 8]; %[5 9 13];%  %[6 10 14] % 8 12 16
            vT = [6 10 14];%;[2 3 4];%  %[6 7 8] % 14 15 16
            
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
            imask = permute(imask,[3 1 2 4]);
            mask = repmat(self.nan_mask,[16,1,1,1]);
            mask2 = repmat(imask,[16,1,1,9]);

            imask = (mask.*mask2)>0;
            
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
                
%                 omeg = [0 0 0 0
%                         0 0 -v(3,indw) v(2,indw)
%                        0  v(3,indw) 0 -v(1,indw)
%                         0 -v(2,indw) v(1,indw) 0];
%                 
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
%             options = optimset('Display','iter','MaxFunEvals',3000,'PlotFcns',@optimplotfval);
%             v = fminsearch(f,v0,options);
            v = fminsearch(f,v0);

        end        
        
        
        function flip_Mueller(self)
            self.Mmueller = flip(self.Mmueller,2);
        end
        

        function checkPMDCorrection(self)
            v = [6 7 8];%[6 10 14]; %
            vT = [14 15 16];%;[2 3 4];%
            
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
            %d_part(3,:,:,:) = -d_part(3,:,:,:);
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
            %d_part(3,:,:,:) = -d_part(3,:,:,:);
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
 % get mR by polar decomposition            
            r_pd = MuellerJonesDecomp(self.Mmueller,1);
            
            mR = makeRot(r_pd);  
            %mR = self.Mmueller([6 7 8 10 11 12 14 15 16],:,:);
            %mR = self.polarDecomposition(self.Mmueller);
            %
%             mR = self.Mueller_gaussfilt(mR,3);
%             mR = euclideanRotation(mR);
             
%                 S_plot = squeeze(mR(1:3,:,:));
%                 S_plot = permute(S_plot,[2 3 1]);
%                 Sj1n = S_plot./repmat(max(sqrt(dot(S_plot,S_plot,3)),1e-9),[1,1,3]);                    
%                 CS1n = stokes2color(Sj1n);
%                 figure;
%                 imshow(CS1n) 
            
%             mR = self.polarDecomposition(self.Mmueller); 
            % shift 1 pixel along depth
            mRdz = circshift(mR,-1,2);
            
            % calculate locat rotation matrix
            mLocal = MatrixMultiply(mR,matrix3inv(mRdz));   
            %aaa = sum(mLocal([1,5,9],:,:,:),1)/2-1/2;
            % find the length of the rotation vector
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
            dopmask = (self.get_depolarization_index()>0.750)& (self.get_depolarization_index()<1);
            dopmask = squeeze(dopmask(1,:));  

            
            filtwidth = 30;
            surfMuellerFilt = imfilter(surface_Mueller,ones(1,filtwidth)/filtwidth);

            
            % some filtering method
            %surface_Mueller = self.Mulluer_fit_filter(surface_Mueller,dopmask);
            
            for ind = 1:size(surfMuellerFilt,2)
                H = M2H(reshape(surfMuellerFilt(:,ind),[4,4]));
                [a,b] = eig(H);
                [~,mp] = max(diag(b));
                surfMuellerFilt(:,ind) = reshape(H2M(a(:,mp)*a(:,mp)'),[16,1]);
            end
            
            [surfrf,surfdf] = MuellerJonesDecomp(surfMuellerFilt);
            
            
            
            %[surfrf,surfdf] = self.unwrapRetDiatSurf(surfrf,surfdf,dopmask);
            
            % I tried to unwrap, not working
                                %%%%%%%%%%%%%%%%%
                    %             ret = ret;
                    %             thit = sqrt(dot(ret,ret,1));
                    %             thit_f = medfilt1(thit,5,[],2);
                    %             ret = (thit_f./thit).*ret;
                    %            ret = self.layer_rotV_align(ret,dopmask);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % reconstruct the Mueller matrix
            
            inv_half_surface_Mueller = self.Jones2Mueller(makeJones(-surfrf/2,-surfdf/2));
            %self.slow_vectors2Mueller(-surfrf./2,-surfdf./2);
            
            if isfield(self.CalStru,'GPU')
                if self.CalStru.GPU == 1
                    inv_half_surface_Mueller= gpuArray(single(inv_half_surface_Mueller));
                else
                end
            else
            end
            
            % apply Mueller matrix
            inv_half_surface_Mueller = repmat(permute(inv_half_surface_Mueller,[1 3 2]),1,dim(2),1);
            self.Mmueller = MatrixMultiply(inv_half_surface_Mueller,MatrixMultiply(self.Mmueller,inv_half_surface_Mueller));
            
        end
        
        function remove_circular_diatenuation(self)
            
            % normalize the Mueller
            
            dim = size(self.Mmueller);
            mask = self.get_depolarization_index()>0.6;
            
            normalize_Mueller = (self.Mmueller(:,:,:)./self.Mmueller(1,:,:)).*permute(mask,[3 1 2]);%
            normalize_Mueller(isnan(normalize_Mueller)) = 0;
            
            % sum the Mueller matrix along the aline
            frame_M = (sum(normalize_Mueller,2));
            mask = (sum(mask,2));

            % calculate the circular diatteuation amount to be corrected
            tanh2z = (-(2.*frame_M(13,:,:))./(frame_M(1,:,:)-frame_M(16,:,:)));
            
            tanh2z(tanh2z>1) = 0.99;
            tanh2z(tanh2z<-1) = -0.99;
            
            z = atanh(tanh2z)/2;
            z = medfilt1(gather(z),30,[],3);
            z(isnan(z)) = 0;
            
            
            % construct the circular diatteuation Mueller matrix 
            [MDC,invMDC] = self.circularDiattenuator(z);

            MDC = repmat(MDC,1,dim(2),1);
            invMDC = repmat(invMDC,1,dim(2),1);
            
            
            % apply the correction
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
        
        
        
        
        
        
        
        function [retFinal,phi] = local_optic_axis(self,iteration_depth)
            %%
            % generate mask to exclude points with low DOP from correcting the local
            % retardation
            dimm = size(self.Mmueller); 
            
            if nargin > 1
              iteration_depth = iteration_depth;
            else
              iteration_depth = dimm(2);
            end

%             dopmask = (self.get_depolarization_index()>0.97)&(self.get_depolarization_index()<1);
%             surface_mask = dopmask(1,:);  
            
            [r_pd,d_pd] = MuellerJonesDecomp(self.Mmueller,1);

            mR_pd = makeRot(r_pd);
            
            
            mR_pddz = circshift(mR_pd,-1,2);
            self.cumulative_r = decomposeRot(mR_pd);  
            
            mLocal_pd = MatrixMultiply(mR_pd,matrix3inv(mR_pddz));     
            r_pd = decomposeRot(mLocal_pd); 
            
%             mLocal_cd = euclideanRotation(mLocal_cd);            
%             self.cumulative_r = decomposeRot(mR);   
%             r_cd = decomposeRot(mLocal_cd);   
                         
            MMM = mR_pd;%self.Mmueller([6,7,8,10,11,12,14,15,16],:,:);
  
            %Mm = self.Mulluer_fit_filter(Mm,surface_mask);
            
            
            
%              Mm = (squeeze(MMM(:,1,:)));
%              Mm = self.Mulluer_fit_filter(Mm,surface_mask);
% %             Mm = squeeze(MMM(:,1,:,:));
%              Mm = euclideanRotation(Mm);
% % %             
%               WFit = decomposeRot(Mm);
% %             WFit(3,:,:) = 0;                        
% % % 
%               WW = self.layer_rotV_align(WFit,ones(size(WFit)));
% % %           
%               N = permute(makeRot(WW/2),[1,3,2]);
             
            %MMM = euclideanRotation(MMM);
            
            % GPU setup
            nanmask = squeeze(sum(isnan(MMM),1));
            
            if isfield(self.CalStru,'GPU')              
                if self.CalStru.GPU == 1                
                    w = gpuArray(single(zeros([3,dimm(2),dimm(3)])));
                else
                    w = zeros([3,dimm(2),dimm(3)]);
                end
            else

            end


            mmask = ones(dimm(2),dimm(3));%dopmask;%

            N = repmat([1;0;0;0;1;0;0;0;1],[1,1,dimm(3)]);
            
            if isfield(self.CalStru,'GPU')              
                if self.CalStru.GPU == 1                
                    N = gpuArray(single(N));
                else
                end
            else
            end            
            
%             N = permute(makeRot(cat(1,zeros(2,numel(WW)),-WW)),[1,3,2]);
      
            % layer-by-layer correction
            
            
            for indz = 2:iteration_depth%refInds(3):roiz(end)

                % D*N(indz)*D*Mtot(indz+1)*N'(indz)
                nloc = MatrixMultiply(bsxfun(@times,N,[1;1;-1;1;1;-1;-1;-1;1]),MatrixMultiply(MMM(:,indz,:,:),N([1,4,7,2,5,8,3,6,9],:,:,:)));
                Wloc = decomposeRot(nloc);   
                w(:,indz,:) = squeeze(Wloc);  
                nlocroot = makeRot(Wloc/2);

                Nnext = MatrixMultiply(nlocroot,N);
                N(:,:,mmask(indz,:)>0) = Nnext(:,:,mmask(indz,:)>0);

            end
            
            %Omega = bsxfun(@times,permute(pa,[2,3,1]),ret);
            %Omegaf = imfilter(permute(w,[2,3,1]),ones(self.fwz,1)/self.fwz);
            %Omegaf = imfilter(permute(w,[2,3,1]),self.zfilt_kernal);
            Omegaf = permute(w,[2,3,1]);
            
            retFinal = squeeze(sqrt(sum(r_pd.^2,1))/2/self.pixel_depth/pi*180);
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
        
        
        function Mueller_medfilt(self,kernel)
            if nargin<2
                kernel = [3 3];
            end
            [m,H,W,A,B] = size(self.Mmueller);
            
            for index = 1:m
                
                S = (self.Mmueller(index,:,:,:));
            	S = reshape(S,H,[]);
                
                S = medfilt2(S,kernel);
                S = reshape(S,1,H,W,A,B);
                self.Mmueller(index,:,:,:) = S;
            end
        end
        
        function M=Mueller_gaussfilt(self,M,sigma)
            if nargin<2
                sigma = 2;
            end
            [m,H,W,A,B] = size(M);
            for index = 1:m
                S = M(index,:,:,:);
            	S = reshape(S,H,[]);
                S = imgaussfilt(S,sigma);
                S = reshape(S,1,H,W,A,B);
                M(index,:,:,:,:) = S;
            end
        end        
        
        function [retuw,diatuw] = unwrapRetDiatSurf(self,ret,diat,mask,init)
            % takes the ret and diat of the surface signal (or any other 1D signal of
            % ret/diat signal) and 'unwraps' it, starting at pixel init (defaults to
            % center)

            if nargin<5
                init = round(size(ret,2)/2);
            end

            % complex valued ret/diat vector
            q = ret + 1i*diat;
            
            q1 = q;
            q2 = q;
            q3 = q;
            q4 = q;
            
            inds = mod(cumsum(cat(2,0,abs(sqrt(sum(diff(q,[],2).^2,1)))>pi)),2)>0;

            q1(:,inds) = q1(:,inds) - 2*pi*q1(:,inds)./sqrt(sum(q1(:,inds).^2,1));
            q2(:,inds) = q2(:,inds) + 2*pi*q2(:,inds)./sqrt(sum(q2(:,inds).^2,1));

            inds = ~inds;
            q3(:,inds) = q3(:,inds) + 2*pi*q3(:,inds)./sqrt(sum(q3(:,inds).^2,1));
            q4(:,inds) = q4(:,inds) - 2*pi*q4(:,inds)./sqrt(sum(q4(:,inds).^2,1));            
            
            dia_amount(1) = sum(sum(imag(q1(:,mask)).^2));
            dia_amount(2) = sum(sum(imag(q2(:,mask)).^2));
            dia_amount(3) = sum(sum(imag(q3(:,mask)).^2));
            dia_amount(4) = sum(sum(imag(q4(:,mask)).^2));
            
            [~,mind] = min(dia_amount);
            qq = {q1,q2,q3,q4};
            
            q =qq{mind};
            
            n = -1:1;
            for ind = 1:3
                iq = q + 2*pi*n(ind)*bsxfun(@rdivide,q,sqrt(sum(q.^2,1)));
                meanRet(ind) = mean(mean(sqrt(sum((iq-self.last_surface_WW).^2,1))));
            end
            [~,ri] = min(meanRet);            
            q = q + 2*pi*n(ri)*bsxfun(@rdivide,q,sqrt(sum(q.^2,1)));
            
            self.last_surface_WW = q;
            retuw = real(q);
            diatuw = imag(q);
            
            
        end
        
        function WW = layer_rotV_align(self,WFit,surf)

%             ss = repmat(cat(2,0,mod(cumsum(sqrt(sum(diff(WFit,[],2).^2,1))>pi,2),2)),[3,1])>0;
%             delta = 2*pi*bsxfun(@rdivide,WFit,sqrt(sum(WFit.^2,1)));
%             WFitc = WFit;
%             WFitc(ss) = WFit(ss) - delta(ss);
            WFitc = WFit;
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
            dimm(1) = 1;
            dopi(isnan(dopi)) = 0;
            dopi(dopi>1) = 0;
            dopi = squeeze(reshape(dopi,dimm));
            
        end
        
        
        function retina_mask = get_retina_mask(self,dop,Isnr,ret)
            
            dopi = self.get_depolarization_index();
            retina_mask1 = (dop>0.4)&(dopi>0.5)&(dopi<1);% 
            retina_mask = (retina_mask1&(Isnr>1.5)&(ret<4));
            SE = strel('disk',1);
            mask = imerode(retina_mask,SE);
            retina_mask(150:end,:)=0;
            
        end

        function retina_mask = get_human_retina_mask(self,dop,Isnr,ret)
            
            dopi = self.get_depolarization_index();
            retina_mask1 =(dopi>0.3);% 0.6
            retina_mask = (retina_mask1&(Isnr>1.1)&(ret<3));
            SE = strel('disk',1);
            retina_mask = imerode(retina_mask,SE);
%             
        end
        
%         function sclera_60um_mask = get_sclera_60um_mask(self,dop,Isnr,ret)
%             
%             dopi = self.get_depolarization_index();
%             retina_mask1 =(dopi>0.6);% 
%             retina_mask = (retina_mask1&(Isnr>1.3)&(ret<5));
%             SE = strel('disk',1);
%             mask = imerode(retina_mask,SE);
%             
%         end
        
        
        function retina_mask = get_PGP_retina_mask(self,dop,Isnr,ret)
            
            dopi = self.get_depolarization_index();
            retina_mask1 = (dop>0.20)&(dopi>0.5)&(dopi<1);% 
            retina_mask = (retina_mask1&(Isnr>1.2)&(ret<5));
            
            SE = strel('disk',1);
            mask = imerode(retina_mask,SE);
            
            retina_mask(120:end,:)=0;
        end        
        
        

        function array = zfiltering(self,array)
            array = permute(array,[2 3 4 1]);
            dim = size(array);
            array = reshape(array,dim(1),dim(2),[]);
            array = imfilter(array,self.zfilt_kernal,'same');
            array = reshape(array,dim);
            array = permute(array,[4 1 2 3]);

        end

        
       function error = checkStokes(self,S1,S2,S3)
        % I copied martin's checking code here
         
            S1=gather(S1);
            S2=gather(S2);
            S3=gather(S3);
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
                a1(ind) = mean(angle1(mask(:,:,ind)),'all'); 
                a2(ind) = mean(angle2(mask(:,:,ind)),'all'); 
                a3(ind) = mean(angle3(mask(:,:,ind)),'all'); 
            end
%             figure(1)
%             plot(a1)
%             hold on 
%             plot(a2)
%             plot(a3)
%             hold off
            
            error = var(a1)+var(a2)+var(a3);
            
            
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

%             figure(2)
%             clf
%             subplot(2,4,1)
%             %imagesc(10*log10(I))
%             title('Intensity')
% 
%             subplot(2,4,5)
%             imagesc(mask(:,:,indw))
%             title(sprintf('Intensity bin %d',indw))
% 
%             subplot(2,4,[2,3,4])
%             imagesc(cat(2,S1n(:,:,indw,1),S1n(:,:,indw,2),S1n(:,:,indw,3)))
%             title(sprintf('QUV bin %d',indw))
% 
%             subplot(2,4,6)
%             imagesc(dotProd12(:,:,indw))
%             title(sprintf('Dot product 12, bin %d',indw))
% 
%             subplot(2,4,7)
%             imagesc(dotProd13(:,:,indw))
%             title(sprintf('Dot product 13, bin %d',indw))
% 
%             subplot(2,4,8)
%             imagesc(dotProd23(:,:,indw))
%             title(sprintf('Dot product 23, bin %d',indw))
% 
% 
%             figure(3)
%             clf
%             subplot(1,3,1)
%             imagesc(HH12)
%             title('Histograms dotproduct 12')
% 
%             subplot(1,3,2)
%             imagesc(HH13)
%             title('Histograms dotproduct 13')
% 
%             subplot(1,3,3)
%             imagesc(HH23)
%             title('Histograms dotproduct 23')


            
            
        end
        
        
%         function mask = get_OAxU(self,level)
%             r = self.cumulative_r;
%             r = permute(r,[2 3 1]);
%             rn = r./repmat(max(sqrt(dot(r,r,3)),1e-9),[1,1,3]);    
%             square = ones(4);
%             square = square./sum(square(:));
%             rnf = imfilter(rn,square,'same'); 
%             map = sqrt(dot(rnf,rnf,3));
%             mask=map>level;
% %             SE = strel('disk',3);
% %             mask = imclose(map>level,SE);
% %             mask = imopen(mask,SE);
%             
%             
%        end
        
        
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


function error = objective_fun_for_symmetrization(v,M,GPU)
    omeg = [0 v(4) v(5) v(6)
            v(4) 0 -v(3) v(2)
           v(5)  v(3) 0 -v(1)
            v(6) -v(2) v(1) 0];
%     omeg = [0 0 0 0
%             0 0 -v(3) v(2)
%            0  v(3) 0 -v(1)
%             0 -v(2) v(1) 0];
    rotation = expm(omeg);
    
    if GPU == 1
        rotation = gpuArray(rotation);

        M = gpuArray(M);        
    end

    Mc = MatrixMultiply(rotation(:),M);
    
    
    McT = transpose4x4(Mc);
    McT([4 8 12 13 14 15],:,:,:) = - McT([4 8 12 13 14 15],:,:,:);

    error = double(gather(sum(abs(Mc-McT),'all')));
end

function error = objective_fun_for_align_vector(v,M1,M2,mask,GPU)

    omeg = [0 v(4) v(5) v(6)
            v(4) 0 -v(3) v(2)
           v(5)  v(3) 0 -v(1)
            v(6) -v(2) v(1) 0];
%         
%        omeg = [0 0 0 0
%             0 0 -v(3) v(2)
%            0  v(3) 0 -v(1)
%             0 -v(2) v(1) 0];     
        
        
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
    
%         M = (Mc+M2)/2;
%         M = reshape(M,16,[]);
%         M = M./M(1,:);
%         T1 = M(1,:).^2+M(2,:).^2+M(3,:).^2+M(4,:).^2;
%         T2 = M(5,:).^2+M(6,:).^2+M(7,:).^2+M(8,:).^2;
%         T3 = M(9,:).^2+M(10,:).^2+M(11,:).^2+M(12,:).^2;
%         T4 = M(13,:).^2+M(14,:).^2+M(15,:).^2+M(16,:).^2;
%         dopi = sqrt((T1+T2+T3+T4-1)/3);
%         dimm(1) = 1;
%         dopi(isnan(dopi)) = 0;
%         dopi(dopi>1) = 0;
% 
%     error = double(gather((-sum(dopi,'all'))));
   
    eM = Mc-M2;
    eM(isnan(eM)) = 0;
    error = double(gather(sum(abs(squeeze(eM).*mask),'all')));

end



