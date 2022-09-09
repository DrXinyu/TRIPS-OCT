classdef CalibrationManager < handle
    % do everything relating to the CalStru
    properties
        folder
        datafile
        fullfilepath
        
        CalStru
        
        whole_depth = 1.5493e+04;
        pn
        aline_per_frame
        calibration_frame


        workingfolder = []
    end
    
    methods
        function self = CalibrationManager(CalStru,folder,datafile,workingfolder)
            

            %CalStru: the original 
            %

            if nargin < 4
                self.workingfolder = pwd;
            else
                self.workingfolder = workingfolder;
            end
            
            current = pwd;
            if ~isempty(self.workingfolder)
                cd(self.workingfolder)
            end
            
            if ~exist('CalStruCache', 'dir')
               mkdir('CalStruCache')
            end
            
            self.CalStru = CalStru;
            self.datafile = datafile;
            self.folder = folder;
            self.fullfilepath = fullfile(folder,datafile);
            self.aline_per_frame = getAlineNum(self.fullfilepath);
            self.calibration_frame = round(self.aline_per_frame/2);
            
            self.pn = self.CalStru.pn;
            
            if (gpuDeviceCount()) > 0
                self.CalStru.GPU = 1;
            else
                self.CalStru.GPU = 0;
            end
            cd(current)
        end
        
        function set_cali_frame(self,findex)
            self.calibration_frame = findex;
        end
        


        function newCalStru = makeCalStru(self)
            current = pwd;
            if ~isempty(self.workingfolder)
                cd(self.workingfolder)
            end

            if exist(fullfile('CalStruCache',strcat(self.datafile,'.mat')),'file')
                
                nc = load(fullfile('CalStruCache',strcat(self.datafile,'.mat')));
                self.CalStru = nc.self.CalStru;
                
                
                if (gpuDeviceCount()) > 0
                    self.CalStru.GPU = 1;
                else
                    self.CalStru.GPU = 0;
                end
                
                newCalStru = self.CalStru;
                
            else
                
                [fringeA_o,fringeB_o] = self.load_frame(self.calibration_frame,2);
                self.search_aux(fringeA_o,fringeB_o);   
                [fringeA,fringeB] = self.k_resample_setup(fringeA_o,fringeB_o);
                self.dispersion_search(fringeA,fringeB);    
                I = self.get_intensity_image(fringeA,fringeB);
                self.galvo_shift_seatch(I);
                self.depth_encoding_shift_search_cdp(fringeA,fringeB);
                self.save_calstru();
                
%                 self.CalStru = structfun(@gather,self.CalStru,'UniformOutput',false);
%                 save(fullfile('CalStruCache',strcat(self.datafile,'.mat')),"self"); 
                newCalStru = self.CalStru;
                
            end
            cd(current)
        end
        
        function newCalStru = makeCalStru_TM(self)
            current = pwd;
            if ~isempty(self.workingfolder)
                cd(self.workingfolder)
            end


            if exist(fullfile('CalStruCache',strcat(self.datafile,'.mat')),'file')
                
                nc = load(fullfile('CalStruCache',strcat(self.datafile,'.mat')));
                self.CalStru = nc.self.CalStru;
                
                
                if (gpuDeviceCount()) > 0
                    self.CalStru.GPU = 1;
                else
                    self.CalStru.GPU = 0;
                end
                
                newCalStru = self.CalStru;
                
            else
                
                [fringeA_o,fringeB_o] = self.load_frame(1,2);
                [fringeA,fringeB] = self.k_resample_setup_TM(fringeA_o,fringeB_o);
                
                self.inverse_dispersion_poly();
                [fringeA_o,fringeB_o] = self.load_frame(self.calibration_frame,2);
                [fringeA,fringeB,~] = k_stabilize(fringeA_o,fringeB_o,self.CalStru,0);
  
                self.dispersion_search(fringeA,fringeB);    
                I = self.get_intensity_image_TM(fringeA,fringeB);
                self.galvo_shift_seatch(I);
                self.save_calstru();
%                 self.CalStru = structfun(@gather,self.CalStru,'UniformOutput',false);
%                 save(fullfile('CalStruCache',strcat(self.datafile,'.mat')),"self"); 
                newCalStru = self.CalStru;     
            end  
            cd(current)
        end        
        
        function newCalStru = makeCalStru_TM_fast(self)
            current = pwd;
            if ~isempty(self.workingfolder)
                cd(self.workingfolder)
            end

            if exist(fullfile('CalStruCache',strcat(self.datafile,'.mat')),'file')
                
                nc = load(fullfile('CalStruCache',strcat(self.datafile,'.mat')));
                self.CalStru = nc.self.CalStru;
                
                
                if (gpuDeviceCount()) > 0
                    self.CalStru.GPU = 1;
                else
                    self.CalStru.GPU = 0;
                end
                
                newCalStru = self.CalStru;
                
            else
                
                [fringeA_o,fringeB_o] = self.load_frame(1,2);
                [fringeA,fringeB] = self.k_resample_setup_TM(fringeA_o,fringeB_o);
                
                self.inverse_dispersion_poly();
                self.CalStru.ReadShift = 100;
                newCalStru = self.CalStru;     
            end  

            cd(current)
        end            
        
        
        
        
        
        function newCalStru = load_other(self,datafile)

            current = pwd;
            if ~isempty(self.workingfolder)
                cd(self.workingfolder)
            end

            nc = load(fullfile('CalStruCache',strcat(datafile,'.mat')));
            self.CalStru = nc.self.CalStru;


            if (gpuDeviceCount()) > 0
                self.CalStru.GPU = 1;
            else
                self.CalStru.GPU = 0;
            end
            [fringeA_o,fringeB_o] = self.load_frame(self.calibration_frame,2);
            [fringeA,fringeB,~] = k_stabilize(fringeA_o,fringeB_o,self.CalStru,0);
            I = self.get_intensity_image_TM(fringeA,fringeB);
            self.galvo_shift_seatch(I);
            newCalStru = self.CalStru;

            cd(current)
        end
        
        
        function search_aux(self,fringeA_o,fringeB_o)
            self.CalStru = find_aux(fringeA_o,fringeB_o,self.CalStru);
        end
        
        
        
        function save_calstru(self,path)
            
            current = pwd;
            if ~isempty(self.workingfolder)
                cd(self.workingfolder)
            end

            if nargin<2
                self.CalStru = structfun(@gather,self.CalStru,'UniformOutput',false);
                save(fullfile('CalStruCache',strcat(self.datafile,'.mat')),"self"); 
            else
                self.CalStru = structfun(@gather,self.CalStru,'UniformOutput',false);
                save(fullfile(path,'CalStruCache',strcat(self.datafile,'.mat')),"self"); 
            end
            
            cd(current)
        end
        
        
        function inverse_dispersion_poly(self)
            
            CPhase = unwrap(-angle(1./self.CalStru.CArray));
            
            xo = (20:length(CPhase)-20)';
            CPV = gather(CPhase(xo));
            normxo = (xo-mean(xo))/std(xo);

            % [p_kv,~,mu] = polyfit(xo,mad,20);
            % mad_fit = polyval(p_kv,(1:length((MAmean)))',[],mu);
            dispersion_fit_x = ((1:length((CPhase)))'-mean(xo))/std(xo);

            dispersion_poly = polyfit(normxo,CPV,4);
            dp_fit = polyval(dispersion_poly,dispersion_fit_x);



            dispersion_fitting.poly = dispersion_poly;
            dispersion_fitting.x = dispersion_fit_x;

            self.CalStru.CArray =  (exp(-(1i.*dp_fit)));
            self.CalStru.dispersion_fitting = dispersion_fitting;
            
        end
        
        
        function [fringeA_o,fringeB_o] = load_frame(self,index,num)
            
            fid1 = fopen(self.fullfilepath,'r', 'b');
            if fid1 < 0
                disp('file open failed');
            end
            if isfield(self.CalStru,'ReadShift')
                offset = self.CalStru.ReadShift;  
            else
                offset = 0;
            end   

            fseek(fid1, self.pn*self.aline_per_frame*index*2-self.pn*offset*2, 'bof');
            B1 = fread(fid1, [self.pn, self.aline_per_frame*num], 'uint16','n');
            
            fclose(fid1);
            fringeA_o = double(B1(1:2:end,:));
            fringeB_o = double(B1(2:2:end,:));
            
            
        end
        
        function [fringeA,fringeB] = detector_calibrate(self,fringeA,fringeB)
%             fringeA = fringeA;
%             fringeB = fringeB;
        end        
           
        
        function [fringeA,fringeB] = k_resample_setup(self,fringeA,fringeB)
            [fringeA,fringeB,self.CalStru] = k_stabilize_setup(fringeA,fringeB,self.CalStru);
        end           
        
        function [fringeA,fringeB] = k_resample_setup_TM(self,fringeA,fringeB)
            [fringeA,fringeB,self.CalStru] = k_stabilize_setup_TM(fringeA,fringeB,self.CalStru);
        end    
%         
        
        
        function dispersion_search(self,fringeA,fringeB)
            self.CalStru = retrieve_dispersion(fringeA,fringeB,self.CalStru);  
        end             
        
        function depth_encoding_shift_search(self,fringeA,fringeB)
            [~,self.CalStru] = alignShift(fringeA(:,1:self.aline_per_frame),fringeB(:,1:self.aline_per_frame),self.CalStru);
        end  
        
        function depth_encoding_shift_search_cdp(self,fringeA,fringeB)
            [~,self.CalStru] = alignShift_compact_depth_encoding(fringeA(:,1:self.aline_per_frame),fringeB(:,1:self.aline_per_frame),self.CalStru);
        end
        
        
        function I = get_intensity_image(self,fringeA,fringeB)
            
            self.CalStru.Window = 1;
            self.CalStru.MinDepth = 1;
            self.CalStru.Shifting = 0;
            self.CalStru.Binning = 0;
            self.CalStru.Complex = 0;
            self.CalStru.Dispersion = 1;            
            self.CalStru.Background = 'mean';
            
            I1 = fringe2image(fringeA,self.CalStru);
            I2 = fringe2image(fringeB,self.CalStru);

            self.CalStru.Shifting = 1;
            
            IS1 = fringe2image(fringeA,self.CalStru);
            IS2 = fringe2image(fringeB,self.CalStru);

            self.CalStru.Shifting = 0;
            I = I1+I2+IS1+IS2;
            
        end
        
        function I = get_intensity_image_TM(self,fringeA,fringeB)
            
            self.CalStru.Window = 1;
            self.CalStru.MinDepth = 1;
            self.CalStru.Shifting = 0;
            self.CalStru.Binning = 0;
            self.CalStru.Complex = 0;
            self.CalStru.Dispersion = 1;            
            self.CalStru.Background = 'mean';
            
            I1 = fringe2image(fringeA,self.CalStru);
            I2 = fringe2image(fringeB,self.CalStru);
            
            I = I1+I2;
            
        end

        function galvo_shift_seatch(self,I)
            [~,self.CalStru] = shift_galvo_seatch(I,self.aline_per_frame,self.CalStru);
        end  
        
        function CalStru = add_frame_stabilization_data(self,frame_stable_manager)
            
            self.CalStru = frame_stable_manager.CalStru;
            self.CalStru.FrameShift = cumsum(frame_stable_manager.shifts);
            self.save_calstru();
            CalStru = self.CalStru;
            
        end

        function CalStru = add_cFactor(self,cFactor)

            self.CalStru.cFactor = cFactor;
            self.save_calstru();
            CalStru = self.CalStru;


        end
        
        function CalStru = add_PS_data(self,PS_manager)
            
            self.CalStru = PS_manager.getCalStru();            
            self.save_calstru();
            CalStru = self.CalStru;
        
        end        
        
                
        function CalStru = add_surface(self,SurfaceManager,location,path)
            
 
            
            face_data = struct;
            face_data.face = SurfaceManager.surface;
            face_data.index = SurfaceManager.findex;
            
            if strcmpi(location,'RNFL')
                self.CalStru.RNFLface = face_data;  
            elseif strcmpi(location,'RPE')
                self.CalStru.RPEface = face_data;                  
            elseif strcmpi(location,'Sclera')
                self.CalStru.Scleraface = face_data;   
            elseif strcmpi(location,'outerSclera')
                self.CalStru.Scleraface = face_data;   
            else
                
            end
            
            
            if nargin<4
                self.save_calstru();
            else
                self.save_calstru(path);
            end
            
            CalStru = self.CalStru;
            
            
            
        end            

        
    end
end

