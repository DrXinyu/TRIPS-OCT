classdef ImageWriter < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
       final_folder
       output_folder 
       frame_index
       prefix
    end
    
    methods
        function self = ImageWriter(output_folder,prefix)
            
            self.final_folder = [output_folder '_images' prefix];
            
            self.output_folder = fullfile('C:',['_images' prefix]);
            
            self.prefix = prefix;
            if ~exist(self.output_folder, 'dir')
                 mkdir(self.output_folder);
            end
%             if ~exist(self.final_folder, 'dir')
%                  mkdir(self.output_folder);
%             end
            
            
        end
        
        function self = setFrameIndex(self,num)
            self.frame_index = num;
        end
        
        function im = match_size(self,tomatch,im)
            im = imresize(im,size(tomatch));
        end
        
        
        
        function image = write_intensity_image(self,int,int_r)
            
            destinationFolder = fullfile(self.output_folder,'intensity_images');
            if ~exist(destinationFolder, 'dir')
                mkdir(destinationFolder);
            end
            intInd = mat2gray(20*log10(gather(int)),int_r);
            image = uint8(255 *intInd);

            fname=sprintf('xz_%s_%d.png',self.prefix,self.frame_index);
            fname = fullfile(destinationFolder,fname);
            options.compress  = 'lzw';
            options.overwrite = true;
            options.color = false;
            saveastiff(image, fname, options);
 
        end
        
        function image = write_intensity_image_with_layer(self,int,int_r,layer_data)
            
            destinationFolder = fullfile(self.output_folder,'segmentation_images');
            if ~exist(destinationFolder, 'dir')
                mkdir(destinationFolder);
            end
            intInd = mat2gray(20*log10(gather(int)),int_r);
            image = uint8(255 *intInd);
            imageR = image;
            imageG = image;
            imageB = image;
            [H, W] = size(image);
 
            [n ~] = size(layer_data);
            color = [255 0 0
                     255 128 0
                     128 255 0
                     0 255 255
                     255 0 255];
  
            for index = 1:n
                
                surface = layer_data(index,:);

                if W ~= length(surface)
                    surface =  round(interp1(1:length(self.flat_surface),self.flat_surface,linspace(1,length(self.flat_surface),W)));
                end
                surface(isnan(surface))= round(mean(surface,'omitnan'));
                surface(surface<1) = 1;
                surface(surface>H) = H;
                indxp = sub2ind([H,W],surface,1:W);   
                imageR(indxp) = color(index,1);
                imageG(indxp) = color(index,2);
                imageB(indxp) = color(index,3);
                
                surface = layer_data(index,:)+1;

                if W ~= length(surface)
                    surface =  round(interp1(1:length(self.flat_surface),self.flat_surface,linspace(1,length(self.flat_surface),W)));
                end
                surface(isnan(surface))= round(mean(surface,'omitnan'));
                surface(surface<1) = 1;
                surface(surface>H) = H;
                indxp = sub2ind([H,W],surface,1:W);   
                imageR(indxp) = color(index,1);
                imageG(indxp) = color(index,2);
                imageB(indxp) = color(index,3);
  
                
            end
            
            image = cat(3,imageR,imageG,imageB);
            

            fname=sprintf('xz_%s_%d.png',self.prefix,self.frame_index);
            fname = fullfile(destinationFolder,fname);
            options.compress  = 'lzw';
            options.overwrite = true;
            options.color = true;
            saveastiff(image, fname, options);
 
        end        
        
        
        
        function image = write_color_ret_image(self,int,ret,int_r,ret_r)
            destinationFolder = fullfile(self.output_folder,'color_ret_images');
            if ~exist(destinationFolder, 'dir')
                mkdir(destinationFolder);
            end
            
            intInd = mat2gray(20*log10(gather(int)),int_r);
            
            dnInd = mat2gray(gather(ret),ret_r);                
            dnRGB = ind2rgb(uint8(dnInd*255),flipud(ametrine));
            dnHSV = rgb2hsv(dnRGB);
            dnHSV(:,:,3) = intInd;
            dnRGBmix = hsv2rgb(dnHSV);   
            image = uint8(255 *dnRGBmix);
            
            fname=sprintf('xz_%s_%d.png',self.prefix,self.frame_index);
            fname = fullfile(destinationFolder,fname);
            options.overwrite = true;
            options.color = true;
            options.compress  = 'lzw';
            saveastiff(image, fname, options);            

        end
        
        function image = write_color_dop_image(self,int,dop,int_r,dop_r,name)
            destinationFolder = fullfile(self.output_folder,['color_fdop_images_' name]);
            if ~exist(destinationFolder, 'dir')
                mkdir(destinationFolder);
            end
            intInd = mat2gray(20*log10(gather(int)),int_r);
            dnInd = mat2gray(gather(dop),dop_r);                
            dnRGB = ind2rgb(uint8(dnInd*255),flipud(ametrine));
            dnHSV = rgb2hsv(dnRGB);
            dnHSV(:,:,3) = intInd;
            dnRGBmix = hsv2rgb(dnHSV);   
            image = uint8(255 *dnRGBmix);
            fname=sprintf('xz_%s_%d.png',self.prefix,self.frame_index);
            fname = fullfile(destinationFolder,fname);
            options.overwrite = true;
            options.color = true;
            options.compress  = 'lzw';
            saveastiff(image, fname, options);            

        end        
        
        
        
       function image = write_color_oa_image(self,int,oa,int_r,colorm)
            destinationFolder = fullfile(self.output_folder,'color_oa_images');
            if ~exist(destinationFolder, 'dir')
                mkdir(destinationFolder);
            end
            intInd = mat2gray(20*log10(gather(int)),int_r);
            oauLimit = [-pi pi];
            oauInd = mat2gray(oa,oauLimit);
            
            
            
            oauRGB = ind2rgb(uint8(oauInd*255),colorm);%hsv(256)
            oauHSV = rgb2hsv(oauRGB);
            oauHSV(:,:,3) = intInd;
            oauRGBmix = hsv2rgb(oauHSV);
            image = uint8(255 *oauRGBmix);
            
            fname=sprintf('xz_%s_%d.png',self.prefix,self.frame_index);
            fname = fullfile(destinationFolder,fname);
            options.overwrite = true;
            options.color = true;
            options.compress  = 'lzw';
            saveastiff(image, fname, options);              
            
        end
        
        function image = write_gray_ret_image(self,ret,ret_r)
            destinationFolder = fullfile(self.output_folder,'int_ret_images');
            if ~exist(destinationFolder, 'dir')
                mkdir(destinationFolder);
            end
            
            retInd = mat2gray(gather(ret),ret_r);
            image = uint8(255 *retInd);

            fname=sprintf('xz_%s_%d.png',self.prefix,self.frame_index);
            fname = fullfile(destinationFolder,fname);
            options.overwrite = true;
            options.color = false;
            options.compress  = 'lzw';
            saveastiff(image, fname, options);
 
        end       

        function image = write_ret_oa_image(self,ret,oa,ret_r,colorm)
            
            destinationFolder = fullfile(self.output_folder,'ret_oa_images');
            if ~exist(destinationFolder, 'dir')
                mkdir(destinationFolder);
            end
            
            retInd = mat2gray((gather(ret)),ret_r);
            oauLimit = [-pi pi];
            oauInd = mat2gray(oa,oauLimit);
            
            
            oauRGB = ind2rgb(uint8(oauInd*255),colorm);%hsv(256)
            oauHSV = rgb2hsv(oauRGB);
            oauHSV(:,:,3) = retInd;
            oauRGBmix = hsv2rgb(oauHSV);
            image = uint8(255 *oauRGBmix);

            fname=sprintf('xz_%s_%d.png',self.prefix,self.frame_index);
            fname = fullfile(destinationFolder,fname);
            options.overwrite = true;
            options.color = true;
            options.compress  = 'lzw';
            saveastiff(image, fname, options);
 
        end      
        
        function image = color_bar(self)
            
            color_space = linspace(0,1,256);
            color_space = repmat(color_space)
            
            oauRGB = ind2rgb(uint8(color_space*255),hsv(256));%hsv(256)
            oauHSV = rgb2hsv(oauRGB);
            oauHSV(:,:,3) = retInd;
            oauRGBmix = hsv2rgb(oauHSV);
            image = uint8(255 *oauRGBmix);
            
        end

    
        function wrap_up(self)
            try
                if exist(self.final_folder, 'dir')
                   rmdir(self.final_folder,'s');
                end
            catch
            end
                
            [status,message,messageId] = copyfile(self.output_folder, self.final_folder, 'f')
            rmdir(self.output_folder,'s');

        end
        
        
    end
end

