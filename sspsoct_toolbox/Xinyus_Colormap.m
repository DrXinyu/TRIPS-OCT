classdef Xinyus_Colormap < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function self = Xinyus_Colormap()
        end
        
        
        function show_color(self,map)
            
            colorbar = repmat(permute(map,[1 3 2]),1,20,1);
            
            figure()
            subplot(1,2,1)
            imshow(colorbar)    
            subplot(1,2,2)

            hsvm = rgb2hsv(colorbar);
            hsvm(:,:,3) = 1;
            hsvm = hsv2rgb(hsvm); 
            imshow(hsvm)  

            
        end
        
        function map = my_turbo(self)
            newturbo = turbo(256);
            newturbo_hsv = rgb2hsv(newturbo);
            newturbo_hsv(:,1) = 0.98*newturbo_hsv(:,1);
            newturbo_hsv(:,2) = 0.9*newturbo_hsv(:,2);
            newturbo_hsv(:,3) = 0.9*newturbo_hsv(:,3);
            map = hsv2rgb(newturbo_hsv);   
            
        end
        
        function map = rich_ametrine(self)
            color(1,:) = [255 255 0];
            color(2,:) = [255 89 0];
            color(3,:) = [208 0 255];
            color(4,:) = [93 0 255];
            color(5,:) = [0 72 255];
            color(6,:) = [0 255 255];
            map = self.rgb_interpolate(color);
            map = flipud(map);
        end
        
        function map = weak_ametrine(self)
            color(1,:) = [252 255 82];
            color(2,:) = [255 82 82];
            color(3,:) = [183 82 255];
            color(4,:) = [125 82 255];
            color(5,:) = [0 82 255];
            color(6,:) = [82 255 255];
            map = self.rgb_interpolate(color);
            
            map = flipud(map);
        end        
        
        function map = rgb_interpolate(self,color)
            [P,R] = size(color);
            marker(1) = 1;
            for pindex = 1:P-1
                marker(pindex+1) = round((255/(P-1))*pindex);
            end
            marker = marker';
            map = interp1(marker,color,(1:255)','linear');
            map = map./255;
            
        end
        
        
        function map = single_ametrine(self)
            
            color(1,:) = [0 0 0];
            color(2,:) = [252 255 82].*0.7;
            color(3,:) = [255 82 82];
            color(4,:) = [183 82 255];
            color(5,:) = [125 82 255];
            color(6,:) = [0 82 255];
            color(7,:) = [82 255 255];
            map = self.rgb_interpolate(color);
            map = flipud(map);
%             
%             hsvmap = rgb2hsv(map);
%             hsvmap(:,3)=(255:-1:1)'./255;
%             
%             map = hsv2rgb(hsvmap);
            
            
        end
        
        
        
        
        function map = ametrine(self)
            map = ametrine;
        end
    end
end

