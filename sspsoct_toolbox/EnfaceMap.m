classdef EnfaceMap < handle

    properties
        map
        frame_index
        map_name
    end
    
    methods
        
        function self = EnfaceMap(map_name)
            
            self.frame_index = 0;
            self.map = [];
            self.map_name = map_name;
            
        end
        
        function add_a_line(self,line)
            
            self.frame_index = self.frame_index+1;
            %ret = medfilt2(ret,[3 3]);
            %nanmask = self.retMask(ret,retina_mask);
            self.map = cat(1,self.map,gather(line));
            
        end

        function map = savemap(self,output_dir,filename)
            
            map = (self.map);
            filename = strcat(filename,self.map_name,'_enface.mat');
            filename = fullfile(output_dir,filename);
            
            save(filename,'map');

        end
        
        
    end
end

