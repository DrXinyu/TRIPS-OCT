classdef FileReader < handle

    properties
        
        buffersize
        data_file
        fid
        buffer_position
        file_bytes
        buffer
        
    end
    
    methods
        
        function self = FileReader(buffersize)
            
            self.buffersize = buffersize;
            self.data_file = 'null';
            self.buffer_position = [0 0];
            
        end
        
        function self = setBufferSize(self,buffersize)
            
            self.buffersize = buffersize;
            
        end        
        
        
        
        
        
        function fopen(self,data_file)
            if strcmp(data_file,self.data_file)         
            else
                if self.fid > 0
                    fclose(self.fid);
                    self.fid = -1;
                end
                file = dir(data_file);
                self.file_bytes = file.bytes;
                self.fid = fopen(data_file,'r', 'b');
                if self.fid < 0
                    disp('file open failed');
                    return
                end
                
                self.data_file = data_file;
            end
        end
        
        function data = fread(self,pointer,block)
            
            if self.is_pointer_in_buffer(pointer,block)
                
                data = self.buffer(self.pointer2buffer_index(pointer):self.pointer2buffer_index(pointer)+block(1)*block(2)-1);
                data = reshape(data,block);
                
            else
                
                fseek(self.fid, pointer, 'bof');              
                self.buffer = fread(self.fid, self.buffersize, 'uint16','n');
                
                self.buffer_position = [pointer, size(self.buffer)*2+pointer];
                
                data = self.buffer(self.pointer2buffer_index(pointer):self.pointer2buffer_index(pointer)+block(1)*block(2)-1);
                data = reshape(data,block);                
            end
           
        end        
        
        function t = is_pointer_in_buffer(self,pointer,block)
            
            if (pointer >= self.buffer_position(1)) && (pointer + block(1)*block(2)*2 < self.buffer_position(2)-1)
                t = 1;
            else
                t = 0;
            end
        end
        
        
        function bindex = pointer2buffer_index(self,pointer)
            
            bindex = (pointer-self.buffer_position(1))/2+1;
            
        end
            
            
     
        
        
        
        
    end
end

