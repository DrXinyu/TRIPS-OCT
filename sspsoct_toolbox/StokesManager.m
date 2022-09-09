classdef RunningAveragingManager < handle
    %Stokes frame process
    properties
        GPU
        win 
        intwin
        total_len
        Stokes_frames
        Intensity_frames
        pointer  
    end
    
    methods
        function self = RunningAveragingManager(frame_num)
            self.total_len = frame_num;
            %self.GPU = GPU;
            self.pointer = 1;
            self.win = gausswin(frame_num,1)./sum(gausswin(frame_num,1));   
            self.intwin = gausswin(frame_num,2)./sum(gausswin(frame_num,2));  
            clear self.Stokes_frames self.Intensity_frames
        end
        
        
        function self = push(self,Stokes_frame,Intensity_frame)

            self.Stokes_frames{self.pointer} = gather(Stokes_frame);
            self.Intensity_frames{self.pointer} = gather(Intensity_frame);
            self.pointer = self.pointer +1;
            if self.pointer == self.total_len +1
                self.pointer = 1;
            end

        end
        
        function [average_Stokes,middle_intensity] = pop(self)
            middle = mod(self.pointer-round(self.total_len/2)-1,self.total_len)+1;
            wind = circshift(self.win,self.pointer-self.total_len-1);
            intwind = circshift(self.intwin,self.pointer-self.total_len-1);
            average_Stokes = self.Stokes_frames{1}.*wind(1);
            average_Inten = self.Intensity_frames{1}.*intwind(1);
            for index = 2:length(self.Stokes_frames)
                average_Stokes = average_Stokes + self.Stokes_frames{index}.*wind(index);
                average_Inten = average_Inten + self.Intensity_frames{index}.*intwind(index);
            end
%             average_Stokes = average_Stokes;
%             average_Inten = average_Inten;
            middle_intensity = average_Inten;
%             if middle > length(self.Intensity_frames)
%                 middle_intensity = self.Intensity_frames{end};
%                 
%             else
%                 middle_intensity = self.Intensity_frames{middle};                    
%             end
        end   
        
    end
end

