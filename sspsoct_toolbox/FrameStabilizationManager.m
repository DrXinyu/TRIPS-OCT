classdef FrameStabilizationManager < handle
    %Stokes frame process
    properties 
        frame_index
        stabilized_window
        last_frame
        shifts
        
    end
    
    methods
        
        function self = FrameStabilizationManager(stabilized_window)
            
            self.stabilized_window = stabilized_window;
            self.frame_index = 0;
            
        end
        
        
        function [signature_frame_s,following_frames_s] = stabilize(self,signature_frame,following_frames)
             
            self.frame_index = self.frame_index+1;
            if self.frame_index == 1
                self.shifts(self.frame_index) = 0;
                move = 0;
            else
                [value,shiftframe]= max(self.frame_correlation(signature_frame,self.last_frame),[],2);
                level = quantile(value,0.75);
                shiftframe = gather(shiftframe);
                self.shifts(self.frame_index) = (mean(shiftframe(value>level)))-11;
                move = 0;
            end
            
            
            
            self.last_frame = signature_frame;
            [signature_frame_s,following_frames_s] = self.action(move,signature_frame,following_frames);
            

        end
        
        
        function [signature_frame_s,following_frames_s] = action(self,move,signature_frame,following_frames)
            signature_frame_s = 0;
            following_frames_s = 0;
        end
        
        function xc = frame_correlation(self,frame1,frame2)
            frame1 = frame1 - mean(frame1);
            frame2 = frame2 - mean(frame2);
            frame1array = repmat(frame1,1,1,21);
            for findex = 1:21
                frame1array(:,:,findex) = circshift(frame1array(:,:,findex),findex-11,1);
            end
            xc = squeeze(sum(abs(frame1array+frame2),1));
        end
        
        
    end
end

