function [average] = galvoAverage(frame)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    c = length(frame(1,:,:));
    frame = circshift(frame,c/4-60,2);

    framefscan = frame(:,1:c/2,:);
    framebscan = fliplr(frame(:,c/2+1:c,:));
    average = (framefscan+framebscan)/2;
    
end

