function [angle,err,angStd] = estimateOrientation(OAtissue,kind)
% estimateOrientation(OAtissue) estimates the least squares orientation of
% the OAtissue. The tissue signal is aligned with the Q-axis.
% kind == 1 uses mean, kind == 2 median

if nargin<2 || kind == 1
    meanfun = @(x)mean(x);
else
    meanfun = @(x)median(x);
end    

angle = atan2(meanfun(2*OAtissue(1,:).*OAtissue(2,:)),meanfun(OAtissue(2,:).^2-OAtissue(1,:).^2))/2;
% pi/2 ambiguity
err1 = meanfun((sin(angle)*OAtissue(1,:) + cos(angle)*OAtissue(2,:)).^2);
err2 = meanfun((sin(angle+pi/2)*OAtissue(1,:) + cos(angle+pi/2)*OAtissue(2,:)).^2);
if err1<err2
    err = err1;
else
    err = err2;
    angle = angle + pi/2;
end
% evaluate center of weight of projection onto this direction to figure out
% the sense of sheathAngle
if mean([cos(angle),-sin(angle)]*OAtissue(1:2,:))<0
    angle = mod(angle+pi+pi,2*pi)-pi;
end

