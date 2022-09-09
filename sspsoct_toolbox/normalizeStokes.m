function [NStokes,DOP] = normalizeStokes( S )
%UNI Summary of this function goes here
%   Detailed explanation goes here
    DOP = sqrt(S(:,1).^2+S(:,2).^2+S(:,3).^2);
    NStokes(:,1) = S(:,1)./DOP;
    NStokes(:,2) = S(:,2)./DOP;
    NStokes(:,3) = S(:,3)./DOP;
end

