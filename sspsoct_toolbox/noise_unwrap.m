function [noise] = noise_unwrap(phase)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    diff_p = diff(phase);
    mask_neg = zeros(size(diff_p));
    mask_neg(diff_p < -2) = pi;
    mask_pos = zeros(size(diff_p));
    mask_pos(diff_p > 2) = -pi;
    
    mask = mask_neg+mask_pos;
    com = zeros(size(phase));
    com(2:end) = cumsum(mask);
    
    noise = phase+com;
end

