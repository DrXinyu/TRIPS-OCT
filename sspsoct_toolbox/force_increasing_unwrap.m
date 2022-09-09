function [ang] = force_increasing_unwrap(ang)
    dif_ang = diff(ang);
    
    mask = zeros(size(dif_ang));
    mask(dif_ang>0) = 0;
    mask(dif_ang<=0.0001) = 2*pi;
    zs = zeros(size(ang));
    ang_com = cumsum([zs(1,:);mask]);
    ang = ang+ang_com;
end

