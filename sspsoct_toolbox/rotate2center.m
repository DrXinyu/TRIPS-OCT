function new_angles = rotate_unwrap(angles)

rotator = -pi:0.1:pi;
std_angle = ones(size(rotator));


    for r_index = 1:length(rotator)
        
        rot_angles = angles - rotator(r_index);
        new_angles = rot_angles;
        new_angles(rot_angles>pi) = rot_angles(rot_angles>pi)-2*pi;
        new_angles(rot_angles<-pi) = rot_angles(rot_angles<-pi)+2*pi;
        std_angle(r_index) = gather(std(new_angles));
        %mean_angle(r_index) = gather(mean(new_angles));
        
    end
    
    [~,m] = min(std_angle);
    %[~,m] = min(abs(mean_angle));
    rot_angles = angles - rotator(m);
    new_angles = rot_angles;
    new_angles(rot_angles>pi) = rot_angles(rot_angles>pi)-2*pi;
    new_angles(rot_angles<-pi) = rot_angles(rot_angles<-pi)+2*pi;
    
    new_angles = new_angles + rotator(m);
    
end

