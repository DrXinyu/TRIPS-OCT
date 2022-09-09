function [c,b,x,tb] = getAlineNum(filename)
    %%% c is the aline number per frame
    %%% b is the bscan number per volume
    %%% x 
    
    pn = 2560*2;
    
    s = dir(filename);         
    filesize = s.bytes;  
    
    if filename(end-1) == 'x' || filename(end) == 'x'
        [filepath,~,~] = fileparts(filename);
        txtname = strcat(filepath,'.txt');
    else
        txtname = strcat(filename,'.txt');        
    end

    

    fileID = fopen(txtname,'r');
    if fileID <0
        c = 500;
        return;
    end
    c = fscanf(fileID, '#AlinesPerFrame@%d');
    if isnumeric(c) && ~(isempty(c))
    else
        c = 500;
    end
    fclose(fileID);
    
    tb = floor(filesize/2/pn/c);
    
    if tb/c ~= round(tb/c)
        b = (round(tb/c))*c;
    else
        b = tb;
    end
    
    x = round(b/c);    
        
        
end

