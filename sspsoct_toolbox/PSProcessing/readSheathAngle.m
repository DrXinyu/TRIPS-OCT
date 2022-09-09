function sheathAngle = readSheathAngle(psStack)

sheathAngle = [];
if ~psStack.stackCheck
    return;
end

%check version    
Nchannels = numel(psStack.format);
ver = psStack.blockData((psStack.NAlines*psStack.Nz)*(Nchannels-1) + (2));

if ver<100
    return;
end
    
% assuming 9 spectral bins
offset = -1;
if ver == 101 
    offset =  2*(9*9 + 2*psStack.NAlines + 1) + 1;
elseif ver == 102
    offset =  2*(10*9 + 2*psStack.NAlines + 1) + 1;
end
if offset>0
  vec1 = psStack.blockData((psStack.NAlines*psStack.Nz)*(Nchannels-1) + offset + (0:psStack.NSlices-1)*(psStack.NAlines*psStack.Nz)*Nchannels);
  vec2 = psStack.blockData(1 + (psStack.NAlines*psStack.Nz)*(Nchannels-1) + offset + (0:psStack.NSlices-1)*(psStack.NAlines*psStack.Nz)*Nchannels);
  sheathAngle = (double(vec1)*2^8 + double(vec2))/2^16*2*pi-pi;
end