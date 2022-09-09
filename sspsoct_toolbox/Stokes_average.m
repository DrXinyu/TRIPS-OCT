function [Sf,dop] = Stokes_average(S,fwx,fwz)

    dim = size(S);
    sdim = length(dim);
    
    I = sqrt(dot(S,S,sdim));
    kernel = (gausswin(fwx,1)./sum(gausswin(fwx,1)))';
    
    Sf = imfilter(S,kernel,'replicate');
    If = imfilter(I,kernel,'replicate');

    kernel = gausswin(fwz,1)./sum(gausswin(fwz,1));
    
    Sf = imfilter(Sf,kernel,'replicate');
    If = imfilter(If,kernel,'replicate');

    % normalization
    QUVf = sqrt(dot(Sf,Sf,sdim));
    dop = QUVf./If;


end

