function exportPSTif(path,fout,N,fwx,dz,logLim,maxRet,opt)
% exports a ps measurement using the pc processing with the provided
% processing parameters.
% The optional parameter provide the possibility to initialize a matlab
% cluster to parallelize the execution

% because the parfor does execute in a non-deterministic way, the output
% files are first written into a temporary folder as individual files, and
% then organized into a stack after compiling all slides.

% input parameters
rS.window = N;
ps.fwx = fwx;
ps.dz = dz;
rs = struct;

% default settings
logF = readLogFile(path);
sliceInd = 1:logF.numImages;

% optional settings
if nargin>7 && ~isempty(opt)
    if isstruct(opt)
        if isfield(opt,'sliceInd')
            sliceInd = opt.sliceInd;
        end
        if isfield(opt,'dzres')
            ps.dzres = opt.dzres;
        end
        if isfield(opt,'bgr')
            rS.bgr = opt.bgr;
            rs.bgr = opt.bgr;
        end
        if isfield(opt,'skipLastPoints')
            rS.skipLastPoints = opt.skipLastPoints;
            rs.skipLastPoints = opt.skipLastPoints;
        end
    end
end


% clear output file
fid = fopen(fout,'w+');
fclose(fid);

tag = sprintf('Format:\tIRD\nInt:\t%d\t%d\nRet:\t0\t%d\nDOP:\t0\t1\nParamsNFwxDz:\t%d\t%d\t%d\n%s',logLim(1),logLim(2),maxRet,N,fwx,dz,datestr(clock));
% loop through files and assemble them in output stack

logOS = logLim(1);
logD = diff(logLim);

for ind = 1:numel(sliceInd)

    [S1,S2] = recstrTom(path,[sliceInd(ind),sliceInd(ind)],rs);
    [Spc1,Spc2] = recstrTom(path,[sliceInd(ind),sliceInd(ind)],rS);
    
    ii = tom2Int(S1,S2);
    pc = PSProcess(Spc1,Spc2,ps);
    
    iiout = uint8(255*(10*log10(ii)-logOS)/logD);
    rrout = uint8(255*pc.rmeancorr/maxRet);
    ddout = uint8(255*pc.dop);
    
    imwrite(iiout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
    imwrite(rrout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
    imwrite(ddout,fout,'tif','WriteMode','append','Compression','none','Description',tag);
    ind
    
end    
    