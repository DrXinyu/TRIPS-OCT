function outStr = writeToUint8(input)
%outStr = writeToUint16(input) converts some of the fields of the output
%structure obtained from OAProcessingCath.m or OAProcessingBT.m into uint16
%format to be saved in the first few lines of the optic axis images channel
%in the pstif format.

%Nb : number of bins
%NAlines : number of Alines
%
%Elements always present:
%wcorr [-pi,pi] 3xNb
%errInit [0,1] Nb
%errEff [0,1] Nb % the effective error of the symmetrization step with the
%possibly provided wcorr
%rc [-pi,pi] 3xNb
%alignErr [0,1] Nb
%alignErrInit [0,10] Nb (only from version 2/102)
%
%Elements only present in catheter data:
%tissueAngle [-pi,pi] NAlines
%tissueRet [0,pi] NAlines
%sheathAngle [-pi,pi] 1
%errSheathAngle [0,1] 1
%ball.OA [-2*pi,2*pi] 2xNAlines
%ball.mret [-pi,pi] 1
%ball.RL1 [-pi,pi] 2
%ball.err [0,1] 1
%cath [0,1024], integers, 3xNAlines
%tissueSurf [0,1024], integers 1xNAlines
%
%Total BT: 10xNb + 1
%Total Catheter: 10xNb + 8xNAlines + 7


% Everything is converted into uint16 and appended into a single array.
% Then it will be split into low and high bytes.

wordStr = 2;% Version number 1: BT; 102: Catheter; Versions 1 and 101 did not include the additional field alignErrInit

range = [-pi,pi];
wordStr = cat(1,wordStr,round((input.wcorr(:)-range(1))/diff(range)*2^16));

range = [0,1];
wordStr = cat(1,wordStr,round((input.errInit(:)-range(1))/diff(range)*2^16));

range = [0,1];
wordStr = cat(1,wordStr,round((input.errEff(:)-range(1))/diff(range)*2^16));

if ~isempty(input.rc)
    range = [-pi,pi];
    wordStr = cat(1,wordStr,round((input.rc(:)-range(1))/diff(range)*2^16));
else
    wordStr = cat(1,wordStr,ones(3*numel(input.errEff),1)/2*2^16);
end

range = [0,1];
wordStr = cat(1,wordStr,round((input.alignErr(:)-range(1))/diff(range)*2^16));

if ~isempty(input.alignErrInit)
    range = [0,10];
    wordStr = cat(1,wordStr,round((input.alignErrInit(:)-range(1))/diff(range)*2^16));
else
    wordStr = cat(1,wordStr,zeros(numel(input.errEff),1));
end
if isfield(input,'tissueAngle')
    wordStr(1) = 102;
    range = [-pi,pi];
    wordStr = cat(1,wordStr,round((input.tissueAngle(:)-range(1))/diff(range)*2^16));

    range = [0,pi];
    wordStr = cat(1,wordStr,round((input.tissueRet(:)-range(1))/diff(range)*2^16));

    range = [-pi,pi];
    wordStr = cat(1,wordStr,round((input.sheathAngle(:)-range(1))/diff(range)*2^16));

    range = [0,1];
    wordStr = cat(1,wordStr,round((input.errSheathAngle(:)-range(1))/diff(range)*2^16));

    range = [-2*pi,2*pi];
    temp = input.ball.OA(1:2,:);
    wordStr = cat(1,wordStr,round((temp(:)-range(1))/diff(range)*2^16));

    range = [-pi,pi];
    wordStr = cat(1,wordStr,round((input.ball.mret(:)-range(1))/diff(range)*2^16));

    range = [-pi,pi];
    wordStr = cat(1,wordStr,round((input.ball.RL1(:)-range(1))/diff(range)*2^16));

    range = [0,1];
    wordStr = cat(1,wordStr,round((input.ball.err(:)-range(1))/diff(range)*2^16));

    wordStr = cat(1,wordStr,round(input.cath(:)));

    wordStr = cat(1,wordStr,round(input.tissueSurf(:)));
end

outStr = uint8(cat(1,floor(wordStr/(2^8))',rem(wordStr,2^8)'));
outStr = outStr(:);

%test = double(outStr(1:numel(outStr)/2))*2^8 + double(outStr(numel(outStr)/2 + 1:end));

