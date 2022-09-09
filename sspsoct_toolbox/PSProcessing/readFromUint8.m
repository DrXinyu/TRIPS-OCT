function oaStr = readFromUint8(byteString,Nb,NAlines)
%oaStr = readFromUint8(byteString) reads the data stored in the first bytes
%of the optic axis channel of the pstif format and converts it to the oaStr
%structure.

%Nb : number of bins
%NAlines : number of Alines
%
%Elements always present:
%wcorr [-pi,pi] 3xNb
%errInit [0,1] Nb
%errFinal [0,1] Nb
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

wordStr = double(byteString(1:2:end))*2^8 + double(byteString(2:2:end));

version = wordStr(1);

pointer = 1;

range = [-pi,pi];
count = 3*Nb;
oaStr.wcorr = reshape(wordStr(pointer + (1:count))/2^16*diff(range)+range(1),[3,Nb]);
pointer = pointer + count;

range = [0,1];
count = Nb;
oaStr.errInit = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
pointer = pointer + count;

range = [0,1];
count = Nb;
oaStr.errFinal = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
pointer = pointer + count;

range = [-pi,pi];
count = 3*Nb;
oaStr.rc = reshape(wordStr(pointer + (1:count))/2^16*diff(range)+range(1),[3,Nb]);
pointer = pointer + count;

range = [0,1];
count = Nb;
oaStr.alignErr = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
pointer = pointer + count;

if version == 2 || version == 102
    range = [0,10];
    count = Nb;
    oaStr.alignErrInit = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
    pointer = pointer + count;
end

if version>100%catheter data
    range = [-pi,pi];
    count = NAlines;
    oaStr.tissueAngle = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
    pointer = pointer + count;

    range = [0,pi];
    count = NAlines;
    oaStr.tissueRet = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
    pointer = pointer + count;

    range = [-pi,pi];
    count = 1;
    oaStr.sheathAngle = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
    pointer = pointer + count;

    range = [0,1];
    count = 1;
    oaStr.errSheathAngle = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
    pointer = pointer + count;

    range = [-2*pi,2*pi];
    count = 2*NAlines;
    oaStr.OAball = reshape(wordStr(pointer + (1:count))/2^16*diff(range)+range(1),[2,NAlines]);
    pointer = pointer + count;

    range = [-pi,pi];
    count = 1;
    oaStr.mret = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
    pointer = pointer + count;

    range = [-pi,pi];
    count = 2;
    oaStr.RL1 = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
    pointer = pointer + count;

    range = [0,1];
    count = 1;
    oaStr.err = wordStr(pointer + (1:count))/2^16*diff(range)+range(1);
    pointer = pointer + count;
    
    count = 3*NAlines;
    oaStr.cath = reshape(wordStr(pointer + (1:count)),[3,NAlines]);
    pointer = pointer + count;

    count = NAlines;
    oaStr.tissueSurf = wordStr(pointer + (1:count));
end

