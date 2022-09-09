function Mout = Jones2MuellerTom(JonesMatIn)
%Mout = Jones2MuellerTom(JonesMatIn) converts the vectorized jones matrix
%JonesMatIn into its corresponding Mueller-Jones matrix Mout.
% The first two dimensions of JonesMatIn are spatial, along the third 
% dimension is the order : t11, t12, t21, t22
% Mout is the linearization of the Mueller matrix, i.e. 
% reshape(Mout(indz,indx,:),[4,4]) yields the Mueller matrix.

Inds = [1,            13,            16,             4;        ...  
   1,            13,           -16,            -4;          ...
  14,            15,             2,             3;          ...
        1i*14,        -1i*15,        0 + 1i*2,        0 - 1i*3;...
   1,           -13,           -16,             4;          ...
   1,           -13,            16,            -4;          ...
 -14,           -15,             2,             3;          ...
        -1i*14,        0 + 1i*15,        0 + 1i*2,        0 - 1i*3;...
  12,             5,             8,             9;          ...
 -12,             5,            -8,             9;          ...
  10,            11,             6,             7;          ...
        1i*10,        -1i*11,        1i*6,        -1i*7;...
        1i*12,        - 1i*5,        -1i*8,        1i*9;...
        -1i*12,       - 1i*5,        1i*8,        1i*9;...
        1i*10,        1i*11,        -1i*6,        -1i*7;...
 -10,            11,             6,            -7];
 
% initialization
KronJ = zeros([size(JonesMatIn,1),size(JonesMatIn,2),16],class(JonesMatIn));

% construct Muller matrices
KronJ(:,:,1) = abs(JonesMatIn(:,:,1)).^2;
KronJ(:,:,2) = JonesMatIn(:,:,1).*conj(JonesMatIn(:,:,3));
KronJ(:,:,3) = JonesMatIn(:,:,3).*conj(JonesMatIn(:,:,1));
KronJ(:,:,4) = JonesMatIn(:,:,3).*conj(JonesMatIn(:,:,3));
KronJ(:,:,5) = JonesMatIn(:,:,1).*conj(JonesMatIn(:,:,2));
KronJ(:,:,6) = JonesMatIn(:,:,1).*conj(JonesMatIn(:,:,4));
KronJ(:,:,7) = JonesMatIn(:,:,3).*conj(JonesMatIn(:,:,2));
KronJ(:,:,8) = JonesMatIn(:,:,3).*conj(JonesMatIn(:,:,4));
KronJ(:,:,9) = JonesMatIn(:,:,2).*conj(JonesMatIn(:,:,1));
KronJ(:,:,10) = JonesMatIn(:,:,2).*conj(JonesMatIn(:,:,3));
KronJ(:,:,11) = JonesMatIn(:,:,4).*conj(JonesMatIn(:,:,1));
KronJ(:,:,12) = JonesMatIn(:,:,4).*conj(JonesMatIn(:,:,3));
KronJ(:,:,13) = JonesMatIn(:,:,2).*conj(JonesMatIn(:,:,2));
KronJ(:,:,14) = JonesMatIn(:,:,2).*conj(JonesMatIn(:,:,4));
KronJ(:,:,15) = JonesMatIn(:,:,4).*conj(JonesMatIn(:,:,2));
KronJ(:,:,16) = JonesMatIn(:,:,4).*conj(JonesMatIn(:,:,4));

Mout = zeros([size(JonesMatIn,1),size(JonesMatIn,2),16],class(JonesMatIn));
for ind = 1:16
    Mout(:,:,ind) = sign(Inds(ind,1))*KronJ(:,:,abs(Inds(ind,1))) + sign(Inds(ind,2))*KronJ(:,:,abs(Inds(ind,2))) + sign(Inds(ind,3))*KronJ(:,:,abs(Inds(ind,3))) + sign(Inds(ind,4))*KronJ(:,:,abs(Inds(ind,4))); 
end
Mout = real(Mout)/2;


