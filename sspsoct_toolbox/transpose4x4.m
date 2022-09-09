function MT = transpose4x4(M)
    dim = size(M);
    if dim(1) ~= 16
        error('first dimnesion of the input should be 16');
    end
    M = reshape(M,16,[]);
    MT = M([1 5 9 13 2 6 10 14 3 7 11 15 4 8 12 16],:);
    MT = reshape(MT,dim);
end
