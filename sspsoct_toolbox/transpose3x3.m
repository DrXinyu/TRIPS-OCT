function MT = transpose3x3(M)
    dim = size(M);
    if dim(1) ~= 9
        error('first dimnesion of the input should be 9');
    end
    M = reshape(M,9,[]);
    MT = M;
    MT(1,:) = M(1,:);
    MT(2,:) = M(4,:);
    MT(3,:) = M(7,:);
    MT(4,:) = M(2,:);
    MT(5,:) = M(5,:);
    MT(6,:) = M(8,:);
    MT(7,:) = M(3,:);
    MT(8,:) = M(6,:);
    MT(9,:) = M(9,:);
    MT = reshape(MT,dim);

end
