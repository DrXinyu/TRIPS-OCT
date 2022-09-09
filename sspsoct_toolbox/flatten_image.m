function new_im_seq = flatten_image(image,flatten_curve)
%flatten image
    image = gather(image);
    IS = size(image);
    L = IS(2);
    D = IS(1);
    
    flaten_s = round(polyval(flatten_curve,1:L)+D);
    flaten_s(flaten_s<(D+1)) = 10+D;
    flaten_s(flaten_s>(2*D-10)) = 2*D-10;
    
    flaten_mask_row = repmat(flaten_s,2*D,1);
    flaten_mask_row = flaten_mask_row - repmat((D:-1:-D+1)',1,L);
    flaten_mask_col = repmat((1:L),2*D,1);
    
    im_seq = reshape(image,D,L,[]);
    [~,~,n] = size(im_seq);
    new_im_seq = ones(D*2,L,n);

    for index = 1:n
        newimage = ones(D*3,L).*1e-9;
        newimage(D+1:2*D,:) = im_seq(:,:,index);
        indxp = sub2ind(size(newimage),flaten_mask_row,flaten_mask_col);
        aaa = newimage(indxp);
        new_im_seq(:,:,index) = aaa;  
        
    end
    IS(1) = 2*D;
    new_im_seq = reshape(new_im_seq,IS);
    

end

