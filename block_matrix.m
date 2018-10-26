function M = block_matrix( LPluckerC1, LPluckerW2, noise1, noise2 )
    %UNTITLED3 Summary of this function goes here
    %   Detailed explanation goes here
    [rP1, cP1] = size(LPluckerC1);
    [rP2, cP2] = size(LPluckerW2);
    [rN1, cN1] = size(LPluckerC1);
    [rN2, cN2] = size(LPluckerW2);
    if (rP1 ~= rP2 | rP1 ~= rN1 | rP2 ~= rN2) 
        error("The number of rows is not consistent");
    end
    if (cP1 ~= cP2 | cP1 ~= cN1 | cP2 ~= cN2) 
        error("The number of columns is not consistent");
    end
    L1 = LPluckerC1 + noise1;
    L2 = LPluckerW2 + noise2;
    M = zeros(cP1, 18);
    
    for i=1:cP1
        l = kron(L2(:,i)', L1(:,i));
        aux = zeros(1,9);
        aux(1:3) = l(1:3);
        aux(4:6) = l(7:9);
        aux(7:9) = l(13:15);
        aux(10:12) = l(4:6) + l(19:21);
        aux(13:15) = l(10:12) + l(25:27);
        aux(16:18) = l(16:18) + l(31:33);
        aux = aux/norm(aux);
        size(aux)
        M(i,:) = aux;
    end
end

