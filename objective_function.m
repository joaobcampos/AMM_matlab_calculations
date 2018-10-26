function [ fun ] = objective_function( M, translation, rotation )
    %OBJECTIVE_FUNCTION_ROTATION Summary of this function goes here
    %   Detailed explanation goes here
    [r c] = size(M);
    v = zeros(18,1);
    v(10) = rotation(1);
    v(11) = rotation(2);
    v(12) = rotation(3);
    v(13) = rotation(4);
    v(14) = rotation(5);
    v(15) = rotation(6);
    v(16) = rotation(7);
    v(17) = rotation(8);
    v(18) = rotation(9);
    
    v(1) = translation(2) * rotation(3) - translation(3) * rotation(2);
    v(2) = translation(3) * rotation(1) - translation(1) * rotation(3);
    v(3) = translation(1) * rotation(2) - translation(2) * rotation(1);
    v(4) = translation(2) * rotation(6) - translation(3) * rotation(5);
    v(5) = translation(3) * rotation(4) - translation(1) * rotation(6);
    v(6) = translation(1) * rotation(5) - translation(2) * rotation(4);
    v(7) = translation(2) * rotation(9) - translation(3) * rotation(8);
    v(8) = translation(3) * rotation(7) - translation(1) * rotation(9);
    v(9) = translation(1) * rotation(8) - translation(2) * rotation(7);
    %fun = 0;
    aux = zeros(1, 18);
    for i=1:r
        aux = aux + M(i,:);
        %ei = M(i,:) * v;
        %fun = fun + ei * ei;
    end
    ei = aux * v;
    fun = ei * ei;
end

