function p1_corr = create_correspondences_for_central_infinite_case(pl_coordinates)
%CREATE_CORRESPONDENCES_FOR_CENTRAL_INFINITE_CASE Summary of this function goes here
%   Detailed explanation goes here

[r c] = size(pl_coordinates);
p1_corr = zeros(r, c);
for i=1:c
    d2x = pl_coordinates(1, i);
    d2y = pl_coordinates(2, i);
    m2z = pl_coordinates(6, i);
    
    m1 = randn(2,1);
    d1z = -(d2x/m2z) * m1(1) - (d2y / m2z) * m1(2);
    
    p1_corr(3, i) = d1z;
    p1_corr(4:5, i) = m1;
end

end

