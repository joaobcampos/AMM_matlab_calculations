function intersections = obtain_intersection_points(pl1,pl2)
%OBTAIN_INTERSECTION_POINTS Summary of this function goes here
%   Detailed explanation goes here
[r1 c1] = size(pl1);
[r2 c2] = size(pl2);
if r1~=r2 || c1 ~= c2
    error('Dimensions mismatcj in obtain intersection points')
end

intersections = zeros(3, c1);
for i=1:c1
    sk1 = skew_symmetric_matrix_from_vector(pl1(1:3, i));
    m1  = pl1(4:6, i);sk1 = skew_symmetric_matrix_from_vector(pl1(1:3, i));
    
    sk2 = skew_symmetric_matrix_from_vector(pl2(1:3, i));
    m2  = pl2(4:6, i);
    
    %%Create full rank matrix A
    A = zeros(3,3);
    A(1:2,:) = sk1(1:2,:);
    A(3,:)   = sk2(1,:);
    v = zeros(3,1);
    v(1:2,1) = m1(1:2,1);
    v(3,1) = m2(1,1);
    intersections(:,i) = inv(A) * v;
end

end
