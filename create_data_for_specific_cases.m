clear all
%%%Central infinite case
% l = (0,0,dz, mx, my,0)'
n_rays = 20

%%Rays in the camera 2 referential
nonzero_elements = 20*rand(3,n_rays);
pl_2 = zeros(6,n_rays);
pl_2(3:5, :) = nonzero_elements;
clear nonzero_elements;

%%Rotation to camera 1 referential
angles = 150* rand(3,1);
R = rotz(angles(3)) * roty(angles(2)) * rotx(angles(1));
t = randn(3,1);
%% Transform 
E = [R zeros(3,3); skew_symmetric_matrix_from_vector(t) * R R];
pl_2_1 = E * pl_2;

%res = pl_2 - [R' zeros(3,3); (skew_symmetric_matrix_from_vector(-R' * t) * R') R'] * pl_2_1
%%Obtain the correspondences in reference 1
p1 = create_correspondences_for_central_infinite_case(pl_2_1);
intersections = obtain_intersection_points(p1,pl_2_1);

%%
% Lines must intersect
LPointC1 = intersections;
LPointC2 = intersections;

figure(1);
hold on;
plot3(intersections(1,:), intersections(2,:), intersections(3,:), 'p', 'color', 'k');
xlabel('-x-'); ylabel('-y-'); zlabel('-z-');
for iter = 1 : n_rays
    points1=[];
    points2=[];
    for lambda = -5:5
        d1 = p1(1:3, iter);
        d1 = d1/norm(d1);
        x1 = intersections(:, iter) + lambda * d1
        points1=[points1 x1];
        
        d2 = pl_2_1(1:3,iter);
        d2 = d2/norm(d2);
        x2 = intersections(:, iter) + lambda * d2
        points2=[points2 x2];
    end
    plot3(points1(1,:), points1(2,:), points1(3,:), ':', 'color','b')
    hold on
    plot3(points2(1,:), points2(2,:), points2(3,:), ':', 'color','r')
    hold on
end
grid on;
hold off;

%% Create equation system

A = []
for i=1:5
    A = [A;return_vector_for_equation_system(kron(pl_2(:,i),p1(:,i)), 'central_infinite')];
end   

n = null(A)
R_sol = zeros(3);
R_sol(3,1) = n(2);
R_sol(3,2) = n(3);
R_sol(1,3) = n(4);
R_sol(2,3) = n(5);

S = [n(2) n(3) 0 0 n(4);
    0 0 n(2) n(3) n(5);
    n(4) 0 n(5) 0 n(2);
    0 n(4) 0 n(5) n(3);
    0 n(5) 0 -n(4) 0;
    -n(5) 0 n(4) 0 0;
    0 0 n(3) -n(2) 0;
    -n(3) n(2) 0 0 0 ];

r = [0; 0; 0; 0; n(2); n(3); n(4); n(5)];

%rest = inv(S) * r;
%R_sol(1,1) = rest(1);
%R_sol(1,2) = rest(2);
%R_sol(2,1) = rest(3);
%R_sol(2,2) = rest(4);
%R_sol(3,3) = rest(5);