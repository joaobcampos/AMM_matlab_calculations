%Known values
r13 = sym('r13');r23 = sym('r23');r31 = sym('r31');r32 = sym('r32');

%Unknown values
x1 = sym('x1');x2 = sym('x2'); x3 = sym('x3'); x4 = sym('x4'); x5 = sym('x5');
alpha = sym('alpha')

%% The rotation matrix will have the form
% R = alpha * [x1  x2  r13;
%              x3  x4  r23;
%              r31 r32 x5 ]
ang=rand(3,1)
R = rotx(ang(1)) * roty(ang(2)) * rotz(ang(3))
angle = acos( (R(1,1) + R(2,2) + R(3,3) - 1) / 2);
w = (1/(2 * sin(angle)) ) * [(R(3,2) - R(2,3)); R(1,3) - R(3,1); R(2,1) - R(1,2)];

%Orthogonal constraints
cd1 = ((alpha * x1).^2) + ((alpha * x3).^2) + ((alpha * r31).^2);
cd2 = ((alpha * x2).^2) + ((alpha * x4).^2) + ((alpha * r32).^2);
cd3 = ((alpha * r13).^2) + ((alpha * r23).^2) + ((alpha * x5).^2);
subs(cd1, [alpha, x1 , x3 , r31], [1, R(1,1), R(2,1), R(3,1)])
subs(cd2, [alpha, x2 , x4 , r32], [1, R(1,2), R(2,2), R(3,2)])
subs(cd3, [alpha, r13, r23, x5 ], [1, R(1,3), R(2,3), R(3,3)])

cd1a = ((alpha * x1).^2) + ((alpha * x2).^2) + ((alpha * r13).^2);
cd2a = ((alpha * x3).^2) + ((alpha * x4).^2) + ((alpha * r23).^2);
cd3a = ((alpha * r31).^2) + ((alpha * r32).^2) + ((alpha * x5).^2);
subs(cd1a, [alpha, x1 , x2 , r13], [1, R(1,1), R(1,2), R(1,3)])
subs(cd2a, [alpha, x3 , x4 , r23], [1, R(2,1), R(2,2), R(2,3)])
subs(cd3a, [alpha, r31, r32, x5 ], [1, R(3,1), R(3,2), R(3,3)])



cd4 = (alpha * x1) * (alpha * x2) + (alpha * x3) * (alpha * x4) + (alpha * r31) * (alpha * r32);
cd5 = (alpha * x1) * (alpha * r13) + (alpha * x3) * (alpha * r23) + (alpha * r31) * (alpha * x5);
cd6 = (alpha * x2) * (alpha * r13) + (alpha * x4) * (alpha * r23) + (alpha * r32) * (alpha * x5);
subs(cd4, [alpha, x1 , x2 , x3 , x4 , r31, r32], [1, R(1,1), R(1,2), R(2,1), R(2,2), R(3,1), R(3,2)])
subs(cd5, [alpha, x1 , r13, x3 , r23, r31, x5 ], [1, R(1,1), R(1,3), R(2,1), R(2,3), R(3,1), R(3,3)])
subs(cd6, [alpha, x2 , r13, x4 , r23, r32, x5 ], [1, R(1,2), R(1,3), R(2,2), R(2,3), R(3,2), R(3,3)])

% vectorial product
cd8 =  (alpha * x3) * (alpha * r32) - (alpha * r31) * (alpha * x4) - (alpha * r13);
cd9 = -(alpha * x1) * (alpha * r32) + (alpha * x2) * (alpha * r31) - (alpha * r23);
cd10 =  (alpha * x1) * (alpha * x4) - (alpha * x2) * (alpha * x3) - (alpha * x5);
disp('cd10')
double(subs(cd10, [x1, x2 , x3 , x4, x5, r32, r23, r13, r31, alpha], ...
    [R(1,1), R(1,2), R(2,1), R(2,2), R(3,3), R(3,2), R(2,3), R(1,3), R(3,1), 1.0]))


%axis
cd11 = 3 + 2 * (alpha * x1) + 2 * (alpha * x4) + 2 * (alpha * x5) - ...
    (alpha * x1) * (alpha * x1) - (alpha * x4) * (alpha * x4) - (alpha * x5) * (alpha * x5) - ...
    2 * (alpha * x1) * (alpha * x4) - 2 * (alpha * x1) * (alpha * x5) - 2 * (alpha * x4) * (alpha * x5) - ...
    (alpha * r32 - alpha * r23) * (alpha * r32 - alpha * r23) - ...
    (alpha * r13 - alpha * r31) * (alpha * r13 - alpha * r31) - ...
    (alpha * x3 - alpha * x2) * (alpha * x3 - alpha * x2)
    
subs(cd11, [x1, x2 , x3 , x4, x5, r32, r23, r13, r31, alpha], ...
    [R(1,1), R(1,2), R(2,1), R(2,2), R(3,3), R(3,2), R(2,3), R(1,3), R(3,1), 1.0])

cd12 = ((alpha * x1) + (alpha * x4) ) * (1 - (alpha * x5)) + ...
    (alpha * r32) * (alpha * r23) + (alpha * r13) * (alpha * r31)
disp('cd12')
subs(cd12, [x1, x2 , x3 , x4, x5, r32, r23, r13, r31, alpha], ...
    [R(1,1), R(1,2), R(2,1), R(2,2), R(3,3), R(3,2), R(2,3), R(1,3), R(3,1), 1.0])


%%
% independent relationships
% alpha x4 = f(alpha, alpha x5)
cd13 = (alpha * x4) + (alpha * r13) * (alpha * r31) / (1 - (alpha * x5)) - ...
    (((alpha*r31)*(alpha*r31) - (alpha*r23)*(alpha*r23))/((alpha * r13)*(alpha*r31) + (alpha*r32)*(alpha*r23)))* alpha * x5
disp('cd13: result')
double(subs(cd13, [x1, x2 , x3 , x4, x5, r32, r23, r13, r31, alpha], ...
    [R(1,1), R(1,2), R(2,1), R(2,2), R(3,3), R(3,2), R(2,3), R(1,3), R(3,1), 1.0]))

% alpha x1 = f(alpha, alpha x5)
cd14 = (alpha * x1) + ((alpha * r32) * (alpha * r23) / (1 - alpha * x5)) + ...
    ( ((alpha * r13) * (alpha * r13) - (alpha * r32) * (alpha * r32)) / ...
    ((alpha * r13)*(alpha*r31)+(alpha*r32)*(alpha*r23)) ) * alpha * x5
disp('cd14: result')
double(subs(cd14, [x1, x2 , x3 , x4, x5, r32, r23, r13, r31, alpha], ...
    [R(1,1), R(1,2), R(2,1), R(2,2), R(3,3), R(3,2), R(2,3), R(1,3), R(3,1), 1.0]))

% alpha x2 = f(alpha, alpha x5)
cd15 = (alpha * x2) + ((r13 / r32)-(((r31 * r13 * r13)-(r31 * r32 * r32))/ ...
    ((r32 * r13 * r31)+(r32 * r32 * r23))) ) *  alpha * x5 - (alpha * r31 * alpha * r23)/(1-alpha * x5)
disp('cd15: result')
double(subs(cd15, [x1, x2 , x3 , x4, x5, r32, r23, r13, r31, alpha], ...
    [R(1,1), R(1,2), R(2,1), R(2,2), R(3,3), R(3,2), R(2,3), R(1,3), R(3,1), 1.0]))

% alpha 3 = f(alpha, alpha x5)
cd16 = (alpha * x3) -(alpha * r32) * (alpha * r13) / (1-alpha*x5) + ...
    ((r23/r31)+((r32*r31*r31)-(r32*r23*r23))/((r13*r31*r31)+(r31*r32*r23))) * alpha * x5
disp('cd16: result')
double(subs(cd16, [x1, x2 , x3 , x4, x5, r32, r23, r13, r31, alpha], ...
    [R(1,1), R(1,2), R(2,1), R(2,2), R(3,3), R(3,2), R(2,3), R(1,3), R(3,1), 1.0]))

% alpha x5 = f(alpha)
cd17 = (alpha * r31) * (alpha * r31) + (alpha * r32) * (alpha * r32) + ...
    (alpha * x5) * (alpha * x5) - 1
disp('cd17: result')
double(subs(cd17, [x1, x2 , x3 , x4, x5, r32, r23, r13, r31, alpha], ...
    [R(1,1), R(1,2), R(2,1), R(2,2), R(3,3), R(3,2), R(2,3), R(1,3), R(3,1), 1.0]))

