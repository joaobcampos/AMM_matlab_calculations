% Gen synthetic data
%This script will generate the symbolic data that will be used in later stages of the program.
%It will confirm that the generated formulas are correct.
% final1 - At this point matrix M is the one that will multiply the vetorialization
% of matrix EgtSideOper. This confirms that the blocks that will be multiplied by each 
% column are the correct ones.
% final2/final3 - At this point we intend to check that the block of matrix M can be put in
% a more compact way. This variable checks that that compacted form is correct.
% final4 - We check the matrices that will multiply by the vetorialization of tr and r will hold.
% At this point, we no longer have EgtSideOper, but the rotation matrix r and tr.
% final5 - Check final4 with real data. We use the function calcER to check that everythong is correct.
% All these test variables must be null.

clear all
clc

NPoints = input('Number of points');
STD_DATA_POINTS = 10;

IPoints = randn( 3, NPoints )*STD_DATA_POINTS;
LDirC1 = randn( 3, NPoints );
LDirC2 = randn( 3, NPoints );
LPointC1 = zeros(size(LDirC1));
LPointC2 = zeros(size(LDirC2));
LPluckerC1 = zeros(6,NPoints);
LPluckerC2 = zeros(6,NPoints);


for iter = 1 : NPoints
    LDirC1(:,iter) = LDirC1(:,iter)/norm(LDirC1(:,iter));
    LDirC2(:,iter) = LDirC2(:,iter)/norm(LDirC2(:,iter));
    LPointC1(:,iter) = randn(1)*STD_DATA_POINTS*LDirC1(:,iter) + IPoints(:,iter);
    LPointC2(:,iter) = randn(1)*STD_DATA_POINTS*LDirC2(:,iter) + IPoints(:,iter);
    
    LPluckerC1(:,iter) = [LDirC1(:,iter); cross(LPointC1(:,iter),LDirC1(:,iter))];
    LPluckerC2(:,iter) = [LDirC2(:,iter); cross(LPointC2(:,iter),LDirC2(:,iter))];
end

%%
% Lines must intersect

figure(1);
hold on;
plot3(IPoints(1,:), IPoints(2,:), IPoints(3,:), 'p', 'color', 'k');
xlabel('-x-'); ylabel('-y-'); zlabel('-z-');

plot3(LPointC1(1,:), LPointC1(2,:), LPointC1(3,:), 'p', 'color', 'b');
plot3(LPointC2(1,:), LPointC2(2,:), LPointC2(3,:), 'p', 'color', 'r');
for iter = 1 : NPoints
    P_ = -5*LDirC1(:,iter) + LPointC1(:,iter);
    P__ = 5*LDirC1(:,iter) + LPointC1(:,iter);
    plot3([P_(1),P__(1)], [P_(2),P__(2)], [P_(3),P__(3)], ':', 'color','b');
    quiver3(LPointC1(1,iter),LPointC1(2,iter),LPointC1(3,iter),LDirC1(1,iter),LDirC1(2,iter),LDirC1(3,iter), 2.0, 'Color', 'b');

    P_ = -5*LDirC2(:,iter) + LPointC2(:,iter);
    P__ = 5*LDirC2(:,iter) + LPointC2(:,iter);
    plot3([P_(1),P__(1)], [P_(2),P__(2)], [P_(3),P__(3)], ':', 'color','r');
    quiver3(LPointC2(1,iter),LPointC2(2,iter),LPointC2(3,iter),LDirC2(1,iter),LDirC2(2,iter),LDirC2(3,iter), 2.0, 'Color', 'r');
end
grid on;
hold off;




%%
% After a rigid transformation


R1 = angle2dcm([randn(1)*180 randn(1)*180 randn(1)*180], 'ZYX');
% R1 = angle2dcm([pi/4 pi/2 pi/4], 'ZYX');
t = randn(3,1)*STD_DATA_POINTS;
sT = [0 -t(3) t(2) ; t(3) 0 -t(1) ; -t(2) t(1) 0 ];
E1 = [R1 zeros(3,3); sT*R1 R1 ];
%E1 = [sT*R R; R zeros(3,3)]

Egt1 = inv(E1);


Rgt = R1';
tgt = -R1'*t;
sTgt = [0 -tgt(3) tgt(2) ; tgt(3) 0 -tgt(1) ; -tgt(2) tgt(1) 0 ];
Egt2 = [ Rgt zeros(3,3); sTgt*Rgt Rgt];



LPluckerW2 = E1 * LPluckerC2;
lL = LPluckerC1(:,1);
lRC2 = LPluckerC2(:,1);
lRW2 = LPluckerW2(:,1);


%Shape of the matrix that makes the side operator go to zero:
% Egt = [(stgT * Rgt) Rgt; Rgt 0]
EgtSideOper = [Egt1(4:6,:); Egt1(1:3,:)];

clearvars -except LPluckerW2 LPluckerC1 EgtSideOper NPoints
