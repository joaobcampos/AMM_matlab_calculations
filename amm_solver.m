function [ rotation, translation ] = amm_solver( iterations, M )
    %AMM_SOLVER Summary of this function goes here
    %   Detailed explanation goes here
    %Consider the first state for both the translation and rotation:
    translation = ones(3,1);
    rotation = reshape(eye(3), 9, 1);
    for i = 1:iterations
        % Minimize the of according to the rotation, knowing the translation
        fun = @(x)objective_function( M, translation, x);
        nonlcon = @orthogonal_restriction;
        options = optimoptions('fmincon', 'Display', 'final','ConstraintTolerance',1e-9, ...
        'FunctionTolerance',1e-9);
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        lb = [];
        ub = [];
        rotation = fmincon(fun,rotation,A,b,Aeq,beq,lb,ub,nonlcon, options);
        % Minimize the of according to the translation, knowing the rotation
        fun = @(x)objective_function( M, x, rotation);
        translation = fminunc(fun,translation);
    end
    rotation = reshape(rotation, 3, 3);
end

