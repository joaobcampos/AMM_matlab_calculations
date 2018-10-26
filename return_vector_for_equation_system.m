function final_vec = return_vector_for_equation_system(parameters, type)
%RETURN_VECTOR_FOR_EQUATION_SYSTEM Summary of this function goes here
%   Detailed explanation goes here
[r c] = size(parameters);
final_vec_regular = zeros(1,18);
final_vec_central_infinite = zeros(1,5);
if (r==1 && c==36 ) || (r==36 && c==1)
    final_vec_regular(1:3) = parameters(1:3); %e1
    final_vec_regular(4:6) = parameters(7:9); %e2
    final_vec_regular(7:9) = parameters(13:15); %e3
    final_vec_regular(10:12) = parameters(4:6) + parameters(19:21);
    final_vec_regular(13:15) = parameters(10:12) + parameters(25:27);
    final_vec_regular(16:18) = parameters(16:18) + parameters(31:33);
    
    final_vec_central_infinite(1) = parameters(15);
    final_vec_central_infinite(2) = parameters(21);
    final_vec_central_infinite(3) = parameters(27);
    final_vec_central_infinite(4) = parameters(16);
    final_vec_central_infinite(5) = parameters(17);
    if strcmp(type, 'central_infinite')
        final_vec = final_vec_central_infinite;
    else
       final_vec = final_vec_regular;
    end
else
    error('The coeficient vector is not correct (wrong dimensions)')
end

