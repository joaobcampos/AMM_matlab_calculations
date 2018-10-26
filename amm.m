% First step: Create the data
create_data
%chama 1 segundo script que gera informação

% Create block matrix
use_noise = input('Is noise to be used?');
noise1 = zeros(6, NPoints);
noise2 = zeros(6, NPoints);
std_var = 0.01
if use_noise ~= 0
    noise1 = std_var * randn(6, NPoints);
    noise2 = std_var * randn(6, NPoints);
end
M = block_matrix( LPluckerC1, LPluckerW2, noise1, noise2 )
%Solve problem
iterations = input('Iterations')
[ rotation, translation ] = amm_solver( iterations, M )

diffR = reshape(rotation, 3, 3) - EgtSideOper(4:6,1:3)
skew_t = EgtSideOper(4:6,1:3) * EgtSideOper(1:3,1:3)';