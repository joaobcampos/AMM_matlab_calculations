function matrix = skew_symmetric_matrix_from_vector(t)
%SKEW_SYMMETRIC_MATRIX_FROM_VECTOR Calculates skew symmetric matrix from
% the translation
matrix = [ 0 -t(3) t(2); t(3) 0 -t(1); -t(2) t(1) 0];
end

