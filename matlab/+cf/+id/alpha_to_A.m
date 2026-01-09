function A = alpha_to_A(alpha)
% Build Ai matrix as in equation (A.8). 
% alpha is 22x1.

a = alpha(:);
A = zeros(8,8);

A(1,:) = [1 0 0 0 0 0 0 0];
A(2,:) = [a(1) 1 0 0 0 0 0 0];
A(3,:) = [a(2) a(6) 1 0 0 0 0 0];
A(4,:) = [0 0 0 1 0 0 0 0];
A(5,:) = [0 0 0 a(12) 1 a(16) 0 0];
A(6,:) = [a(3) a(7) 0 a(13) a(17) 1 0 0];
A(7,:) = [a(4) a(8) a(10) a(14) a(18) a(20) 1 0];
A(8,:) = [a(5) a(9) a(11) a(15) a(19) a(21) a(22) 1];

end

