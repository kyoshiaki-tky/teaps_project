function u = MB_constraint(u)
% describe mass balance constraints
% u: vector, parameter number x 1
%    or matrix, parameter number x number of parameter set
% ex.) GOOD expression: u(3)= 2-u(4)
% ex.) BAD  expression: u(3)+u(4)=2 

u(3)=1;

end
