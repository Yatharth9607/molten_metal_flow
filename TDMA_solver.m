% Matrix solver (TDMA)

% parameters required:
% N: Number of unknowns
% A: Coefficient matrix
% b: RHS/constant array prior to calling TDMA

% return output:
% b: solution array after calling TDMA

function [b] = TDMA_solver(N, A, b)
    
    % forward substitution
    A(1,3) = -A(1,3)/A(1,2);
    b(1) = b(1)/A(1,2);
    for i = 2:N
        A(i,3) = -A(i,3)/(A(i,2) + A(i,1)*A(i-1,3));
        b(i) = (b(i) - A(i,1)*b(i-1))/(A(i,2) + A(i,1)*A(i-1,3));
    end
    
    % backward elimination
    for i = N-1:-1:1
        b(i) = A(i,3)*b(i+1) + b(i);
    end

end
    