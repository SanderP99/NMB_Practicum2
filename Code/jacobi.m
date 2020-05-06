function [V,D]= jacobi(A,tol)

% Aanname dat A vierkant is
n = size(A,1);

upperTriangularIndexList = nonzeros(triu(reshape(1:numel(A),size(A)),1));
k = upperTriangularIndexList';
err = sqrt(sum(A(k).^2));
V = eye(n);
row = 1;
col = 2;
while err > tol 
    a = A(row, row);
    b = A(col, col);
    d = A(row, col);
    
    % Theta bepalen
    if a == b
        theta = pi / 4;
    else
        theta = 0.5 * atan(2 * d / (b - a));
    end
    
    % Jacobi matrix opstellen
    J = eye(n);
    c = cos(theta);
    s = sin(theta);
    J(row, row) = c;
    J(col, col) = c;
    J(row, col) = s;
    J(col, row) = -s;
    
    % Iteratiestap
    A = J' * A * J;
    V = V * J;
    
    % Fout vinden
    err = sqrt(sum(A(k).^2));
    
    % row en col ophogen
    if row ~= n - 1
        if col == n
            row = row + 1;
            col = row + 1;
        else
            col = col + 1;
        end
    else
        if col == n
            row = 1;
            col = 2;
        else
            col = col + 1;
        end
    end
end
D = A;
