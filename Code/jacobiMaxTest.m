function [V,D]= jacobiMaxTest(A,tol)

% Aanname dat A vierkant is
n = size(A,1);
C = A;

upperTriangularIndexList = nonzeros(triu( reshape(1:numel(A), size(A)) ,1));
k = upperTriangularIndexList';
error = sqrt(sum(A(k).^2));
errorList = [error];

V = eye(n);
row = 1;
col = 2;
while error > tol 
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
    
    A = J' * A * J;
    V = V * J;
    
    % Fout vinden
    k = upperTriangularIndexList';
    error = sqrt(sum(A(k).^2));
    
    
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
            errorList = [errorList, error];
        else
            col = col + 1;
        end
    end
end
errorList = [errorList, error];
D = A;
figure('Renderer', 'painters', 'Position', [10 10 900 600], 'PaperPositionMode', 'auto');
set(gca, 'Units', 'normalized','FontUnits','points','FontWeight','normal','FontSize',14,'FontName','Times')

semilogy(0:1:8,errorList);
hold on
plot(size(upperTriangularIndexList,1):size(upperTriangularIndexList,1):size(errorList,2),errorList(size(upperTriangularIndexList,1):size(upperTriangularIndexList,1):end), 'o','MarkerSize',10, 'Color', 'red');
plot(size(errorList,2)-1,errorList(end), 'x','MarkerSize',10, 'Color', 'red');

hold off
xlabel('$n$', 'interpreter', 'latex')
ylabel('Frobenius norm $\|R_A\|_F$', 'interpreter', 'latex')
print  ('-r500', '-depsc2', 'fout-jacobiSweep.eps')