matrix_size = input('Enter size of matrix: '); % rozmiar macierzy, która będzie wylosowana
A = generateSPDmatrix(matrix_size) % dana macierz, domyślne losowana (losuję macierz symetryczną, dodatnie określoną)
b = randn(matrix_size, 1) % wektor danych wyrazów wolnych, domyślnie losowany
X = zeros(matrix_size, 1) % wektor początkowy, domyślnie wektor zerowy
e = 1e-7; % dokładność

jac=Jacobi(A,b,X,e);
disp('ttttt')
disp('ttttt')
disp('ttttt')
sei=Seidel(A,b,X,e);
output = [num2str(jac),' - jacobi method ; ',num2str(sei),' - saidel method'];
disp(output)

function A = generateSPDmatrix(n)
% Generate a dense n x n symmetric, positive definite matrix

A = rand(n,n); % generate a random n x n matrix

% construct a symmetric matrix using either
A = 0.5*(A+A');
% The first is significantly faster: O(n^2) compared to O(n^3)

% since A(i,j) < 1 by construction and a symmetric diagonally dominant matrix
%   is symmetric positive definite, which can be ensured by adding nI
A = A + n*eye(n);

end

% iteruje iteracją prostą, dopóki nie będzie osiągnięta dokładność
function j = Jacobi(A,b,X,e)

N = length(b);
P = X+1;
j=1;
while abs((X-P))>e
    P = X;
    for i=1:N
        X(i)= (b(i)/A(i,i)) - ((A(i,[1:i-1,i+1:N])*P([1:i-1,i+1:N])))/A(i,i);
    end
    
    X
    fprintf('**Iteration number %d**\n',j)
    j = j+1;
end
end

% iteruje iteracją Seidla, dopóki nie będzie osiągnięta dokładność
function j = Seidel(A,b,P,e)

N = length(b);
X = zeros(N,1);
Y = X+1;
j=1;
while abs((Y-X))>e
    Y = X;
    for i=1:N
        X(i)= (b(i)/A(i,i)) - ((A(i,[1:i-1,i+1:N])*P([1:i-1,i+1:N])))/A(i,i);
        P(i)= X(i);
    end
    
    X
    fprintf('**Iteration number %d**\n',j)
    j = j+1;
end
end