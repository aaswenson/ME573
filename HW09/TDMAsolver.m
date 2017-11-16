function x = TDMAsolver(A,b_mat)
%a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
n = length(b_mat); % n is the number of rows
u = triu(A); d = diag(A); l = tril(A);
% Modify the first-row coefficients
l(1) = l(1) / d(1);    % Division by zero risk.
b_mat(1) = b_mat(1) / d(1);   
 
for i = 2:n-1
    temp = d(i) - u(i) * l(i-1);
    l(i) = l(i) / temp;
    b_mat(i) = (b_mat(i) - u(i) * b_mat(i-1))/temp;
end
 
b_mat(n) = (b_mat(n) - u(n) * b_mat(n-1))/( d(n) - u(n) * l(n-1));
 
% Now back substitute.
x(n) = b_mat(n);
for i = n-1:-1:1
    x(i) = b_mat(i) - l(i) * x(i + 1);
end