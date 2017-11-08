
steps = 10;
d = 1 / (steps - 1);
alpha = 0.5;
% Make The A - matrices
e = ones(steps,1);
A = spdiags([-alpha*e 2*(1+alpha)*e -alpha*e], -1:1, steps, steps);

tril(A);

b = 0:d:1;

f = A \b'

f_test = thomas_solver(A,b)
