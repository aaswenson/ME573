function f = thomas_solver(A, b)
    
    upper_diag = triu(A);
    diag_mat = diag(A);
    lower_diag = tril(A);
    f = zeros(length(b));
    gamma = zeros(length(b) - 1);
    beta = diag_mat(1);
    f(1) = b(1) / beta;

    for j = 2:length(b)
        
        gamma(j) = lower_diag(j-1) / beta;
        beta = diag_mat(j) - upper_diag(j) * gamma(j);
      
    f(j) = (b(j) - upper_diag(j)*f(j-1)) / beta;
    end
    
    for j=1:(length(b)-1)
        k = length(b) - j;
        f(k) = f(k) - gamma(k+1)*f(k+1);
    end
end