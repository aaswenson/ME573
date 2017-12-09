function f = thomas(A, b)

n = length(b);

for k = 1:n-1
    i=k+1;
    I = A(i,k) / A(k, k);
    
    for j = k:k+1
        A(i,j) = A(i,j) - I*A(k,j);
    end
    
    b(i) = b(i) - I*b(k);
end

for k=n:-1:1
    f(k) = b(k);
    for j=k+1:min(n,k+1)
       f(k) = f(k) - A(k,j)*f(j);
    end
    
    f(k) = f(k) / A(k,k);

end