function W2 = wigtst3(x,p,func)
N = length(func);
y = x;
gg = zeros(N,N);
for i = 1:N,
    for j = 1:N,
        for k = 1:N-1,
            fk=func(min(N,i+k))*func(max(1,i-k))*exp(2i*p(j)*x(k));
            fk1=func(min(N,i+k+1))*func(max(1,i-k-1))*exp(2i*p(j)*x(k+1));
            gg(i,j) = gg(i,j) + abs(x(2)-x(1))*(fk + fk1)/2.;
        end
    end
end
W2 = abs(gg/pi);
return
