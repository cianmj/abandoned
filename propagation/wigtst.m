function W2 = wigtst(x,p,func)
N = length(func);
y = x;
Mx = x*ones(1,N);
My = y*ones(N,1);
dxy = abs(x(2)-x(1));
ip = min(max(1,int32((Mx+My)/dxy + N/2)),N);
im = min(max(1,int32((Mx-My)/dxy + N/2)),N);
ff = transpose(func(ip).*func(im))*exp(2i*transpose(y)*p);
W2 = real(ff*x(N));
return
