

function [n,x,sb]=msphj2(n,x,sb)
%This program computes the spherical Bessel functions jn(x)and jn'(x)using subroutine
%                SPHJ
%       Input :  x --- Argument of jn(x)
%                n --- Order of jn(x,n = 0,1,���,� 250)
%       Output:  SJ(n)--- jn(x)
%                DJ(n)--- jn'(x)

nm=[];dj=[];
dj=zeros(1,100+1);
%n=5;
%x=10.0;
[n,x,nm,sj,dj]=sphj(n,x,nm,sj,dj);
end



function [n,x,nm,sj,dj]=sphj(n,x,nm,sj,dj,varargin);
% Compute spherical Bessel functions jn(x)and
%                their derivatives
%       Input :  x --- Argument of jn(x)
%                n --- Order of jn(x,n = 0,1,���)
%       Output:  SJ(n)--- jn(x)
%                DJ(n)--- jn'(x)
%                NM --- Highest order computed

nm=fix(fix(n));
if(abs(x)== 1.0d-100);
for  k=0:n;
sj(k+1)=0.0d0;
dj(k+1)=0.0d0;
end;  k=fix(n)+1;
sj(0+1)=1.0d0;
dj(1+1)=.3333333333333333d0;
return;
end;
sj(0+1)=sin(x)./x;
sj(1+1)=(sj(0+1)-cos(x))./x;
if(n >= 2);
sa=sj(0+1);
sb=sj(1+1);
m=msta1(x,200);
if(m < n);
nm=fix(m);
else;
m=msta2(x,fix(n),15);
end;
f0=0.0d0;
f1=1.0d0-100;
for  k=m:-1:0;
f=(2.0d0.*k+3.0d0).*f1./x-f0;
if(k <= nm)sj(k+1)=f; end;
f0=f1;
f1=f;
end;  k=0-1;
if(abs(sa)> abs(sb))cs=sa./f; end;
if(abs(sa)<= abs(sb))cs=sb./f0; end;
for  k=0:nm;
sj(k+1)=cs.*sj(k+1);
end;  k=fix(nm)+1;
end;
dj(0+1)=(cos(x)-sin(x)./x)./x;
for  k=1:nm;
dj(k+1)=sj(k-1+1)-(k+1.0d0).*sj(k+1)./x;
end;  k=fix(nm)+1;
return;
end



function [msta1Result]=msta1(x,mp,varargin);
% Determine the starting point for backward
%                recurrence such that the magnitude of
%                Jn(x)at that point is about 10^(-MP)
%       Input :  x     --- Argument of Jn(x)
%                MP    --- Value of magnitude
%       Output:  MSTA1 --- Starting point

a0=abs(x);
n0=fix(1.1.*a0)+1;
f0=envj(n0,a0)-fix(mp);
n1=n0+5;
f1=envj(n1,a0)-fix(mp);
for  it=1:20;
nn=n1-(n1-n0)./(1.0d0-f0./f1);
f=envj(nn,a0)-fix(mp);
if(abs(nn-n1)< 1)break; end;
n0=n1;
f0=f1;
n1=nn;
f1=f;
end;
msta1Result=fix(nn);
return;
end



function [msta2Result]=msta2(x,n,mp,varargin);
% Purpose: Determine the starting point for backward
%                recurrence such that all Jn(x)has MP
%                significant digits
%       Input :  x  --- Argument of Jn(x)
%                n  --- Order of Jn(x)
%                MP --- Significant digit
%       Output:  MSTA2 --- Starting point
%
a0=abs(x);
hmp=0.5d0.*fix(mp);
ejn=envj(fix(n),a0);
if(ejn <= hmp);
obj=fix(mp);
n0=fix(1.1.*a0);
else;
obj=hmp+ejn;
n0=fix(n);
end;
f0=envj(n0,a0)-obj;
n1=n0+5;
f1=envj(n1,a0)-obj;
for  it=1:20;
nn=n1-(n1-n0)./(1.0d0-f0./f1);
f=envj(nn,a0)-obj;
if(abs(nn-n1)< 1)break; end;
n0=n1;
f0=f1;
n1=nn;
f1=f;
end;
msta2Result=fix(nn+10);
return;
end
function [envjResult]=envj(n,x,varargin);
envjResult=0.5d0.*log10(6.28d0.*fix(n))-fix(n).*log10(1.36d0.*x./fix(n));
return;
end

