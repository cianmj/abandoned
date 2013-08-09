n=[];x=[];sj=[];
sj=zeros(1,100+1);
Nxx=40;Nyy=60;Lmax=10;
dxx=40.0/Nxx; dyy=40.0/Nyy;
xx(1)=-20.0;
for ii=1:Nxx-1
    xx(ii+1)=xx(ii)+dxx;
    yy(1)=-20.0;
    for jj=1:Nyy-1
        yy(jj+1)=yy(jj)+dyy;
        r=sqrt(xx(ii+1).^2+yy(jj+1).^2);
        theta=atan(-xx(ii+1)/yy(jj+1));
        pw0(ii,jj)=exp(1i*r*cos(theta));
        pw1(ii,jj)=0.0;
        n=Lmax;x=double(xx(ii+1));
        [n,x,sj]=sphj2(n,x,sj);
        for ll=1:Lmax
            P=legendre(ll,cos(theta));
            pw1(ii,jj)=pw1(ii,jj)+1i.^ll*(2*ll+1)*P(1)*sj(ll);
        end
    end
end
mesh(real(pw0))