clear all
importfile('nwav_in.txt');
x=nwav_in(1:end,1);
y=nwav_in(1:end,2);
zr=nwav_in(1:end,3);
zi=nwav_in(1:end,4);
xlin=linspace(min(x),max(x),400);
ylin=linspace(min(y),max(y),400);
[X,Y]=meshgrid(xlin,ylin);
Zr = griddata(x,y,zr,X,Y,'cubic');
Zi = griddata(x,y,zi,X,Y,'cubic');