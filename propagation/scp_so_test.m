%% Finite Square Well: Split-Operator Method
%% 
clear all; close all;
a = 48;  %% Length
M = 1/2; %% Mass
N = 512;
x = linspace(-a/2,a/2,N); x = x';
k = N*linspace(-1/2,1/2,N); k = k';
dt = 1e-3; %% Time step
%%-------------------------------------------------------------------------
%%
%%% Potential
%V0 = 200;
%V = zeros(length(x),1) - V0; % 1*((2*x).^2 - (0.6*a)^2);  %  
b = a/16;
%V(x<-b) = 0;
%V(x>+b) = 0;
kk = 10;
V = 1/2 * kk * x.^2;
%p = -diff(V,1);
p = -gradient(V);
V(1:5) = V(1:5)-1i*1e3;  		%% Absorption at simulation boundary
V(end-5:end) = V(end-5:end)-1i*1e3;	%% Absorption at simulation boundary
%%%
%%-------------------------------------------------------------------------
%%
%%% initial wave packet
%% Gaussian
Phi0 = exp(-(.5*(x+5)).^2);
%
WG = mywigner(x,p,Phi0) * a/N;
figure(2);set(gcf,'position',[37 208 538 532]);
mesh(x,p,WG);hold off; 
pause(2);
%return
GK = fftshift(exp(-(1i*p.^2*dt/(4*M)))); %% dt/2 kinetic energy propagator
GK2 = fftshift(exp(-(1i*p.^2*dt/(2*M)))); %% dt kinetic energy propagator
GV = exp(-1i*dt*V); %% Potential spatial interaction

iPhi = fft(Phi0);
Phi = ifft(iPhi.*GK);
Phi = GV.*Phi;
NPt = 40;
Pt = zeros(1,NPt);

T = dt*NPt;
t = linspace(0,T,NPt);
%wt = (1-cos(2*pi*t/length(t)));
%uns = 0;
for nrn = 1:NPt
    iPhi = fft(Phi);
    Phi = ifft(iPhi.*GK2);
    if mod(nrn,1) == 0
        figure(1);
        subplot(3,1,1);plot(x, real(Phi),'r'); 
        title(['\Phi_x  t=',num2str(t(nrn))]);
        xlabel('x');ylabel('Amp');axis tight;
        subplot(3,1,2);plot(x, abs(Phi),'r'); 
        title(['\Phi_x  t=',num2str(t(nrn))]);
        xlabel('x');ylabel('Amp');axis tight;
        subplot(3,1,3);plot(k,fftshift(real(iPhi)),'r-');
        title('\Phi_k'); 
        xlabel('k');ylabel('k-space Amp');axis tight;
        pause(0.2);
    end
    %%
    Phi = GV.*Phi;    
end

return

iPhi = fft(Phi);
Phi = ifft(iPhi.*GK);
%%
estep = 1;  %% Sampling period
Po = Pt(1:estep:length(Pt));
T = dt*NPt;
t = linspace(0,T,length(Po));
E = (1/dt)*(linspace(-pi,pi,length(Po)));
%%
Pe = fftshift(fft(((1-cos(2*pi*t/T)).*Po/T)));
%%
figure(2);subplot(2,1,1);plot(t,real(Po));
title('Correlation Function ');xlabel('Time');
figure(2);subplot(2,1,2);plot(E,log(fliplr(abs(Pe))),'r');
title('Energy Spectrum');xlabel('Energy');ylabel('Power');
axis([-210 0 -17 5]);
pause(1);

%%-------------------------------------------------------------------------
%% Analytic method: For Even Solutions (Even Wave functions)
%%
z0 = b*sqrt(2*M*V0);
z = 0:0.01:20*pi;
y1 = tan(z);
y2 = sqrt((z0./z).^2 - 1);
figure(3);subplot(2,1,1);plot(z,y1,z,y2);
hold on;
plot(z,0*z,'r');
axis([0 45 0 35]);
title('tan(z)  =  [(z_0/z)^2 - 1]^{1/2}');
crss_n = [1.5 4.5 7.6 10.8 13.83 16.9 20.0 23.0 26.1 29.1 32.2 35.2 38.2 41.1]; 
%% ^-- get these values by looking at the graph (approx)
g =  inline('tan(z) - sqrt((z0/z).^2 - 1)','z','z0');
for nrn = 1:14
zn(nrn) = fzero(@(z) g(z,z0),crss_n(nrn));
end
figure(3);subplot(2,1,1);hold on;plot(zn,tan(zn),'rx');
q = zn/b;
Em = ((q.^2)/(2*M))-V0;
%%
for nrn = 1:length(Em),
    figure(3);subplot(2,1,2);hold on; 
    plot([Em(nrn),Em(nrn)],[-17,6]); 
end
%%
figure(3);subplot(2,1,2);
plot(E,log(fliplr(abs(Pe))),'r');hold on;
title('Energy Spectrum (Blue: Even solutions)');
xlabel('Energy');ylabel('Power');
axis([-210 0 -17 5]);
%%
%%-------------------------------------------------------------------------

