function psi = wavepacket(x,x0,sig,amp)
psi = amp * exp(-(x-x0).^2/(2*sig));
