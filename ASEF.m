function PASEF = ASEF(x,pump,signal,asef,aseb)
%This is a function to calculate the forward ASE power, and this function will be
%called by the RK4 method
global Gamma_s;
global alpha_s;
global sigma_se;
global sigma_sa;
global hc;%planck constant 
global fs;  
global delta_v;%bandwidth
delta_v =3.1e+12;%GHz
x1=x;
pump1=pump;
signal1=signal;
asef1=asef;
aseb1=aseb;
%call the population inversion function to calculate the population of ions
%for given pump power,signal power and ase power
[N1,N2,N3] = PpIn(x1,pump1,signal1,asef1,aseb1);
PASEF= (asef1.*Gamma_s).*(sigma_se*N2 - sigma_sa*N1)+2.*sigma_se*N2*Gamma_s*hc*fs.*delta_v -(alpha_s.*asef1 );
end

