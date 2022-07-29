function Pp = pump(x,pump,signal,asef,aseb)
%This is a function to calculate the pump power, and this function will be
%called by the RK4 method
global Gamma_p;
global alpha_p;
global sigma_pa; 
global sigma_pe;
x1=x;
pump1=pump;
signal1=signal;
asef1=asef;
aseb1=aseb;
%call the population inversion function to calculate the population of ions
%for given pump power,signal power and ase power
[N1,N2,N3] = PpIn(x1,pump1,signal1,asef1,aseb1);
n1=N1;
n2=N2;
n3=N3;
%The equation can be found in the report
Pp= (-1)*(pump.*Gamma_p).*(sigma_pa*n1 - sigma_pe*n3)- (alpha_p.*pump );
end

