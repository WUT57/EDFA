function [Gain,Gainfinal,X1,PP1,SP1] = func(FiberLength,pump_initial)
%This is the main alogrithm function that use RK4 as iteration method
format long
global Nt;%Er ion concnetration
global lambda_p;%pump wavelength
global lambda_s;%signal
global sigma_pa;%pump absorption cross section
global sigma_pe;%pump emission--
global sigma_se;%signal emission--
global sigma_sa;%signal absoption--
global A_21;%emission rate 
global AR;%fiber cros section area 
global Gamma_s;%signal overpal ratio 
global Gamma_p;%pump-- 
global hc;%planck 
global c_light; 
global fp;%pump freq
global fs;%signal freq
global alpha_s;%Background loss for signal 0.25dB/km
global alpha_p;%background loss for pump
%now start define parameter values
Nt= 6*(10^24);%Er ion density, notice that 100ppm is approximatly 5.7e1024
Gamma_s =0.4; 
Gamma_p =0.4; 
lambda_p =1480*(10^-9);
lambda_s =1550*(10^-9);
sigma_pa =0.75*(10^-25);
sigma_pe =0.19*(10^-25);
sigma_se =3.8*(10^-25);
sigma_sa =2.4*(10^-25);
hc=6.626*(10^-34);
A_21 =100; 
AR =1.25*(10^-11);
c_light=3.0*(10^8);
fp =c_light/lambda_p;  
fs =c_light/lambda_s; 
alpha_s=5.756e-5;%bakcgroud loss
alpha_p=5.756e-5;
%background loss, about 0.25dB/km, this value showed here is
%in unit of /m, the conversion can be found at:
%https://electronics.stackexchange.com/questions/204331/convert-db-km-to-km
n=200;% step number
h=FiberLength/n;%step size, 
x_initial = 0;%initial fiber length
signal_initial = 1*(10^-7);% inital signal power in W(-40dB) 
asef=0;%inital forward ase
aseb=0;%backward ase
x = x_initial ;
pumpi = pump_initial ;
signali = signal_initial ;
X1=[x];
PP1=[pumpi];
SP1=[signali];
ASE_F1=[asef] ;
ASE_B1=[aseb] ;
%Now lets do the RK4 methods, this is described in th report.
for a=0:n-1      
   x=x+h;
   
   k1=pump(x,pumpi,signali,asef,aseb);
   m1=signal(x,pumpi,signali,asef,aseb);
   f1=ASEF(x,pumpi,signali,asef,aseb); 
   b1=ASEB(x,pumpi,signali,asef,aseb);  
      
   k2=pump((x+0.5*h),(pumpi+0.5*k1*h),(signali+0.5*m1*h), (asef+0.5*f1*h),(aseb+0.5*b1*h)) ;
   m2=signal((x+0.5*h),(pumpi+0.5*k1*h),(signali+0.5*m1*h),(asef+0.5*f1*h),(aseb+0.5*b1*h));
   f2=ASEF((x+0.5*h),(pumpi+0.5*k1*h),(signali+0.5*m1*h), (asef+0.5*f1*h),(aseb+0.5*b1*h)) ;
   b2=ASEB((x+0.5*h),(pumpi+0.5*k1*h),(signali+0.5*m1*h), (asef+0.5*f1*h),(aseb+0.5*b1*h)) ;
  
   k3=pump((x+0.5*h),(pumpi+0.5*k2*h),(signali+0.5*m2*h), (asef+0.5*f2*h),(aseb+0.5*b2*h)) ;
   m3=signal((x+0.5*h),(pumpi+0.5*k2*h),(signali+0.5*m2*h), (asef+0.5*f2*h),(aseb+0.5*b2*h)) ;
   f3=ASEF((x+0.5*h),(pumpi+0.5*k2*h),(signali+0.5*m2*h), (asef+0.5*f2*h),(aseb+0.5*b2*h));
   b3=ASEB((x+0.5*h),(pumpi+0.5*k2*h),(signali+0.5*m2*h), (asef+0.5*f2*h),(aseb+0.5*b2*h));
   
   k4=pump((x+h),(pumpi+k3*h),(signali+m3*h), (asef+f3*h),(aseb+b3*h));
   m4=signal((x+h),(pumpi+k3*h),(signali+m3*h), (asef+f3*h),(aseb+b3*h));
   f4=ASEF((x+h),(pumpi+k3*h),(signali+m3*h), (asef+f3*h),(aseb+b3*h));
   b4=ASEB((x+h),(pumpi+k3*h),(signali+m3*h), (asef+f3*h),(aseb+b3*h));
   
   pumpi=pumpi+(1/6)*(k1+2*(k2+k3)+k4)*h;   
   signali=signali+(1/6)*(m1+2*(m2+m3)+m4)*h;   
   asef= asef +(1/6)*(f1+2*(f2+f3)+f4)*h;  
   aseb= aseb +(1/6)*(b1+2*(b2+b3)+b4)*h;   
    %here we return the value for each segment/loop into an array 
   X1=[X1,x];%position fiber length
   PP1=[PP1,pumpi];%pump power
   SP1=[SP1,signali];%signal power
   ASE_F1=[ASE_F1, asef];%forward ase
   ASE_B1=[ASE_B1, aseb];%backward ase
end
Gain=10*log10(SP1./signal_initial) ;
%Gain=ASE_F1;
Gainfinal=10*log10(signali/signal_initial) ; 
end
