function [Gain,Gainfinal,X1,PP1,SP1] = func1(FiberLength,pump_initial)

format long

global Nt;   %Er ion density in ion/m^3 
%  50ppm =2*(10^24) per cubic meter
%
global lambda_p    %m
global lambda_s    ; %m
global sigma_pa ; %m^2
global sigma_pe   ; %m^2
global sigma_pe2  ; %m^2
global sigma_se   ; %m^2
global sigma_sa  ; %m^2
global A_21    ; %s^-1
%global A_32   ; %s^-1
global AR   ; %m^2
global Gamma_s   ; % Signal to core overlap
global Gamma_p  ; % Pump to core overlap
global hc            ; % J/s ; Planck's constant
global c_light  ; % m/s ; Velocity of light
global fp    ;  %Hz, Pump frequency
global fs   ;  %Hz, Signal frequency
global alpha_s; % Background loss at signal wavelength (per m)=0.25dB/km
global alpha_p ; % Background loss at pump wavelength (per m)=0.25dB/km
%FiberLength= 20 ; %m
n=200 ;          % No of Steps / No. of loops
h=(FiberLength-0)/n;           % h is step size, 
% (dy/dx)= y'= K(x,y,z);   where y(t=a)=yo and a<x<b
% (dz/dx)= z'= M(x,y,z);   where z(t=a)=yo and a<x<b
x_initial = 0  ;         % x=zl longitudinal distance
%y_initial = 20*10^-3  ; % y= Pp , Pump power (W) ; Pp=20mW
yc=0 ;   % yc=Pp(counter):Counterpropagating Pump(W)
signal_initial = 1*(10^-7) ;       % z= Ps , Signal power (W);
%Using Ps= -40dBm =(1*10^-4)mW
%
ASE_F1= 0 ;  % Forward ASE  
asef=ASE_F1 ;
aseb=0;        % Backward ASE
Nt= 3*(10^24);   %Er ion density in ion/m^3 
%  50ppm =2*(10^24) per cubic meter
%
lambda_p =980*(10^-9)  ;  %m
lambda_s =1550*(10^-9) ;   %m
sigma_pa =0.75*(10^-25)  ; %m^2 2.787671*(10^-21)*(10^-4)
sigma_pe =0.19*(10^-25)  ; %m^2 0.8105639*(10^-21)*(10^-4)
sigma_pe2 =0.8105639*(10^-21)*(10^-4) ; %m^2
sigma_se =3.8*(10^-25)  ; %m^2 4.118853*(10^-21)*(10^-4)
sigma_sa =2.4*(10^-25)  ; %m^2 2.91056*(10^-21)*(10^-4)
A_21 =100   ; %s^-1
%A_32 =10^9  ; %s^-1
AR =1.25*(10^-11)   ; %m^2 1.633*(10^-11) 
Gamma_s =0.4  ; % Signal to core overlap 0.74
Gamma_p =0.4 ; % Pump to core overlap 0.4
hc=6.626*(10^-34)             ; % J/s ; Planck's constant
c_light=3.0*(10^8)  ; % m/s ; Velocity of light
fp =c_light/lambda_p   ;  %Hz, Pump frequency
fs =c_light/lambda_s   ;  %Hz, Signal frequency
alpha_s=5.756e-5 ; % Background loss at signal wavelength (per m)=0.25dB/km
alpha_p=5.756e-5  ; % Background loss at pump wavelength (per m)=0.25dB/km
%

x = x_initial ;
pumpi = pump_initial ;
signali = signal_initial ;
%
LL1=[];
X1=[x];
PP1=[pumpi];
SP1=[signali];
ASE_F1=[asef] ;
ASE_B1=[aseb] ;
%
L1=0;
for w=0:n-1             % No. of loops
   LL1=[LL1,L1]; 
   L1=L1+1; 
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
   
   pumpi=pumpi+(1/6)*(k1+2*(k2+k3)+k4)*h;   % Pump values
   signali=signali+(1/6)*(m1+2*(m2+m3)+m4)*h;   % Signal values
   asef= asef +(1/6)*(f1+2*(f2+f3)+f4)*h;   % Forward ASE values
   aseb= aseb +(1/6)*(b1+2*(b2+b3)+b4)*h;   

   X1=[X1,x];
   PP1=[PP1,pumpi];
   SP1=[SP1,signali];
   ASE_F1=[ASE_F1, asef] ;
   ASE_B1=[ASE_B1, aseb] ;
   UF1=ASE_F1 ;
end
Gain=10*log10(SP1./signal_initial) ;


Gainfinal=10*log10(signali/signal_initial) ; 
end
