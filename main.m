clear
clc           % Clear the screen
format long
tic 
%
%-----PARAMETER INITIALIZATION------------------------------
%  NOTE :  Enter/Change your parameters here.
%
global Nt;   %Er ion density in ion/m^3 
%  50ppm =2*(10^24) per cubic meter
%
global lambda_p   ;  %m
global lambda_s    ; %m
global sigma_pa ; %m^2
global sigma_pe   ; %m^2
global sigma_pe2  ; %m^2
global sigma_se   ; %m^2
global sigma_sa  ; %m^2
global A_21    ; %s^-1
global A_32   ; %s^-1
global AR   ; %m^2
global Gamma_s   ; % Signal to core overlap
global Gamma_p  ; % Pump to core overlap
global h            ; % J/s ; Planck's constant
global c_light  ; % m/s ; Velocity of light
global v_p    ;  %Hz, Pump frequency
global v_s   ;  %Hz, Signal frequency
global alpha_s; % Background loss at signal wavelength (per m)=0.25dB/km
global alpha_p ; % Background loss at pump wavelength (per m)=0.25dB/km
FiberLength= 20 ; %m
n=140 ;          % No of Steps / No. of loops
hSTEP=(FiberLength-0)/n;           % h is step size, 
% (dy/dx)= y'= K(x,y,z);   where y(t=a)=yo and a<x<b
% (dz/dx)= z'= M(x,y,z);   where z(t=a)=yo and a<x<b
x_initial = 0  ;         % x=zl longitudinal distance
y_initial = 20*10^-3  ; % y= Pp , Pump power (W) ; Pp=20mW
yc=0 ;   % yc=Pp(counter):Counterpropagating Pump(W)
z_initial = 1*(10^-7) ;       % z= Ps , Signal power (W);
%Using Ps= -40dBm =(1*10^-4)mW
%
ASE_F1= 0 ;  % Forward ASE  
uf=ASE_F1 ;
ub=0;        % Backward ASE
Nt= 4*(10^24);   %Er ion density in ion/m^3 
%  50ppm =2*(10^24) per cubic meter
%
lambda_p =1480*(10^-9)   ;  %m
lambda_s =1550*(10^-9) ;   %m
sigma_pa =2.787671*(10^-21)*(10^-4)  ; %m^2
sigma_pe =0.8105639*(10^-21)*(10^-4)  ; %m^2
sigma_pe2 =0.8105639*(10^-21)*(10^-4) ; %m^2
sigma_se =4.118853*(10^-21)*(10^-4)  ; %m^2
sigma_sa =2.91056*(10^-21)*(10^-4)  ; %m^2
A_21 =100   ; %s^-1
A_32 =10^9  ; %s^-1
AR =1.633*(10^-11)   ; %m^2
Gamma_s =0.74  ; % Signal to core overlap
Gamma_p =0.77 ; % Pump to core overlap
h=6.626*(10^-34)             ; % J/s ; Planck's constant
c_light=3.0*(10^8)  ; % m/s ; Velocity of light
v_p =c_light/lambda_p   ;  %Hz, Pump frequency
v_s =c_light/lambda_s   ;  %Hz, Signal frequency
alpha_s=5.76e-5 ; % Background loss at signal wavelength (per m)=0.25dB/km
alpha_p=5.76e-5  ; % Background loss at pump wavelength (per m)=0.25dB/km
%

x = x_initial ;
y = y_initial ;
z = z_initial ;
%
LL1=[];
XX1=[x];
YY1=[y];
ZZ1=[z];
ASE_F1=[uf] ;
ASE_B1=[ub] ;
%
L1=0;
for w=0:n-1             % No. of loops
   LL1=[LL1,L1]; 
   L1=L1+1; 
   x=x+hSTEP;
   
   k1=pump(x,y,z,uf,ub);
   m1=signal(x,y,z,uf,ub);
   f1=ASEF(x,y,z,uf,ub); 
   b1=ASEB(x,y,z,uf,ub);  
      
   k2=pump((x+0.5*hSTEP),(y+0.5*k1*hSTEP),(z+0.5*m1*hSTEP), (uf+0.5*f1*hSTEP),(ub+0.5*b1*hSTEP)) ;
   m2=signal((x+0.5*hSTEP),(y+0.5*k1*hSTEP),(z+0.5*m1*hSTEP),(uf+0.5*f1*hSTEP),(ub+0.5*b1*hSTEP));
   f2=ASEF((x+0.5*hSTEP),(y+0.5*k1*hSTEP),(z+0.5*m1*hSTEP), (uf+0.5*f1*hSTEP),(ub+0.5*b1*hSTEP)) ;
   b2=ASEB((x+0.5*hSTEP),(y+0.5*k1*hSTEP),(z+0.5*m1*hSTEP), (uf+0.5*f1*hSTEP),(ub+0.5*b1*hSTEP)) ;
  
   k3=pump((x+0.5*hSTEP),(y+0.5*k2*hSTEP),(z+0.5*m2*hSTEP), (uf+0.5*f2*hSTEP),(ub+0.5*b2*hSTEP)) ;
   m3=signal((x+0.5*hSTEP),(y+0.5*k2*hSTEP),(z+0.5*m2*hSTEP), (uf+0.5*f2*hSTEP),(ub+0.5*b2*hSTEP)) ;
   f3=ASEF((x+0.5*hSTEP),(y+0.5*k2*hSTEP),(z+0.5*m2*hSTEP), (uf+0.5*f2*hSTEP),(ub+0.5*b2*hSTEP));
   b3=ASEB((x+0.5*hSTEP),(y+0.5*k2*hSTEP),(z+0.5*m2*hSTEP), (uf+0.5*f2*hSTEP),(ub+0.5*b2*hSTEP));
   
   k4=pump((x+hSTEP),(y+k3*hSTEP),(z+m3*hSTEP), (uf+f3*hSTEP),(ub+b3*hSTEP));
   m4=signal((x+hSTEP),(y+k3*hSTEP),(z+m3*hSTEP), (uf+f3*hSTEP),(ub+b3*hSTEP));
   f4=ASEF((x+hSTEP),(y+k3*hSTEP),(z+m3*hSTEP), (uf+f3*hSTEP),(ub+b3*hSTEP));
   b4=ASEB((x+hSTEP),(y+k3*hSTEP),(z+m3*hSTEP), (uf+f3*hSTEP),(ub+b3*hSTEP));
   
   y=y+(1/6)*(k1+2*(k2+k3)+k4)*hSTEP;   % Pump values
   z=z+(1/6)*(m1+2*(m2+m3)+m4)*hSTEP;   % Signal values
   uf= uf +(1/6)*(f1+2*(f2+f3)+f4)*hSTEP;   % Forward ASE values
   ub= ub +(1/6)*(b1+2*(b2+b3)+b4)*hSTEP;   

   XX1=[XX1,x];
   YY1=[YY1,y];
   ZZ1=[ZZ1,z];
   ASE_F1=[ASE_F1, uf] ;
   ASE_B1=[ASE_B1, ub] ;
   UF1=ASE_F1 ;
end
Gain=10*log10(ZZ1./z_initial) ;

figure(1)
plot(XX1,Gain,'linewidth',1.5);
xlabel('Fiber length (m)') % x-axis label
ylabel('pump1 (dB))') % y-axis label







