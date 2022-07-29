function [N1,N2,N3] = PpIn(x,pump,signal,asef,aseb)
%This function is used to calculate the ion population of each energy
%levels
global pumpc;%The counter pumping power
global Nt;%doping concentration
global hc; 
global fp;%pump freq  
global fs;%signal frequency 
global Gamma_s;%signal overlap ratio 
global Gamma_p;%pump overlap ratio
global sigma_pa;
global sigma_pe;
global sigma_se;
global sigma_sa; 
global A_21; 
global AR;%fiber cross section area
A_21 =100;%emission rate
pumpc=0;%to calculate the bidirectional pump configuration,
%this value should be change.
Nt= 3*(10^24);
%The following equation are described in the report or you can find it in
%Corak et al.
R12=(sigma_pa*Gamma_p)/(hc*fp*AR)*(pump+pumpc);
R21=(sigma_pe*Gamma_p)/(hc*fs*AR)*(pump+pumpc);
W12=(sigma_sa*Gamma_s)/(hc*fp*AR)*(signal+asef+aseb);
W21=(sigma_se*Gamma_s)/(hc*fs*AR)*(signal+asef+aseb);
%calculate the population by equations in report
N1= Nt*(R21+W21+A_21)/(R21+R12+W21+W12+A_21);
N2= Nt*(R12+W12) /(R21+R12+W21+W12+A_21);
N3= Nt- N1- N2;


