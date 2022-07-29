%The main code for ploting the performance characteristics
%first we plot the gain versus fiber length
%we vary the fiber length with pump power constant
%The pump power change from 10mW to 200 mW
%This can also be used to plot the ASE power, just change the output in the
%end of func.m from gain to ASEpower
FiberLength1=50;
FiberLength2=40;
FiberLength100=60;
pump_initial1=10*10^-3;
pump_initial2=15*10^-3;
pump_initial3=20*10^-3;
pump_initial4=40*10^-3;
pump_initial5=60*10^-3;
pump_initial6=80*10^-3;
pump_initial7=100*10^-3;
pump_initial8=200*10^-3;

[Gain1,Gainfinal,X1,PP1,SP1] = func(FiberLength100,pump_initial1);
[Gain2,Gainfinal,X12,PP1,SP1] = func(FiberLength100,pump_initial2);
[Gain3,Gainfinal,X13,PP1,SP1] = func(FiberLength100,pump_initial3);
[Gain4,Gainfinal,X14,PP1,SP1] = func(FiberLength100,pump_initial4);
[Gain5,Gainfinal,X15,PP1,SP1] = func(FiberLength100,pump_initial5);
[Gain6,Gainfinal,X16,PP1,SP1] = func(FiberLength100,pump_initial6);
[Gain7,Gainfinal,X17,PP1,SP1] = func(FiberLength100,pump_initial7);
[Gain8,Gainfinal,X18,PP1,SP1] = func(FiberLength100,pump_initial8);
%plot the figures. 
figure(1)
plot(X1,Gain1,'linewidth',1.5);
hold on
plot(X12,Gain2,'linewidth',1.5)
hold on
plot(X13,Gain3,'linewidth',1.5);
hold on
plot(X14,Gain4,'linewidth',1.5);
hold on
plot(X15,Gain5,'linewidth',1.5);
hold on
plot(X16,Gain6,'linewidth',1.5);
hold on
plot(X17,Gain7,'linewidth',1.5);
hold on
plot(X18,Gain8,'linewidth',1.5);
xlabel('Fiber length (m)') 
%The y-axis label and title can be changed to plot the ASEpower
ylabel('Gain (dB))') 
%ylabel('ASE Power (W))') 
legend('10mW','15mW','20mW','40mW','60mW','80mW','100mW','200mW');
title('EDFA simulated Gain vs Fiber length of 1480nm');
%title('EDFA simulated ASE Power vs Fiber length of 1480nm');


%here we plot the gain vs pump power
%The pump power is set from 0 to 100mW, every 2mW it will run a loop and
%return a value for the final gain(the gain at the end of the fiber)
pump_initials=0:2e-3:100e-3;%vary inital pump power

FiberLengtha=10;
FiberLengthb=20;
FiberLengthc=30;
FiberLengthd=40;
FiberLengthe=50;
for i=1:1:length(pump_initials)
    [Gain,Gainfinala,X1,PP1,SP1] = func(FiberLengtha,pump_initials(i));
    [Gain,Gainfinalb,X1,PP1,SP1] = func(FiberLengthb,pump_initials(i));
    [Gain,Gainfinalc,X1,PP1,SP1] = func(FiberLengthc,pump_initials(i));
    [Gain,Gainfinald,X1,PP1,SP1] = func(FiberLengthd,pump_initials(i));
    [Gain,Gainfinale,X1,PP1,SP1] = func(FiberLengthe,pump_initials(i));
        %store the final gin value for each loop in gainfinalabcde
    gainfinala(i)=Gainfinala;
    gainfinalb(i)=Gainfinalb;
    gainfinalc(i)=Gainfinalc;
    gainfinald(i)=Gainfinald;
    gainfinale(i)=Gainfinale;
end
%plot the figures for gain vs pump
figure(2)
plot(pump_initials,gainfinala,'linewidth',1.5);
hold on
plot(pump_initials,gainfinalb,'linewidth',1.5);
hold on
plot(pump_initials,gainfinalc,'linewidth',1.5);
hold on
plot(pump_initials,gainfinald,'linewidth',1.5);
hold on
plot(pump_initials,gainfinale,'linewidth',1.5);
xlabel('Pump power (W)') 
ylabel('Gain (dB))') 
legend('10m','20m','30m','40m','50m');
title('EDFA simulated Gain vs Pump power of 1480nm');



