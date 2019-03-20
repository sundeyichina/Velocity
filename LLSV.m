clc
clear
%comment
laserline = 26.5;%laser line length  /mm
num = 100; %number of point in the velocity curve
calibra = [];
%% velocity signal
xlrange = 'A22..B100021';
signal5 = csvread('Velocity.csv',21,0,xlrange);
time = signal5(:,1)';
velocity = signal5(:,2)';
Fvel = fopen('LLSV.txt','w+');
for i = 1:length(time)
    fprintf(Fvel,'%14.6E %14.6E\n',time(i),velocity(i));
end
fclose(Fvel);

%% Velocity data
% Read velocity signal

velocity_signal = textread('LLSV.txt');
Time_velocity = velocity_signal(:,1);
Voltage = velocity_signal(:,2);

figure
Time_velocity = Time_velocity*10^3;         % unit:ms
Voltage = Voltage*10^3;                   % unit:mV
plot(Time_velocity,Voltage,'k'); hold on;
title('$Projectile\:Velocity$','interpreter','latex','fontsize',20,'fontweight','bold');
xlabel('$Time\:(ms)$','interpreter','latex','fontsize',20,'fontweight','bold');
ylabel('$Velocity\:Signal(mV)$','interpreter','latex','fontsize',20,'fontweight','bold');
 xlim([-0.3,0.2]);
% ylim([0,400]);
% text(-0.1,300,'velocity = 107 m/s','Interpreter','latex','EdgeColor','red','fontsize',20,'fontweight','bold','BackgroundColor',[0 1 0]);
set(gca,'fontsize',20);
grid on
%% Get the useful segment of
tstart = -0.192;
tend = 0.02626;
index1 = find(Time_velocity == tstart);
index2 = find(Time_velocity == tend);
time = Time_velocity(index1 : index2) - Time_velocity(index1);
voltage = -(Voltage(index1 : index2) - Voltage(index1));
figure
plot(time,voltage,'k');
xlabel('$Time\:(ms)$','interpreter','latex','fontsize',20,'fontweight','bold');
ylabel('$Voltage\:(mV)$','interpreter','latex','fontsize',20,'fontweight','bold');
grid on
%% curve fit the function between displacement disp and voltage v
disp_voltage = fit(calibra(:, 1), calibra(:, 2),'poly2');
Disp = disp_voltage([calibra(1,1):0.02:calibra(length(calibra),1)]);
figure
plot([calibra(1,1):0.02:calibra(length(calibra),1)],Disp,'b-','linewidth', 1.5); hold on
plot(calibra(:, 1), calibra(:, 2),'ro','linewidth', 1.5)
title('$Displacement \: calibration$','interpreter','latex','fontsize',20,'fontweight','bold');
xlabel('$Voltage\:(mV)$','interpreter','latex','fontsize',20,'fontweight','bold');
ylabel('$Displacement\:(mm)$','interpreter','latex','fontsize',20,'fontweight','bold');
grid on

%% fit the disp - time curve 
displacement = disp_voltage(voltage);
t = time;
 x0 = [5 5 5];
 Disp = @(x, t) 1/2*x(1).*t.^2 + x(2).*t + x(3);
 x = lsqcurvefit(Disp, x0, t, displacement);
disp_t = fit(t, displacement,'poly2');
fitted_disp = disp_t(t);
velo = diff(fitted_disp)./diff(t);
 
 figure
 H = plot(t, displacement,'g','linewidth', 1.5);
 hold on
 [AX, H1, H2] = plotyy(t, Disp(x, t), t(1:length(t)-1), velo,'plot');
 set(AX(1),'XColor','k','YColor','b');
 set(AX(2),'XColor','k','YColor','r');
 HH1=get(AX(1),'Ylabel');
set(HH1,'String','$Displacement\:(mm)$','interpreter','latex','fontsize',20,'fontweight','bold');
set(HH1,'color','b');
HH2=get(AX(2),'Ylabel');
set(HH2,'String','$Velocity\:(m/s)$','interpreter','latex','fontsize',20,'fontweight','bold');
set(HH2,'color','r');
set(H1,'LineStyle','-','linewidth', 2);
set(H1,'color','b');
set(H2,'LineStyle','-','linewidth', 1.5);
set(H2,'color','r');
legend([H, H1,H2],{'Displacement';'Fitted displacement';'Velocity'},'fontsize', 13);
xlabel('$Time\:(ms)$','interpreter','latex','fontsize',20,'fontweight','bold');
grid on








