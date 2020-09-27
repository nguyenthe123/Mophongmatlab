%% Buoi_1_BTVN :Do hoa, tim khong diem, cuc tri ham 1 bien
% W : ndt1/duythe.276051@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.1 : Ve do thi
% 1.1a x(t) = 15cos(pi/2*t-pi/2) , t=[0 15]
% figure(1) Ve x,v,a bang fplot, dung legend
% figure(2) Ve khong gian pha, nhan xet nang luong dao dong
clc;clear all;close all;
fprintf('\nBai 1.1: Ve do thi x,v,a ham\n\nx = 15cos(pi/2*t-pi/2)\n');
x1 =@(t) 15*cos(pi/2*t - pi/2); d = [0 15];
figure(1); hold on;
fplot(x1,d);
xlabel('Time');
ylabel('Coordinates');
syms x % Find velocity
u = 15*cos(pi/2*x-pi/2);
v = diff(u,1)
a = diff(u,2) % Find acceleration
fplot(v,d);
fplot(a,d);
legend('x=15cos(pi/2*t-pi/2)','Velocity','Acceleration');
%figure(2) Plot phase space
d1 = 0:0.1:15;
x2 = eval(['@(x)',vectorize(u)]);
v1 = eval(['@(x)',vectorize(v)]);
p1 = -1*v1(d1); %(p=m*v, m=1)
x0 = x2(d1);
figure(2);hold on;
plot(x0,p1,'-b');
xlabel('Coordinate');
ylabel('Momentum');
legend('Phase space');
% => Energy is not change so it is harmonic oscillation
% 1.1b
clear all;
figure(3);hold on;
subplot(211);
t = 0:0.1:15;
x = 15.*exp(-0.45.*t).*cos(pi/2.*t - pi/2);
plot(t,x,'-r');
xlabel('Time');ylabel('Coordinates');
clear all;
subplot(212);
ezplot('15.*exp(-0.45.*t).*cos(pi/2.*t - pi/2)',[0 15]);
xlabel('Time');ylabel('Coordinates');hold on;
ylim([-5 10]);
title('Graph of 15.*exp(-0.45.*t).*cos(pi/2.*t - pi/2)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 1.2 RLC series U=220V, f=50Hz, C = (10^-4/(2pi))F, L=(1/pi)H
% R = x;
clear all;
fprintf('\nBai 1.2: Mach RLC\n');
%1.2a Xac dinh R de Pmax, find Pmax, so sanh Rm voi |Zc-Zl|
U=220;f=50;C=(10^-4/(2*pi));L=(1/pi);w=2*pi*f;
Zl=w*L;Zc=1/(w*C);
P=@(R) U^2.*R./(R.^2+(Zl-Zc)^2)
figure(4);hold on;
ezplot(P,[0 800])
P1=@(R) -U^2.*R./(R.^2+(Zl-Zc)^2);
[x0max P10max]=ginput;
for k=1:length(P10max)
    [xm(k,1) Pm(k,1)]=fminsearch(P1,x0max);
end
fprintf('Cuc dai Pmax = %f tai Rm = %.0f \n ',-Pm,xm);
fprintf('Nhan xet Rm = %.0f = Zl-Zc \n',xm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2b With P<Pmax, exist R1#R2 root, R1.R2=(Zl-Zc)^2. Bang pp Graph
fprintf('\n\nBang phuong phap do thi ta se giai phuong trinh P-P*=0\n')
q1=input('Nhap P* = q < P, P* = ','s');q=str2num(q1);
G=@(R) U^2.*R./(R.^2+(Zl-Zc)^2) - q
figure(5);hold on;
ezplot(G,[0 1500]);
[x0 fx0] = ginput;
fprintf(' Pt P-P*=0 se tro thanh q.R^2 - 220^2.R + q.100^2 = 0 \n');
fprintf(' De thay Delta = 220^4 - 4.q^2.100^2 = %.2f >0 \n',220^4-4*q^2*100^2);
fprintf(' Vay phuong trinh luon co 2 nghiem phan biet\n');
xn=fsolve(G,x0);xn1=100^2/xn;
fprintf('2 nghiem cua pt la : %.4f (Giai bang do thi) va %.4f\n',xn,xn1);
fprintf('Tu day ta thay R1*R2 = 100^2 = (Zl-Zc)^2 (dpcm) ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.2c Lam lai cau 1 co them r=50
fprintf('\n\n\nCuon cam co them r=50\n');r=50;
P=@(R) U^2.*R./((R+r).^2+(Zl-Zc)^2);figure(6);hold on;
ezplot(P,[0 800])
P1=@(R) -U^2.*R./((R+r).^2+(Zl-Zc)^2);
[x0max P10max]=ginput;
for k=1:length(P10max)
    [xm(k,1) Pm(k,1)]=fminsearch(P1,x0max);
end
fprintf('\nCuc dai Pmax = %f tai Rm = %.0f \n ',-Pm,xm);
fprintf('\nNhan xet Rm = %.0f \n',xm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 1.3
%1.3a Moi lien he dao dong dieu hoa va chuyen dong tron deu
clear all;
fprintf('\nBai 1.3: Lien he dao dong dieu hoa va chuyen dong tron deu: \n');
A = 15;
w = pi/2;
alpha = -pi/2;
t = linspace(0,15,1000);
x = A*cos(w*t+alpha);
d1 = linspace(0,15,1000);
syms e
u = 15*cos(pi/2*e-pi/2)
v = diff(u,1)
x2 = eval(['@(e)',vectorize(u)]);
v1 = eval(['@(e)',vectorize(v)]);
p1 = -1*v1(d1); %(p=m*v, m=1)
x0 = x2(d1);
figure(7);hold on;
h5=plot(x0(1),p1(1),'b.');hold on;
h1=plot(t(1),x(1),'b--',t(1),x(1),'r.');hold on;
h2=plot(0,x(1),'r')
h3=plot([0 t(1)],[x(1) x(1)],'r:');hold off
set(h5(1),'erase','non');
set(h1(1),'erase','non');
set(h1(2),'erase','xor','markersize',30);
set(h2,'erase','xor','markersize',30);
set(h3,'erase','xor');
axis([min(x0) max(x0) min(p1) max(p1)]);
d=0;
for k=2:length(t)
set(h5(1),'xdata',x0(1:k),'ydata',p1(1:k));
set(h1(1),'xdata',t(1:k),'ydata',x(1:k));
set(h1(2),'xdata',t(k),'ydata',x(k));
set(h2,'xdata',0,'ydata',x(k));
set(h3,'xdata',[0 t(k)], 'ydata',[x(k) x(k)]);
pause(0.05)
end
%1.3b
%%%%%%%%% ??????????????? %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%