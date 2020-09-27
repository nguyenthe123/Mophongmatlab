%% Buoi_2_ BTVN : Giai hpt, tim cuc tri ham 2 bien
% W : ndt1/ duythe.276051@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 2.1 Tim nghiem gan dung pt phi tuyen thuc
% x^2sin(y) + cos(2x)y = 0.5
% 2^x.y - xy^2 - 1 = 0
% -2<=x<=2 , 0<=y<=2
clc;clear all;close all;
fprintf('Bai 2.1: Tim nghiem gan dung pt phi tuyen thuc: ');
h1=ezplot('x^2*sin(y) + cos(2*x)*y = 0.5', [-2,2,0,2]);hold on;
ezplot('2^x - x*y^2 - 1 = 0', [-2,2,0,2]);hold off;
f=@(x) [((x(1).^2)*sin(x(2)) + cos(2.*x(1)).*x(2) - 0.5) , (power(2,x(1)) - x(1).*(x(2)).^2 - 1)];
[x0,fx0] = ginput;
for k=1:length(x0)
    [xn(k,:),fxn(k,:)] = fsolve(f,[x0(k),fx0(k)])
end
xyn = [xn(:),fxn(:)]
fprintf('/n Ta thay 2 so dau tien la 1 cap nghiem (x,y)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 2.2 Tim nghiem gan dung pt phi tuyen phuc
% z^2 + cos(ln(|z| +1) + 3) - 2^z - i = 0
% -2 <= Re,Im <= 2 bang pp Graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
fprintf('\nBai 2.2: Tim nghiem gan dung pt phi tuyen phuc');
fz = @(z) (z.^2 + cos((log(abs(z)) + 1) + 3) - 2.^z - i);
f1 = @(x,y) real(fz(x+1i*y)); f2=@(x,y) imag(fz(x+1i*y));
h1 = ezplot(f1,[-2,2]);hold on;
ezplot(f2,[-2,2]);hold off;
[x0,y0] = ginput ; z0=x0+i*y0;
for k=1:length(z0)
    [zn(k),zyn(k)] = fsolve(fz,z0(k));
end
disp('Nghiem z = '); disp(zn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 2.3
% z = f(x,y) = (2x^2+3y^2)e^(-(x^2+y^2))
% 2.3a -5<=x,y<=5
% 2.3b Tinh truong vector gradient cua z va ve chung tren ham quiver
% 2.3c Tinh cac cuc tri cua ham nay
clear all;
fprintf('\nBai 2.3: Ve gradient va tim cuc tri: ');
[x,y]=meshgrid(-5:.15:5);
z=(2*x^2+3*y^2)*exp(-(x^2+y^2));
[px,py]=gradient(z,0.15,0.15);
figure(1);surf(x,y,z);
figure(2);contour(x,y,z,30);hold on;grid on;
quiver(x,y,px,py),hold off, axis image
% Min fx
f=@(x) ((2.*x(1).^2+3.*x(2).^2).*exp(-(x(1).^2+x(2).^2)));
[x10,y10]=ginput;
disp('Cuc tieu');
for k=1:length(x10)
    [xct(k,:),yct(k,:)] = fminsearch(f,[x10(k),y10(k)]);
end
xct1_2_yct=[xct,yct]
% Max fx
fcd = @(x) -((2.*x(1).^2+3.*x(2).^2).*exp(-(x(1).^2+x(2).^2)));
[x20,y20]=ginput;
disp('Cuc dai');
for k=1:length(x20)
    [xcd(k,:),ycd(k,:)] = fminsearch(fcd,[x20(k),y20(k)]);
end
xcd_1_2_ycd=[xcd,-ycd]