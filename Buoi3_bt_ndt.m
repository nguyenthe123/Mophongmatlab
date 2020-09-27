%% Buoi_3_ BTVN : Cac bai toan giai tich co ban
% W : ndt1/ duythe.276051@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compose, finverse, limit, dif, int,quad,dblquad,taylor, symsum
%% Bai 3.1 Ham hop va ham nguoc
clc;clear all;close all
syms x y;
% 3.1.1a
disp('Bai 3.1.1a : ')
f(x) = x^2
g(x) = 2^x
fprintf('f(f(x)); g(g(x)); g(f(x)); f(g(x)) la: \n');
disp(compose(f,f)); disp(compose(g,g)); disp(compose(g,f)); disp(compose(f,g));
%compose(f,f) %f(f(x))
%compose(g,g) %g(g(x))
%compose(g,f) %g(f(x))
%compose(f,g) %f(g(x))
% 3.1.1b
fprintf('\nBai 3.1.1b : ')
f = sign(x)
g = 1/x
fprintf('f(f(x)); g(g(x)); g(f(x)); f(g(x)) la: \n');
fprintf('%s; %s; %s; %s\n',compose(f,f),compose(g,g),compose(g,f),compose(f,g));
% 3.1.2ab
disp('Bai 3.1.2ab : ');
f = (1-x)/(1+x)
g = (exp(x) - exp(-x))/2
fprintf('Ham nguoc f, g la: %s; %s \n',finverse(f));disp(finverse(g));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 3.2: Tinh cac gioi han
clear all;syms x
% 3.2 
disp('Bai 3.2: Tinh cac gioi han');
disp('Cac ham so abcdef va cac gioi han tai 2,1,+inf,+inf,a,+0 la:' );
a = ((x^2 - x - 2)^20)/((x^3 -12*x +16)^10)
limit(a,x,2)
b = (x^100 - 2*x + 1)/(x^50 - 2*x + 1)
limit(b,x,1)
c = (sqrt(x+sqrt(x+sqrt(x))))/(sqrt(x+1))
limit(c,x,inf,'left')
d = (sqrt(x) + power(x,1/3) + power(x,1/4))/(sqrt(2*x+1))
limit(d,x,inf,'left')
syms a 
e = (sin(x) - sin(a))/(x-a)
limit(e,x,a)
f = power(cos(sqrt(x)),x)
limit(f,x,0,'right')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 3.3: Tinh dao ham rieng va dao ham ham hop
clear all;syms x y
% 3.1.a Tinh dao ham rieng
a = ((x^2 - x -2)^20)/((x^3 - 12*x + 16)^10)
fprintf('Bai 3.1.a: Tinh dao ham rieng\n\t(A)1 : Tinh bang diff\n\t(A)2 : Tinh bang jacobian.\n');
A1 = [diff(a,x), diff(a,y)]
A2 = jacobian(a,[x y])
b = log((sqrt(x^2+y^2)-x)/(sqrt(x^2+y^2)+x))
B1 = [diff(b,x), diff(b,y)]
B2 = jacobian(b,[x y])
c = exp(x*y)*cos(x)*sin(y)
C1 = [diff(c,x), diff(c,y)]
C2 = jacobian(c,[x y])
% 3.1.b Tinh dao ham cua ham hop
fprintf('\nTinh dao ham ham hop:\n\tZ1: Tinh bang diff\n\tZ2: Tinh bang jacobian');
u = cos(x);
v = sqrt(x^2+y^2);
z = exp(u^2-2*v^2)
Z1 = [diff(z,x), diff(z,y)]
Z2 = jacobian(z,[x y])
fprintf('\nHam z khac');
u = x*y;
v = x/y;
z = log(u^2 + v^2)
Z1 = [diff(z,x), diff(z,y)]
Z2= jacobian(z,[x y])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 3.4: Tinh nguyen ham va tich phan
clear all;
% Nguyen ham
fprintf('\nTinh nguyen ham va tich phan\n');
syms x
fprintf('Nguyen ham: \n');
a = (x+1)/sqrt(x)
Antidiff_a = 0 + int(a,'x')
syms a
b = 1/(a^2 - x^2)
Antidiff_b = 0 + int(b,'x')
% Tich phan
fprintf('Tich phan:\n\tI1 theo quad\n\tI2 theo int \n');
disp('c = e^(sin^2(x)),0,pi/2')
c1 = @(x) exp(power(sin(x),2));
c = exp(power(sin(x),2));
I1 = quad(c1,0,pi/2)
I2 = int(c,'x',0,pi/2);
I2 = double(I2)
clear all; syms x;
disp('d = cos(x)/sqrt(1+x^4)');
d1 = @(x) cos(x)./sqrt(1+x.^4);
d = cos(x)/sqrt(1+x^4);
I1 = quad(d1,10,18)
I2 = int(d,'x',10,18);
I2=double(I2)
% Tich phan 2,3 lop
clear all;
disp('e = (y^3.e^y)/(x^2+y^2), 0<=x<=1, -4<=y<=2');
e = @(x,y) (y.^3.*exp(y))./(x.^2 + y.^2);
I = dblquad(e,0,1,-4,2)
disp('f = 1/sqrt(x^2 + y^2 + (z-2)^2),-1<=x<=1,-1.5<=y<=1.5,-1<=z<=1');
f = @(x,y,z) 1./(sqrt(x.^2 + y.^2 + (z-2).^2));
I = triplequad(f,-1,1,-1.5,1.5,-1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bai 3.5: Tinh khai trien Taylor quanh diem 0
clear all;
fprintf('\nBai 3.5: Tinh khai trien Taylor quanh diem 0 cac ham so f g h\n');
syms x
f = exp(x)
disp('Khai trien taylor bac 4 cua f quanh diem 0 la :');
taylor(f,x,0,'Order',4)
g = sinh(x)
disp('Khai trien taylor bac 4 cua g quanh diem 0 la :');
taylor(g,x,0,'Order',4)
h = cosh(x)
disp('Khai trien taylor bac 4 cua h quanh diem 0 la: ');
taylor(h,x,0,'Order',4)