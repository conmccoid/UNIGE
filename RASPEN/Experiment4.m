%% Comparing FPI

g1 = @(z) exp(-z);
g2 = @(z) -log(z);

xx = linspace(-5,100,10000);
yy = linspace(-pi,pi,1000);
[X,Y] = meshgrid(xx,yy);
z = X + 1i*Y;

%% Basins of attraction (analytic)

dx = xx(xx>0);
dy = yy(abs(yy)<pi/2);
[DX,DY] = meshgrid(dx,dy);
D1 = DX + 1i*DY;
D1 = D1(:);

% D1 = [ (0:0.001:10)-1i*pi/2,...
%     10+1i*(-pi/2:0.001:pi/2),...
%     (0:0.001:10)+1i*pi/2,...
%     1i*(-pi/2:0.001:pi/2) ];

g11 = g1(D1);
g12 = g1(g11);
g13 = g1(g12);
g21 = g2(D1);
g22 = g2(g21);
g23 = g2(g22);

figure(1)
plot(real(g23),imag(g23),'.',...
    real(g22),imag(g22),'.',...
    real(g21),imag(g21),'.',...
    real(D1),imag(D1),'.',...
    real(g11),imag(g11),'.',...
    real(g12),imag(g12),'.',...
    real(g13),imag(g13),'.')
axis([-5,10,-pi,pi])
xlabel('Re(z)')
ylabel('Im(z)')
% legend('g_2(D)','D','g_1(D)')

%% Basins of attraction (numeric)
xx = linspace(-10,10,1000);
yy = linspace(-pi,pi,1000);
[X,Y] = meshgrid(xx,yy);
z = X + 1i*Y;

for i = 1:20
    z = g1(z);
end

figure(2)
plot(real(g23),imag(g23),'.')
hold on
contour(xx,yy,log10(abs(z-lambertw(0,1))),-7:0)
hold off

%% Precon. Newton

f1 = @(z) (1 + z)./(1 + exp(z));
f2 = @(z) z.*(1 - log(z))./(z + 1);
f1inv = @(z) z - 1 - lambertw(0,-z.*exp(z-1));
% f2inv = @(z) -z./lambertw(0,-z.*exp(z-1));
f2inv = @(z) exp(-f1inv(z));

xx = linspace(-10,10,1001);
yy = linspace(-10,10,1001);
[X,Y] = meshgrid(xx,yy);
z = X + 1i*Y;
% z = z(:);

contourf(xx,yy,log10(abs(f2(f2inv(z)) - z)))

%% Basins of attraction (analytic)

f1prime = @(z) (1 - z.*exp(z))./(1 + exp(z)).^2;
f2prime = @(z) (z + log(z))./(1 + z).^2;

xx = linspace(-10,10,1001);
yy = linspace(-10,10,1001);
[X,Y] = meshgrid(xx,yy);
z = X + 1i*Y;

figure(1)
contour(xx,yy,abs(f1prime(z)),[1,1],'r')
hold on
contour(xx,yy,abs(f2prime(z)),[1,1],'b')
hold off

figure(2)
alpha = 1;
zz = z(f1prime(z)<=1 & abs(z)<alpha);
plot(real(f1(zz)),imag(f1(zz)),'.')
hold on
contour(xx,yy,abs(f1prime(z)),[1,1])
contour(xx,yy,abs(z),[alpha,alpha])
hold off
axis([-10,10,-10,10])

figure(3)
zz = z(abs(z)<1);
D1 = f1inv(zz);
D2 = f1inv(D1);
D3 = f1inv(D2);
d1 = f1(zz);
d2 = f1(d1);
d3 = f1(d2);
plot(real(D3),imag(D3),'.',...
    real(D2),imag(D2),'.',...
    real(D1),imag(D1),'.',...
    real(zz),imag(zz),'.',...
    real(d1),imag(d1),'.',...
    real(d2),imag(d2),'.',...
    real(d3),imag(d3),'.')
axis([-10,10,-10,10])
xlabel('Re(z)')
ylabel('Im(z)')

%% Basins of attraction (numerical)

xx = linspace(-10,10,1001);
yy = linspace(-10,10,1001);
[X,Y] = meshgrid(xx,yy);
z = X + 1i*Y;

for i = 1:20
    z = f1(z);
end

plot(real(D3),imag(D3),'.')
hold on
contour(xx,yy,log10(abs(z - lambertw(0,1))))
hold off

%% Highly resolved basins of attraction

rr = linspace(0,1,1e4);
th = linspace(0,2*pi,1e4);
[R,T] = meshgrid(rr,th);
z = R.*exp(1i*T);
z = z(:);

D1 = f1inv(z);
D2 = f1inv(D1);
D3 = f1inv(D2);

plot(real(D3),imag(D3),'.',...
    real(D2),imag(D2),'.',...
    real(D1),imag(D1),'.',...
    real(z),imag(z),'.')
axis([-10,10,-10,10])
xlabel('Re(z)')
ylabel('Im(z)')