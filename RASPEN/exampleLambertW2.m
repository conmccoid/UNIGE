x=1;                                                           % try 1, 10, 50
N=20;                                                          % no. of iter.
n=1000;                                                        % resolution
m=1000;                                                        % fig. range
a=-10;b=10;c=-10;d=10;
%a=9;b=9.4;c=2.8;d=3.2;
xx=linspace(a,b,n); yy=linspace(c,d,m);
[X,Y]=meshgrid(xx,yy);
P1=X+1i*Y;P2=P1;P3=P1;P4=P1;P1n=P1;P2n=P1;P3n=P1;P4n=P1;
for i=1:N
  P1=x*exp(-P1);                                               % 1st from book
  %P1n=(P1n.^2+x*exp(-P1n))./(P1n+1);                          % ???
  P1n=(exp(-P1n)*x.*(P1n+1))./(x*exp(-P1n)+1);                   % Newton on P1  
  P2=(x+P2)./(1+exp(P2));                                      % 2nd from book
  %P2n=((P2n.^2+P2n*x+x).*exp(P2n)+x)./((x+P2n-1).*exp(P2n)-1); % ???
  P2n=(P2n.^2+P2n*x+exp(-P2n)*x+x)./(1+x+P2n+exp(P2n));         % Newton on 2nd
  P3=P3+x-P3.*exp(P3);                                         % 3rd from book
  %P3n=(P3n.^2.*exp(P3n)+x)./(exp(P3n).*(1+P3n)-1);            % ??? 
  P3n=(P3n.^2+exp(-P3n)*x)./(1+P3n);
  P4=P4-x+P4.*exp(P4);                                        % 4th from book
  P4n=(P4n.^2+exp(-P4n)*x)./(1+P4n);
%   figure(2)
%   subplot(2,4,1); image(xx,yy,abs(P1)); colormap hot;
%   subplot(2,4,2); image(xx,yy,abs(P2)); colormap hot;
%   subplot(2,4,3); image(xx,yy,abs(P3)); colormap hot;
%   subplot(2,4,4); image(xx,yy,abs(P4)); colormap hot;
%   subplot(2,4,5); image(xx,yy,abs(P1n)); colormap hot;
%   subplot(2,4,6); image(xx,yy,abs(P2n)); colormap hot;
%   subplot(2,4,7); image(xx,yy,abs(P3n)); colormap hot;
%   subplot(2,4,8); image(xx,yy,abs(P4n)); colormap hot;
%   pause
end

figure(1)
  subplot(2,4,1); contourf(xx,yy,abs(P1))
  subplot(2,4,2); contourf(xx,yy,abs(P1n))
  subplot(2,4,3); contourf(xx,yy,abs(P2))
  subplot(2,4,4); contourf(xx,yy,abs(P2n))
  subplot(2,4,5); contourf(xx,yy,abs(P3))
  subplot(2,4,6); contourf(xx,yy,abs(P3n))
  subplot(2,4,7); contourf(xx,yy,abs(P4))
  subplot(2,4,8); contourf(xx,yy,abs(P4n))

% colormaps jet parula hot cool spring gray bone copper pink

%% Testing changes in theoretical basin of attraction

n=1000;
m=1000;
a=-10;b=10;c=-10;d=10;
xx=linspace(a,b,n); yy=linspace(c,d,m);
[X,Y]=meshgrid(xx,yy);
Z = X + 1i*Y;

% g1 = @(z) exp(-z);
% g2 = @(z) exp(-z).*( exp(-z) - z )./( exp(-z) + 1 ).^2;

% g1 = @(z) (1 - z.*exp(z))./(1 + exp(z)).^2;
% g2 = @(z) -(exp(z).*(z-1) - z - 3).*(1 - exp(z).*z) ./ (z + exp(z) + 2).^2;

g1 = @(z) z + 1 - z.*exp(z);
g2 = @(z) z - (z.*exp(z) - 1)./( (z+1).*exp(z) );

figure(3)
contour(xx,yy,abs(g1(Z)),[1,1],'r')
hold on
contour(xx,yy,abs(g2(Z)),[1,1],'b')
hold off
xlabel('Re(z)')
ylabel('Im(z)')
title('Theoretical basins of attraction')

%%
figure(1)
subplot(2,4,5)
hold on
contour(xx,yy,abs(g1(Z)),[1,1],'r')
hold off
subplot(2,4,6)
hold on
contour(xx,yy,abs(g2(Z)),[1,1],'b')
hold off