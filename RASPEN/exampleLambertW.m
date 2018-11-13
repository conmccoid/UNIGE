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
%   P1n=(P1n.^2+x*exp(-P1n))./(P1n+1);                           % Newton on P1  
  P1n = (P1n + 1)./(exp(P1n) + 1);
  P2=(x+P2)./(1+exp(P2));                                      % 2nd from book
  P2n=((P2n.^2+P2n*x+x).*exp(P2n)+x)./((x+P2n-1).*exp(P2n)-1); % Newton on 2nd
  P3=P3+x-P3*exp(P3);                                          % 3rd from book
  P3n=(P3n.^2.*exp(P3n)+x)./(exp(P3n).*(1+P3n)-1);
  P4=P4-x+P4*exp(P4);                                          % 4th from book
  P4n=(P4n.^2.*exp(P4n)+x)./(exp(P4n).*(1+P4n)+1);
  figure(1)
  subplot(2,4,1); image(xx,yy,abs(P1)); colormap hot;
  subplot(2,4,2); image(xx,yy,abs(P2)); colormap hot;
  subplot(2,4,3); image(xx,yy,abs(P3)); colormap hot;
  subplot(2,4,4); image(xx,yy,abs(P4)); colormap hot;
  subplot(2,4,5); image(xx,yy,abs(P1n)); colormap hot;
  subplot(2,4,6); image(xx,yy,abs(P2n)); colormap hot;
  subplot(2,4,7); image(xx,yy,abs(P3n)); colormap hot;
  subplot(2,4,8); image(xx,yy,abs(P4n)); colormap hot;
  
%   figure(2)
%   subplot(1,2,1)
%   contourf(xx,yy,log10(abs(P1n - lambertw(0,1))))
%   colorbar
%   subplot(1,2,2)
%   contourf(xx,yy,log10(abs(P1 - lambertw(0,1))))
%   colorbar
%   pause
end

% colormaps jet parula hot cool spring gray bone copper pink

  figure(2)
  subplot(1,2,1)
  contourf(xx,yy,log10(abs(P1n - lambertw(1,1))))
  colorbar
  subplot(1,2,2)
  contourf(xx,yy,log10(abs(P1 - lambertw(0,1))))
  colorbar
