% test the advancing front program
global fig;

[N1,T1]=NewMesh(1);
% N1 = N_P; T1 = T_P;
[N1,T1]=RefineMesh(N1,T1);
[N1,T1]=RefineMesh(N1,T1);
N1(:,5:end)=N1(:,5:end)+(0.1*rand(size(N1(:,5:end)))-0.05);
[N2,T2]=NewMesh(2);
% N2 = N_S; T2 = T_S;
[N2,T2]=RefineMesh(N2,T2);
[N2,T2]=RefineMesh(N2,T2);
N2(:,5:end)=N2(:,5:end)+(0.1*rand(size(N2(:,5:end)))-0.05);

fig=figure(1);
clf
set(fig,'DoubleBuffer','on');
PlotMesh(N1,T1,'b');
PlotMesh(N2,T2,'r');
M=InterfaceMatrix(N1,T1,N2,T2);