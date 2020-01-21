load triangulos.mat 
fig=figure(1);
clf
set(fig,'DoubleBuffer','on');
PlotMesh(N_P,T_P,'b');
PlotMesh(N_S,T_S,'r');
perturb = eps * (rand(2,3)-0.5);
M=InterfaceMatrix(N_P,T_P,N_S + perturb,T_S);