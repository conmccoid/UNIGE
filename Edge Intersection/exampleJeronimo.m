load triangulos.mat 
fig=figure(1);
clf
set(fig,'DoubleBuffer','on');
PlotMesh(N_P,T_P,'b');
PlotMesh(N_S,T_S,'r');
M=InterfaceMatrix(N_P,T_P,N_S,T_S);