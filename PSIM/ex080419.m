%% Examples for Seminar 08.04.19

%% u'' - u = cos(pi*x/2) with homo. Dirichlet BCs

NN = round(logspace(log10(2),log10(1024),50));
errm = zeros(size(NN));
erri = errm;
errL = errm;
for k = 1:length(NN)
    N  = NN(k);
    [Dm,xC] = chebdifmat(N,2,1);
    I  = eye(N+1);
    A  = Dm(:,:,end) - I;
    A([1,end],:) = I([1,end],:);
    F  = cos(pi*xC/2);
    F([1,end]) = 0;
    U  = A \ F;
    
    R  = IOMcc([1,0,-1],[1 0],[1 0],[1, N+1],N);
    UI = R * F;

    xL = linspace(-1,1,N+1);
    dx = 1./abs(xL(2) - xL(1));
    DL = ones(N+1,1);
    DL = dx^2 * [ DL -2*DL DL ];
    DL = spdiags(DL,[-1,0,1],N+1,N+1);
    AL = DL - I;
    AL([1,end],:) = I([1,end],:);
    FL = cos(pi*xL'/2);
    FL([1,end]) = 0;
    UL = AL \ FL;

    Em = -cos(pi*xC /2) / ((pi/2)^2 + 1);
    EL = -cos(pi*xL'/2) / ((pi/2)^2 + 1);

    errm(k) = max(abs(U - Em));
    erri(k) = max(abs(UI- Em));
    errL(k) = max(abs(UL- EL));

end

figure(1)
loglog(NN,errm,'ro--',NN,errL,'b*--','linewidth',2)
xlabel('N')
ylabel('Max error on collocation points')
legend('Chebyshev collocation','Finite differences')
set(gca,'fontsize',26,'linewidth',2)

figure(2)
loglog(NN,errm,'ro--',NN,errL,'b*--',NN,erri,'ks--','linewidth',2)
xlabel('N')
ylabel('Max error on collocation points')
legend('Chebyshev collocation','Finite differences','IOMcc')
set(gca,'fontsize',26,'linewidth',2)