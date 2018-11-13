%% Spectral analysis of first order example

% u'(x) = y u(x), u(1) = y u(xi)

N = 64;
[D,x] = chebdifmat(N,1,1);

I = eye(N+1);

% for v = 1:N+1
for v = N+1
    A = D;
    A(v,:) = I(1,:);
    z = zeros(1,2*N + 1);
    
    for k = -N:N
        z(k+N+1) = lambertw(k,x(v) - 1)./(x(v) - 1);
    end
    
    [vA, eigA] = eig(A);
    eigA = diag(eigA);
    
    B = PSIM(N,1,v,1,[]);
    eigB = 1./eig(B);
    
    figure(1)
    subplot(1,2,1)
    plot(real(eigA),imag(eigA),'bo',real(eigB),imag(eigB),'k^',real(z),imag(z),'r.')
    xlabel('Re(\lambda)')
    ylabel('Im(\lambda)')
    legend('Numerical','PSIM','Exact')
    title([ 'x_i = ' num2str(x(v)) ])
    axis([-50,500,-500,500])
    drawnow
%     M1{v} = frame2im(getframe);
    
    subplot(1,2,2)
    plot(real(eigA),imag(eigA),'bo',real(eigB),imag(eigB),'k^',real(z),imag(z),'r.')
    xlabel('Re(\lambda)')
    ylabel('Im(\lambda)')
    legend('Numerical','PSIM','Exact')
    title([ 'x_i = ' num2str(x(v)) ])
    axis([-1,2,-50,50])
    drawnow
%     M2{v} = frame2im(getframe);
    
%     [A1,map1] = rgb2ind(M1{v},256);
%     [A2,map2] = rgb2ind(M2{v},256);
%     if v ==1
%         imwrite(A1,map1,'SA1','gif','LoopCount',Inf)
%         imwrite(A2,map2,'SA2','gif','LoopCount',Inf)
%     else
%         imwrite(A1,map1,'SA1','gif','WriteMode','append')
%         imwrite(A2,map2,'SA2','gif','WriteMode','append')
%     end

    figure(2)
    ind = N+1;
    vAexact = exp(eigA(ind).*x);
    const = vA(1,ind)./vAexact(1);
    vAexact = const * vAexact;
    subplot(1,2,1)
    plot(x,real(vA(:,ind)),'ro',x,real(vAexact),'b--')
    xlabel('x')
    ylabel('Eigenvectors of A')
    subplot(1,2,2)
    plot(x,imag(vA(:,ind)),'ro',x,imag(vAexact),'b--')
    xlabel('x')
    ylabel('Eigenvectors of A')
    
    figure(3)
    plot(real(vA(:,N/2:end)),imag(vA(:,N/2:end)),'-.')
    xlabel('Re(v)')
    ylabel('Im(v)')
end