function  [DM,xc] = chebmatb(N,M,nst);
% CHEBMATB  Compute Chebyshev spectral differentiation matrices
%
%   Use text book formulas 
%      
%     nst=0:  no correction
%     nst=1:  Negative sum trick NST  (default)
%
% DM = chebmatb(N,M,nst)    [Each DM(:,:,ell) is an N+1 by N+1 matrix], or
% [DM,xc] = chebmatnst(N,M,nst)
%
% N: Degree of approximation  (N+1  points)
% M: Maximum order of differentiation (default=1)
% nst: optional argument, controlling "negative sum trick"
%
% D = first order differentiation,  
%
if (nargin<3) nst=1; if (nargin<2) M=1; end; end;
k=[0:N]';  x=cos(pi*k/N);                 % Chebyshev points
L=logical(eye(N+1));
C = toeplitz((-1).^k);                    % C is the matrix with 
C(1,:) = C(1,:)*2; C(N+1,:) = C(N+1,:)*2; % entries (-1)^(j+k)*c(k)/c(j)
C(:,1) = C(:,1)/2; C(:,N+1) = C(:,N+1)/2;
X=repmat(x,1,N+1);                        % Off-diagonal elements
E=X-X'; E(L)=ones(N+1,1); 
Z=1./E; Z(L)=zeros(N+1,1);                % Z has elements 1/(x(k)-x(j)) 
D=eye(N+1);
for ell=1:M,
   D = ell*Z.*(C.*repmat(diag(D),1,N+1) - D);
   for k=1:N+1,
      ind=[1:k-1,k+1:N+1]; y=D(k,ind);
      [yabs,ii]=sort(abs(D(k,ind)));
      D(k,k) = -sum(y(ii));               % NST with sorting
   end;
   Dnew=D;
%    Dnew(floor(N/2)+2:end,:)=-Dnew(ceil(N/2):-1:1,end:-1:1); % Flipping trick
   DM(:,:,ell)=Dnew;
end;
if (nst==0),
   dig = x(2:N)./(ones(N-1,1)-x(2:N).^2); corner=(2*N^2 + 1)/6;
   dig = [corner -dig'/2 -corner]';          % Diagonal entries
   D=DM(:,:,1); D(L)=dig; DM(:,:,1)=D; E=D;
   for ell=2:M,
      dig=sum(E'.*D);  
      E=DM(:,:,ell); E(L)=dig; DM(:,:,ell)=E;
   end;
end;   %if
if nargout > 1, xc=x; 
% xc(floor(N/2)+2:end)=-xc(ceil(N/2):-1:1); %Flipping trick
end;
%END
