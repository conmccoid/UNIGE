function B = PSIM(N,n,v,U1,U2)
%Computes the Birkhoff pseudospectral integration matrix
%   B = PSIM(N,n,v,U1,U2) is the nth order integration matrix for N+1
%   Chebyshev points with rows indexed by v removed and boundary conditions
%   perscribed by U1 (at x = 1) and U2 (at x = 2).
%   The length of v is n.
%
%   B = PSIM(N,n,[a,b,c,...],U1,U2) replaces the rows indexed
%   by a, b, c, ... with the boundary conditions specified in U1 and U2.
%   The rows are replaced in order, ie. if v = [a b c] and U1 has 2 rows,
%   rows a and b will be replaced with boundary conditions at x = 1, while
%   row c will be replaced with boundary conditions at x = -1.
%
%   ex. For Dirichlet boundary conditions on both boundaries, U1 = U2 = [1 0].
%   For Neumann boundary conditions, U = [0 1]. For mixed
%   boundary conditions, such as au(1)+bu'(1)=0 and cu(-1)+du'(-1)=0,
%   U1 = [a b] and U2 = [c d]. If only one side of the boundary is
%   specified, one of U1 or U2 can be left empty.

%---Chebyshev polynomials and integrals---%
% Polynomials (eq. 1.5)
T = cos( (0:N+2)'*(0:N)*pi/N );

% Integrals (eq. 2.11 and 2.12)
dTm=zeros(N+3,N+1,n+1);
dTm(:,:,1)=T;
for m=1:n
    dTm(1,:,m+1)=dTm(2,:,m);
    dTm(2,:,m+1)=dTm(3,:,m)/4;
    for k=3:N+3-m
        dTm(k,:,m+1)=((dTm(k+1,:,m)/k) - (dTm(k-1,:,m)/(k-2)))/2;
    end
end
dTm = permute(dTm,[2,1,3]);

%---Beta Coefficients---%
% Chebyshev weights (eq. 1.3 with modification)
c=ones(1,N+1);
c(1)=2; c(end)=2;
c = N*(c'*c);

% Coefficients (eq. 2.9)
b = T(1:N+1,:) - T(1:N+1,v)*( T(N-n+2:N+1,v)\T(N-n+2:N+1,:) );
b=2*b./c;

%---Boundary conditions---%
% Matrix system (eq. 2.16, 2.17)
Cm1 = zeros(n,n);
Cm1(1,:) = ones(1,n);
for k = 2:n
    Cm1(k,k:n) = ( ( (k:n) - 1).^2 - (k - 2).^2 ) / (2*k - 3);
    Cm1(k,k:n) = Cm1(k,k:n).*Cm1(k-1,k:n);
end
Cm2 = Cm1.*(-1).^( bsxfun(@plus,0:n-1,(0:n-1)') );

% Rhs (eq. 2.17)
dTm = permute(dTm,[3,2,1]);
P1 = dTm(n+1:-1:2,1:N+1,1);
P2 = dTm(n+1:-1:2,1:N+1,end);

% Free parameters (eq. 2.19, 2.23)
if isempty(U1)
    R = inv( U2*Cm2 );
    r = R*( U2*P2 );
elseif isempty(U2)
    R = inv( U1*Cm1 );
    r = R*( U1*P1 );
else
    R = inv( [ U1*Cm1 ; U2*Cm2 ] );
    r = R*[ U1*P1 ; U2*P2 ];
end

%---Construction of the PSIM---%
dTm = permute(dTm,[3,2,1]);
B = ( dTm(:,1:N+1,end) - dTm(:,1:n,1)*r ) * b ;     % (eq. 2.13)
B(:,v) = dTm(:,1:n,1) * R;                          % (eq. 2.21, 2.22)

end