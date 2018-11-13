function Inverse = IOM(P,U1,U2,v)
%Computes the inverse operator matrix
%   R = IOM(P,U1,U2,v) is the inverse integration matrix for the linear
%   operator with homogeneous solutions defined in P, with rows indexed by
%   v removed and boundary conditions prescribed by U1 (at x = 1) and U2
%   (at x = 2). The length of v defines the order of the linear operator.
%
%   R = IOM(P,U1,U2,[a,b,c,...]) replaces the rows indexed by a, b, c, ...
%   with the boundary conditions specified in U1 and U2. The rows are
%   replaced in order, ie. if v = [a b c] and U1 has 2 rows, rows a and b
%   will be replaced with boundary conditions at x = 1, while row c will be
%   replaced with boundary conditions at x = -1.
%
%   Let the length of v be n. If P is of size (N+1) by n by n, then the
%   (j-1)th page of P gives the jth derivatives of the homogeneous
%   derivatives. If P is of size (N+1) by n by 1, the derivatives can be
%   calculated by multiplying by the Chebyshev differentiation matrices. In
%   either case, each column of P is a homogenous solution to the linear
%   operator of concern evaluated on the N+1 Chebyshev nodes.
%
%   If P is (N+1) by (N+1) then it is assumed that P is the collocation
%   matrix representing the linear operator. The homogeneous solutions are
%   then approximated by the null space of P.

%---Homog. Solns & their derivatives---%
% Size of inputs P and v
[N,type1,type2] = size(P);      % defines type of homog. solns
m = length(v);                  % defines order of linear operator

if type2==m                     % P provides derivatives
    [T,dT,cmat] = formCheb(N-1);
    Pm = P;
elseif type2==1                 % P does not provide derivatives
    [T,dT,cmat,D] = formCheb(N-1);
    
    if type1==N                 % P represents the linear operator
        I = eye(N);
        P = P \ I(:,v);
    end
    
    Pm = zeros(N,m,m);
    Pm(:,:,1) = P;
    for l = 2:m
        Pm(:,:,l) = D*Pm(:,:,l-1);
    end
end

% Particular homog. solns (eq. 3.20)
Coeffs = Wronskian(Pm,v);
if abs(det(Coeffs))<=eps
    disp('Error: choice of V does not support a fundamental set of solutions')
    Inverse = zeros(N,N);
    return
end

for i = 1:m
    Pm(:,:,i) = Pm(:,:,i)*Coeffs;
end

% Beta scalars (eq. 3.17)
W = Wronskian(Pm,1:N);

%---Interpolants---%
% G coefficient functions, homog. solns
B = zeros(N,N,m);
Inverse = zeros(N,N,m);
for k = 1:m
    B(:,:,k) = T' - T(v(k),:)'*T(:,N)'/T(v(k),N);     % (eq. 3.13)
    B(:,:,k) = cmat.*B(:,:,k);
    B(:,:,k) = dT*B(:,:,k);
    Inverse(:,:,k) = ( Pm(:,k,1)*W(k,:) ).*B(:,:,k);  % (eq. 3.18)
end
Inverse = sum(Inverse,3);

%---Boundary conditions---% (eq. 3.19)
Pm = permute( Pm,[3,2,1] );
B  = permute( B,[3,2,1] );

if isempty(U2)
    Lhs = U1*Pm(:,:,1);
    Rhs = Lhs*( W.*B(:,:,1) );
elseif isempty(U1)
    Lhs = U2*Pm(:,:,end);
    Rhs = Lhs*( W.*B(:,:,end) );
else
    BCplus  = U1*Pm(:,:,1); 
    BCminus = U2*Pm(:,:,end);
    Lhs = [ BCplus ; BCminus ];
    Rhs = [ BCplus*( W.*B(:,:,1) ) ; BCminus*( W.*B(:,:,end) ) ];
end

Pm = permute( Pm,[3,2,1] );
Inverse = Inverse - Pm(:,:,1)*( Lhs\Rhs );
Inverse(:,v) = Pm(:,:,1)*( Lhs\eye(m) );

end

function W = Wronskian(Fm,x) % (eq. 3.23)
% Wronskian - calculates the factors W_k(x)/W(x) for the functions Fm
% Inputs: Fm --- N+1 by m by m matrix containing the functions and their
%                derivatives
%         x  --- vector of indices at which to calculate the factors
M = length(x);
[~,m,~] = size(Fm);
Fm = permute(Fm,[3,2,1]);
W = zeros(m,M);

if m==1
   W = 1./permute(Fm(:,:,x),[1,3,2]); 
else
    for j = 1:M
        Num = det(Fm(:,:,x(j)));
        for k = 1:m
            ind = [ 1:k-1 k+1:m ];
            W(k,j) = (-1)^(k+m)*det( Fm(1:m-1,ind,x(j)) ) / Num;
        end
    end
end

end

function [T,dT,cmat,D] = formCheb(N)
% formCheb - forms the Chebyshev polynomials and integrals as well as the
%       differentiation matrix and the Chebyshev points

% Differentiation matrix (eq. 1.2)
if nargout==4
    D = chebdifmat(N,1,1);
end

% Polynomials (eq. 1.5)
T = cos( (0:N)'*(0:N+2)*pi/N );

% Integrals (eq. 2.11)
dT = zeros(N+1,N+1);
dT(:,1)=T(:,2);
dT(:,2)=T(:,3)/4;
for k=3:N+1
    dT(:,k)=(1/2)*((T(:,k+1)/k) - (T(:,k-1)/(k-2)));
end

T = T(1:N+1,1:N+1);

% Weights (eq. 1.3 with modification)
c = ones(1,N+1);
c(1) = 2*c(1); c(end) = 2*c(end);
cmat = c'*c; cmat = (2/N)./cmat;

end