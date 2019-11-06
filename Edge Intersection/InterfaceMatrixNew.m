function M=InterfaceMatrixNew(Na,Ta,Nb,Tb)
% INTERFACEMATRIX projection matrix for nonmatching triangular grids 
%   M=InterfaceMatrix(Na,Ta,Nb,Tb); takes two triangular meshes Ta
%   and Tb with associated nodal coordinates in Na and Nb and
%   computes the interface projection matrix M

% bl and bil need to be changed so that the first pair of triangles has an
% intersection - this seems like a real crutch
% replacing the new meshes in example with those from Jeronimo requires bl
% to be 15 (or so?)

bl =1;                           % bl: list of triangles of Tb to treat
bil=1;                           % bil: list of triangles Ta to start with
bd=zeros(size(Tb,1)+1,1);        % bd: flag for triangles in Tb treated 
bd(end)=1;                       % guard, to treat boundaries
bd(1)=1;                         % mark first triangle in b list.
M=sparse(size(Nb,2),size(Na,2));
while isempty(bl)~=1
    bc=bl(1); bl=bl(2:end);      % bc: current triangle of Tb 
    al=bil(1); bil=bil(2:end);   % triangle of Ta to start with
    ad=zeros(size(Ta,1)+1,1);    % same as for bd
    ad(end)=1;
    ad(al)=1;
    n=[0 0 0];                   % triangles intersecting with neighbors
    while ~isempty(al)
        ac=al(1); al=al(2:end);  % take next candidate
        [P,nc,Mc]=Intersect(Nb(:,Tb(bc,1:3)),Na(:,Ta(ac,1:3)));
        if ~isempty(P)           % intersection found
            M(Tb(bc,1:3),Ta(ac,1:3))=M(Tb(bc,1:3),Ta(ac,1:3))+Mc;
            t=Ta(ac,3+find(ad(Ta(ac,4:6))==0)); 
            al=[al t];           % add neighbors 
            ad(t)=1;
            n(nc>0)=ac;          % ac is starting candidate for neighbor  
        end
    end
    tmp=find(bd(Tb(bc,4:6))==0); % find non-treated neighbors
    idx=find(n(tmp)>0);          % take those which intersect
    t=Tb(bc,3+tmp(idx));
    bl=[bl t];                   % and add them
    bil=[bil n(tmp(idx))];       % with starting candidates Ta
    bd(t)=1;
end
end

function [P,n,M]=Intersect(X,Y)
% INTERSECT intersection of two triangles and mortar contribution
%   [P,n,M]=Intersect(X,Y); computes for the two given triangles X and
%   Y (point coordinates are stored column-wise, in counter clock
%   order) the points P where they intersect, in n the indices of
%   which neighbors of X are also intersecting with Y, and the local
%   mortar matrix M of contributions of the P1 elements X on the P1
%   element Y. The numerical challenges are handled by including
%   points on the boundary and removing duplicates at the end.

[P,Q,R,n] = Geometry(X,Y);
% plot(P(1,:),P(2,:),'k^',Q(1,:),Q(2,:),'rx',R(1,:),R(2,:),'bo')
% triX = [ X , X(:,1)];
% triY = [ Y , Y(:,1)];
% hold on
% plot(triX(1,:),triX(2,:),'-',triY(1,:),triY(2,:),'-')
% hold off
% pause
if size(P,2)>1                         % if two or more interior points
    n=[1,1,1];                         % the triangle is candidate for all
end                                    % neighbors
P=[P Q R];

P=SortAndRemoveDoubles(P);             % sort counter clock wise
M=zeros(3,3);
if size(P,2)>0
    for j=2:size(P,2)-1                % compute interface matrix
        M=M+MortarInt(P(:,[1 j j+1]),X,Y);
    end
    patch(P(1,:),P(2,:),'m')           % draw intersection for illustration
%  H=line([P(1,:) P(1,1)],[P(2,:),P(2,1)]);
%  set(H,'LineWidth',3,'Color','m');
    pause(0)
end
end

function [P,Q,R,n] = Geometry(X,Y)
% GEOMETRY figures out the intersection of the triangles X and Y. It does
% so by first putting all variables in terms of two of the lines of the
% triangle Y, decides if the points of X lie within or outside Y, then
% calculates the intersections of the lines of X with those of Y and again
% deciding if they are inside or outside Y.

ind=[0,0,0];
v0=Y(:,2)-Y(:,1);  % vector between y1 and y2
v1=Y(:,3)-Y(:,1);  % vector between y1 and y3
d00=v0'*v0;        % norm of v0
d01=v0'*v1;
d11=v1'*v1;        % norm of v1
id =d00*d11 - d01*d01; % this value is used to normalize quantities

%---Interior points---%
U = zeros(2,3);
for i=1:3 
    r=X(:,i)-Y(:,1); % vector between ith vertex of triangle X and y1
    d02=v0'*r;       % projection of v2 in v0 direction
    d12=v1'*r;       % projection of v2 in v1 direction
    U(1,i)=(d11*d02-d01*d12)/id;  % amount of v2 in v0 direction
    U(2,i)=(d00*d12-d01*d02)/id;  % amount of v2 in v1 direction
    % v2 = u*v0 + v*v1
    % if u and v are positive and their sum less than 1 then v2 points
    % inside the simplex created by the vectors v0 and v1
    if U(1,i)>=0 && U(2,i)>=0 && sum(U(:,i))<=1            % also include nodes on the boundary
        ind(i) = 1;
    end
end
P = X(:,ind==1);

%---Intersections---%
if nargout>1
    ind = zeros(1,9);
    Q = zeros(2,9);
    n = [0,0,0];
    indR = n;
    for i=1:3
        x1=U(1,i);
        x2=U(2,i);
        y1=U(1,mod(i,3)+1);
        y2=U(2,mod(i,3)+1);
        
        % with v1
        if sign(x1)~=sign(y1)       % check if intersection occurs
            y0 = y2 - y1*(x2-y2)/(x1-y1);
            if y0>=0 && y0<=1
                Q(:,i) = Y(:,1) + y0*v1;
                ind(i) = 1;
                n(i) = 1;
            end
            if exist('y0old','var') % check if Y is interior to X
                if sign(y0old)~=sign(y0)
                    indR(1) = 1;
                end
                if sign(y0old-1)~=sign(y0-1)
                    indR(3) = 1;
                end
            end
            y0old = y0;
        end
        
        % with v0
        if sign(x2)~=sign(y2)
            x0 = y1 - y2*(x1-y1)/(x2-y2);
            if x0>=0 && x0<=1
                Q(:,3+i) = Y(:,1) + x0*v0;
                ind(3+i) = 1;
                n(i) = 1;
            end
            if exist('x0old','var') % check if Y is interior to X
                if sign(x0old)~=sign(x0)
                    indR(1) = 1;
                end
                if sign(x0old-1)~=sign(x0-1)
                    indR(2) = 1;
                end
            end
            x0old = x0;
        end
        
        % with the line connecting v0 and v1
        if sign(x1+x2-1)~=sign(y1+y2-1)
            x = ( (x1-y1) + y1*(x2-y2) - y2*(x1-y1) )/( (x1-y1) + (x2-y2) );
            y = ( (x2-y2) + y2*(x1-y1) - y1*(x2-y2) )/( (x1-y1) + (x2-y2) );
            if x*y>=0
                Q(:,6+i) = Y(:,1) + x*v0 + y*v1;
                ind(6+i) = 1;
                n(i) = 1;
            end
            % new addition meant to increase robustness when nonzero
            % vertices of reference triangle may be excluded wrongly
            if exist('xold','var') && exist('yold','var')
                if sign(y-x-1)~=sign(yold-xold-1)
                    indR(3) = 1;
                end
                if sign(y-x+1)~=sign(yold-xold+1)
                    indR(2) = 1;
                end
            end
            xold = x; yold = y;
        end

    end
    Q = Q(:,ind==1);  % the intersections
    R = Y(:,indR==1); % the points of Y that are inside X
end
end

function P=SortAndRemoveDoubles(P)
% SORTANDREMOVEDOUBLES sort points and remove duplicates
%   P=SortAndRemoveDoubles(P); orders polygon corners in P counter
%   clock wise and removes duplicates

ep=10*eps;                           % tolerance for identical nodes
m=size(P,2); 
if m>0                              
    c=sum(P,2)/m;                      % order polygon corners counter 
    ao=zeros(1,m);
    for i=1:m                          % clockwise
        d=P(:,i)-c; ao(i)=angle(d(1)+1i*d(2));
    end
    [~,id]=sort(ao); 
    P=P(:,id);
    i=1;j=2;                           % remove duplicates
    while j<=m
        if norm(P(:,i)-P(:,j))>ep
            i=i+1;P(:,i)=P(:,j);j=j+1;
        else
            j=j+1;
        end
    end
    P=P(:,1:i);
end
end

function M=MortarInt(T,X,Y)
% MORTARINT computes mortar contributions 
%   M=MortarInt(T,X,Y); computes for triangles X and Y with nodal coordinates 
%   stored columnwise the integrals of all P1 element shape function combinations 
%   between triangles X and Y on triangle T.
%   The result is stored in the 3 by 3 matrix M 


Jd=-T(1,1)*T(2,3)-T(1,2)*T(2,1)+T(1,2)*T(2,3)+T(1,1)*T(2,2)+...
   T(1,3)*T(2,1)-T(1,3)*T(2,2);

a=T(1,1); 
d=T(2,1);
e=-T(2,1)+T(2,2);
b=-T(1,1)+T(1,2);
c=-T(1,1)+T(1,3);
f=-T(2,1)+T(2,3);

T1=1/2*(X(1,1)*(X(2,2)-X(2,3))-X(1,2)*(X(2,1)-X(2,3))+X(1,3)*(X(2,1)-X(2,2)));
T2=1/2*(Y(1,1)*(Y(2,2)-Y(2,3))-Y(1,2)*(Y(2,1)-Y(2,3))+Y(1,3)*(Y(2,1)-Y(2,2)));

A11=1/(2*abs(T1))*((X(2,2)-X(2,3))*b-(X(1,2)-X(1,3))*e);
B11=1/(2*abs(T1))*((X(2,2)-X(2,3))*c-(X(1,2)-X(1,3))*f);
C11=1/(2*abs(T1))*((X(2,2)-X(2,3))*a-(X(1,2)-X(1,3))*d+(X(1,2)*X(2,3)-X(2,2)*X(1,3)));
A12=1/(2*abs(T1))*((X(2,3)-X(2,1))*b+(X(1,1)-X(1,3))*e);
B12=1/(2*abs(T1))*((X(2,3)-X(2,1))*c+(X(1,1)-X(1,3))*f);
C12=1/(2*abs(T1))*((X(2,3)-X(2,1))*a+(X(1,1)-X(1,3))*d-(X(1,1)*X(2,3)-X(2,1)*X(1,3)));
A13=1/(2*abs(T1))*((X(2,1)-X(2,2))*b-(X(1,1)-X(1,2))*e);
B13=1/(2*abs(T1))*((X(2,1)-X(2,2))*c-(X(1,1)-X(1,2))*f);
C13=1/(2*abs(T1))*((X(2,1)-X(2,2))*a-(X(1,1)-X(1,2))*d+(X(1,1)*X(2,2)-X(2,1)*X(1,2)));
A21=1/(2*abs(T2))*((Y(2,2)-Y(2,3))*b-(Y(1,2)-Y(1,3))*e);
B21=1/(2*abs(T2))*((Y(2,2)-Y(2,3))*c-(Y(1,2)-Y(1,3))*f);
C21=1/(2*abs(T2))*((Y(2,2)-Y(2,3))*a-(Y(1,2)-Y(1,3))*d+(Y(1,2)*Y(2,3)-Y(2,2)*Y(1,3)));
A22=1/(2*abs(T2))*((Y(2,3)-Y(2,1))*b+(Y(1,1)-Y(1,3))*e);
B22=1/(2*abs(T2))*((Y(2,3)-Y(2,1))*c+(Y(1,1)-Y(1,3))*f);
C22=1/(2*abs(T2))*((Y(2,3)-Y(2,1))*a+(Y(1,1)-Y(1,3))*d-(Y(1,1)*Y(2,3)-Y(2,1)*Y(1,3)));
A23=1/(2*abs(T2))*((Y(2,1)-Y(2,2))*b-(Y(1,1)-Y(1,2))*e);
B23=1/(2*abs(T2))*((Y(2,1)-Y(2,2))*c-(Y(1,1)-Y(1,2))*f);
C23=1/(2*abs(T2))*((Y(2,1)-Y(2,2))*a-(Y(1,1)-Y(1,2))*d+(Y(1,1)*Y(2,2)-Y(2,1)*Y(1,2)));


M(1,1)=1/24*Jd*(2*A11*A21+B11*A21+A11*B21+2*B11*B21+4*C11*A21+4*A11*C21+4*C11*B21+4*B11*C21+12*C11*C21);
M(1,2)=1/24*Jd*(2*A11*A22+B11*A22+A11*B22+2*B11*B22+4*C11*A22+4*A11*C22+4*C11*B22+4*B11*C22+12*C11*C22);
M(1,3)=1/24*Jd*(2*A11*A23+B11*A23+A11*B23+2*B11*B23+4*C11*A23+4*A11*C23+4*C11*B23+4*B11*C23+12*C11*C23);
M(2,1)=1/24*Jd*(2*A12*A21+B12*A21+A12*B21+2*B12*B21+4*C12*A21+4*A12*C21+4*C12*B21+4*B12*C21+12*C12*C21);
M(2,2)=1/24*Jd*(2*A12*A22+B12*A22+A12*B22+2*B12*B22+4*C12*A22+4*A12*C22+4*C12*B22+4*B12*C22+12*C12*C22);
M(2,3)=1/24*Jd*(2*A12*A23+B12*A23+A12*B23+2*B12*B23+4*C12*A23+4*A12*C23+4*C12*B23+4*B12*C23+12*C12*C23);
M(3,1)=1/24*Jd*(2*A13*A21+B13*A21+A13*B21+2*B13*B21+4*C13*A21+4*A13*C21+4*C13*B21+4*B13*C21+12*C13*C21);
M(3,2)=1/24*Jd*(2*A13*A22+B13*A22+A13*B22+2*B13*B22+4*C13*A22+4*A13*C22+4*C13*B22+4*B13*C22+12*C13*C22);
M(3,3)=1/24*Jd*(2*A13*A23+B13*A23+A13*B23+2*B13*B23+4*C13*A23+4*A13*C23+4*C13*B23+4*B13*C23+12*C13*C23);
end