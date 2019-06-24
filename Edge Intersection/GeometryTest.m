function [P,Q,n] = GeometryTest(X,Y)
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
    for i=1:3
        r = U(:,mod(i,3)+1) - U(:,i);
        A = [r, zeros(2,1)];
        x1=U(1,i);
        x2=U(2,i);
        y1=U(1,mod(i,3)+1);
        y2=U(2,mod(i,3)+1);
        
        % with v1
        if sign(x1)~=sign(y1)       % check if intersection occurs
            A(:,2) = [0;-1]; b = A \ -U(:,i);
            if b(2)>=0 && b(2)<=1
                Q(:,i) = Y(:,1) + v1*(x2 + b(1)*r(2));
                ind(i) = 1;
                n(i) = 1;
            end
        end
        
        % with v0
        if sign(x2)~=sign(y2)
            A(:,2) = [-1;0]; b = A \ -U(:,i);
            if b(2)>=0 && b(2)<=1
                Q(:,3+i) = Y(:,1) + v0*(x1 + b(1)*r(1));
                ind(3+i) = 1;
                n(i) = 1;
            end
        end
        
        % with the line connecting v0 and v1
        if sign(x1+x2-1)~=sign(y1+y2-1)
            A(:,2) = [-1;1]; b = A \ ([0;1] - U(:,i));
            if b(2)>=0 && b(2)<=1
                Q(:,6+i) = Y(:,1) + [v0, v1]*(U(:,i) + b(1)*r);
                ind(6+i) = 1;
                n(i) = 1;
            end
        end

    end
    Q = Q(:,ind==1);
end
end