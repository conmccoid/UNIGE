function [P,Q,n] = Geometry(X,Y)
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
        x1=U(1,i);
        x2=U(2,i);
        y1=U(1,mod(i,3)+1);
        y2=U(2,mod(i,3)+1);
        
        x = -1; y = -1;

        % with v1
        if sign(x1)~=sign(y1)       % check if intersection occurs
            y = y2 - y1*(x2-y2)/(x1-y1);
        elseif x1==0 && max(x2,y2)>1 && min(x2,y2)<1
            % if the lines are not parallel then the only intersections
            % that can occur are coincidental lines - if the entire line is
            % coincident then the line is covered by Interior points above;
            % if not, the only intersection not covered is at the point
            % y=1, x=0.
            y = 1;
        end
        if y>=0 && y<=1
            Q(:,i) = Y(:,1) + y*v1;
            ind(i) = 1;
            n(i) = 1;
        end

        % with v0
        if sign(x2)~=sign(y2)
            x = y1 - y2*(x1-y1)/(x2-y2);
        elseif x2==0 && max(x1,y1)>1 && min(x1,y1)<1
            x = 1;
        end
        if x>=0 && x<=1
            Q(:,3+i) = Y(:,1) + x*v0;
            ind(3+i) = 1;
            n(i) = 1;
        end

        % with the line connecting v0 and v1
        if sign(x1+x2-1)~=sign(y1+y2-1)
            x = ( (x1-y1) + y1*(x2-y2) - y2*(x1-y1) )/( (x1-y1) + (x2-y2) );
            y = ( (x2-y2) + y2*(x1-y1) - y1*(x2-y2) )/( (x1-y1) + (x2-y2) );
            if x>=0 && y>=0 && x<=1 && y<=1
                Q(:,6+i) = Y(:,1) + x*v0 + y*v1;
                ind(6+i) = 1;
                n(i) = 1;
            end
        end

    end
    Q = Q(:,ind==1);
end
end