% Script to generate anatomy_cortical_pairwise_distance_distribution.mat
parpool(2);

%%%% Set up triangular mesh of cortex (units in mm)
load('/lustre04/scratch/nbrake/data/anatomy_nyhead_model.mat');
X = struct();
X.vertices = sa.cortex75K.vc;
X.faces= sa.cortex75K.tri;

A0 = zeros(size(X.faces,1),1);
norm_vec = zeros(size(X.faces));
for i = 1:size(X.faces,1)
    tri = X.vertices(X.faces(i,:),:);
    A0(i) = 1/2*norm(cross(tri(1,:)-tri(2,:),tri(3,:)-tri(1,:))); % Area of each triangular face
    norm_vec(i,:) = mean(sa.cortex75K.normals(X.faces(i,:),:)); % Normal to face
end
X.face_norm = norm_vec./vecnorm(norm_vec,2,2);
X.face_area = A0;

%%%% Get total mesh area within a radius R of mesh point i
% Sample N starting positions
N = 10e3;
idcs = randperm(size(X.faces,1),N);
% Compute area within 20 radii on log scale, ranging from 0.01 mm to 20 mm. Beyong 20 mm, correlation function
% drops essentially to zero, so larger radii would be wasted computations.
rValues = 10.^linspace(-2,2,50);

signed_area = zeros(length(rValues),length(idcs));
total_area = zeros(length(rValues),length(idcs));
h = waitbar(0);
parfor i = 1:N
    % Compute area
    A = areaCDF(X,idcs(i),rValues);

    % Total area within radius
    total_area(:,i) = A(:,2);

    % Area of each triangular face scaled by inner product between normal vector and startpoint normal.
    % This is used to compute EEG scaling from dipole correlations.
    signed_area(:,i) = A(:,1);
end

save('/lustre04/scratch/nbrake/data/simulation_analyzed/cortical_area/pairwise_distance.mat','A','B','idcs','rValues');

function A = areaCDF(X,idcs0,rValues)
% This function takes a cortical triangular mesh, an index for a starting face, and a vector rValues
% For each radius in rValues, all the other triangular faces are identified that intersect with a ball of radius r,
% centered at the starting face. The area of all these triangular faces is computed
    c0 = mean(X.vertices(X.faces(idcs0,:),:));
    norm0 = X.face_norm(idcs0,:);
    V = X.vertices - c0; % center the entire mesh grid on starting face
    alpha = sum(X.face_norm.*norm0,2);

    included = zeros(size(X.face_area));
    A = zeros(length(rValues),2);
    for j = 1:length(rValues)
        R = rValues(j);
        idcs = setdiff(1:size(X.faces,1),find(included));
        if(j>1)
            A(j,:) = A(j-1,:);
        end
        for i = idcs
            tri = V(X.faces(i,:),:);
            if(sum(vecnorm(tri,2,2)<(R+2)))
                A_temp = calculate_intersection_area(tri,R);
                if(A_temp./X.face_area(i)==1)
                    included(i)=1;
                end
                A(j,1) = A(j,1)+A_temp*alpha(i);
                A(j,2) = A(j,2)+A_temp;
            end
        end
    end
end
function A = calculate_intersection_area(tri,R)
% Given vertices for a triangle, tri, and a raidus, R, this function
% returns the area of intersection between the triangle and a sphere
% centered at zero of radius R.

    % Entire triangle is within sphere
    if(all(vecnorm(tri,2,2)<R))
        A = 1/2*norm(cross(tri(1,:)-tri(2,:),tri(3,:)-tri(1,:)));
        return;
    end

    % Get plane of the triangle
    params = null([tri,ones(3,1)])';
    N = params(1:3);
    rho = -params(4)/norm(N,2);

    % Entire triangle is outside sphere  (type 1)
    if(abs(rho)>R)
        A = 0;
        return;
    end

    % Calculate circle of the sliced sphere
    R0 = sqrt(R^2-rho^2);
    c1 = rho*N/norm(N,2);

    % Calculate intersection points
    pairs = [1,2;2,3;3,1];
    A = sum(tri(pairs(:,1),:).^2,2) - R^2;
    C = sum((tri(pairs(:,1),:)-tri(pairs(:,2),:)).^2,2);
    B = sum(tri(pairs(:,2),:).^2,2) - A - C - R^2;
    in = zeros(6,3);
    T = zeros(6,3);
    for k = 1:3
        t = roots([C(k),B(k),A(k)]);
        in(2*k-1:2*k,:) = [tri(pairs(k,1),:).*(1-t) + tri(pairs(k,2),:).*t];
        if(sum(t<0) == 2 || sum(t>1) == 2)
            % The line intersects, but segment doesn't
            % Combine this case with imaginary roots
            t = t*sqrt(-1);
        end
        T(2*k-1:2*k,:) = [pairs(k,:),t(1);pairs(k,:),t(2)];
    end

    % No intersection points
    if(all(imag(T(:,3))~=0))
        if(norm(mean(tri),2)<R)
            % Entire sphere is inside triangle
            A = pi*R0^2;
        else
            % Entire triangle is outside sphere (type 2)
            A = 0;
        end
        return
    end

    % Compute vertices that define intersected area
    in(T(:,3)<0,:) = tri(T(T(:,3)<0,1),:);
    in(T(:,3)>1,:) = tri(T(T(:,3)>1,2),:);
    in(imag(T(:,3))~=0,:) = [];
    T(imag(T(:,3))~=0,:) = [];
    [T,idcs] = sortrows(T);
    in = in(idcs,:);

    % Change coordinates to intersecting plane
    [x1,x2] = getOrthBasis(params(1:3));
    in_proj = [sum((in-c1).*x1,2),sum((in-c1).*x2,2)];
    in_proj(end+1,:) = in_proj(1,:);

    % Get area of polygon
    A = sum((in_proj(2:end,1)+in_proj(1:end-1,1)).*(diff(in_proj(:,2))))/2;

    % Get sign of area to determine direction of circular segments
    clockwise = (A<0);

    % Convert to polar coordinates
    theta = mod(atan2(in_proj(:,2),in_proj(:,1)),2*pi);

    % Angles that define start and end points of each circular segments
    arcTheta = [];
    for k = 2:2:length(theta)
        if(theta(k)~=theta(k+1))
            if(clockwise)
                if(theta(k+1)>theta(k))
                    theta(k+1) = theta(k+1)-2*pi;
                end
            else
                if(theta(k+1)<theta(k))
                    theta(k+1) = theta(k+1)+2*pi;
                end
            end
            arcTheta(end+1,:) = theta(k:k+1);
        end
    end
    % If entire intersection is one circular segment
    if(size(arcTheta,1)==1)
        arcTheta = mod(arcTheta,2*pi);
        alph1 = mod(0.5*(arcTheta(1)+arcTheta(2)),2*pi);
        alph2 = mod(alph1-pi,2*pi);
        tri_proj = [sum((tri-c1).*x1,2),sum((tri-c1).*x2,2)];
        mu = mean(tri_proj);
        N1 = norm([cos(alph1)*R0,sin(alph1)*R0] - mu,2);
        N2 = norm([cos(alph2)*R0,sin(alph2)*R0] - mu,2);
        if(N1<N2)
            alph0 = alph1;
        else
            alph0 = alph2;
        end
        arcTheta = sort(arcTheta);
        if(alph0<arcTheta(1) || alph0>arcTheta(2))
            arcTheta(2) = arcTheta(2)-2*pi;
        end
    end
    % Add area of circular segments to that of polygon
    A = abs(A);
    if(~isempty(arcTheta))
        alph = abs(arcTheta(:,1)-arcTheta(:,2));
        A = A+sum(R0^2/2*(alph-sin(alph)));
    end
end
function [x1,x2] = getOrthBasis(x0)
    x0 = x0/norm(x0,2);
    if(prod(size(x0))==3)
        if(x0>0.9)
            x1 = [0,1,0];
        else
            x1 = [1,0,0];
        end
        x1 = x1-x0*sum(x1.*x0);
        x1 = x1/norm(x1,2);
        x2 = cross(x1,x0);
    else
        idcs = x0(:,1)>0.9;
        x1 = zeros(size(x0));
        x1(idcs,2) = 1;
        x1(~idcs,1) = 1;

        x1 = x1-x0.*sum(x1.*x0,2);
        x1 = x1./vecnorm(x1,2,2);
        x2 = cross(x1,x0);
    end
end