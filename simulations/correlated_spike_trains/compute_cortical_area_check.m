% Script to generate anatomy_cortical_pairwise_distance_distribution.mat

%%%% Set up triangular mesh of cortex (units in mm)
[sa,X] = network_simulation_beluga.getHeadModel;
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
X.total_face_area = A0;
X.captured_face_area = 0*A0;

%%%% Get total mesh area within a radius R of mesh point i
% Sample N starting positions
idcs0 = randperm(size(X.faces,1),1);
% Compute area within 20 radii on log scale, ranging from 0.01 mm to 20 mm. Beyong 20 mm, correlation function
% drops essentially to zero, so larger radii would be wasted computations.
rValues = linspace(0,500,10);

% This function takes a cortical triangular mesh, an index for a starting face, and a vector rValues
% For each radius in rValues, all the other triangular faces are identified that intersect with a ball of radius r,
% centered at the starting face. The area of all these triangular faces is computed
c0 = mean(X.vertices(X.faces(idcs0,:),:));
norm0 = X.face_norm(idcs0,:);
V = X.vertices - c0; % center the entire mesh grid on starting face
alpha = sum(X.face_norm.*norm0,2);

included = zeros(1,size(X.total_face_area,1));
A = zeros(length(rValues),2);
for j = 1:length(rValues)
    R = rValues(j)
    idcs = find(~included);
    if(j>1)
        A(j,:) = sum(X.total_face_area(find(included)));
    end
    for i = idcs
        tri = V(X.faces(i,:),:);
        if(sum(vecnorm(tri,2,2)<(R+2)))
            % vecnorm(tri,2,2)
            A_temp = calculate_intersection_area(tri,R);
            % If entire face has been computed, ignore for future balls
            if(abs(A_temp - X.total_face_area(i))<1e-3)
                included(i)=1;
            end
            A(j,1) = A(j,1)+A_temp*alpha(i);
            A(j,2) = A(j,2)+A_temp;
        end
    end
end