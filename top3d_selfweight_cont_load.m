%for fcc: x = 30, y = 9 , z = 18
% function top3d_selfweight_cont_load(nelx,nely,nelz,volfrac,penal,rmin)
% USER-DEFINED LOOP PARAMETERS
nelx=30;
nely=12;
nelz=24;
volfrac=0.5;
penal=3;
rmin=1.2;
disp 'time'
maxloop = 100;    % Maximum number of iterations
tolx = 0.01;      % Terminarion criterion
displayflag = 0;  % Display structure flag
% USER-DEFINED MATERIAL PROPERTIES
E0 = 70000000000;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
nu = 0.3;         % Poisson's ratio
% % USER-DEFINED LOAD DOFs
% [il,jl,kl] = meshgrid(nelx/2, 0, 0:nelz);                 % Coordinates
% loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
% loaddof = 3*loadnid(:) - 1;                             % DOFs

% USER-DEFINED LOAD DOFs (continuous load on top and bottom planes)
[il_top, jl_top, kl_top] = meshgrid(0:nelx, 0:nely, nelz);  % Top plane (z = nelz)
[il_bottom, jl_bottom, kl_bottom] = meshgrid(0:nelx, 0:nely, 0);  % Bottom plane (z = 0)

% Calculate node IDs for top and bottom planes
loadnid_top = kl_top * (nelx + 1) * (nely + 1) + il_top * (nely + 1) + (nely + 1 - jl_top);
loadnid_bottom = kl_bottom * (nelx + 1) * (nely + 1) + il_bottom * (nely + 1) + (nely + 1 - jl_bottom);

% Define load DOFs (for example, apply loads in the z-direction, which is typically the third DOF)
loaddof_top = 3 * loadnid_top(:);         % z-direction DOFs for the top plane
loaddof_bottom = 3 * loadnid_bottom(:);    % z-direction DOFs for the bottom plane
loaddof = unique([loaddof_top; loaddof_bottom]);  % Combine and remove duplicates

% Define the force vector with continuous load applied to top and bottom planes

% % USER-DEFINED SUPPORT FIXED DOFs
% [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
% fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
% fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs

% USER-DEFINED SUPPORT FIXED DOFs (both end faces parallel to y-z plane)
[iif_left, jf_left, kf_left] = meshgrid(0, 0:nely, 0:nelz);  % Left face (x = 0)
[iif_right, jf_right, kf_right] = meshgrid(nelx, 0:nely, 0:nelz);  % Right face (x = nelx)

% Calculate node IDs for both faces
fixednid_left = kf_left * (nelx + 1) * (nely + 1) + iif_left * (nely + 1) + (nely + 1 - jf_left);
fixednid_right = kf_right * (nelx + 1) * (nely + 1) + iif_right * (nely + 1) + (nely + 1 - jf_right);

% Combine fixed DOFs from both faces
fixeddof_left = [3 * fixednid_left(:); 3 * fixednid_left(:) - 1; 3 * fixednid_left(:) - 2];
fixeddof_right = [3 * fixednid_right(:); 3 * fixednid_right(:) - 1; 3 * fixednid_right(:) - 2];
fixeddof = unique([fixeddof_left; fixeddof_right]);  % Combine and remove duplicates

% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-70*9.82,ndof,1);
% F = sparse(loaddof, 1, -1 / numel(loaddof), ndof, 1);  % Distribute load uniformly
% USER-DEFINED MATERIAL PROPERTIES (including density for self-weight calculation)
density = 2700;     % Material density
g = 9.81;        % Gravitational acceleration (m/s^2)

% Initialize the force vector for self-weight (same size as F)
F_self_weight = sparse(ndof, 1);

% Loop through each element to compute the self-weight contribution
for k = 2:nelz-1  % Ignore surface elements in the z-direction
    for i = 2:nelx-1  % Ignore surface elements in the x-direction
        for j = 2:nely-1  % Ignore surface elements in the y-direction
            % Element volume (assuming unit cube elements)
            element_volume = 1 * 1 * 1;  % Modify if element size is different

            % Calculate the element mass
            element_mass = 3*0.9*1.7*density/nelx/nely/nelz;

            % Calculate the self-weight force on each node of the element
            node_force = element_mass * g / 8;  % Divide by 8 for each node

            % Get the global node indices for this element
            element_nodes = [
                (k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j,   % Node 1
                (k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j,       % Node 2
                (k-1)*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1, % Node 3
                (k-1)*(nelx+1)*(nely+1) + i*(nely+1) + j+1,     % Node 4
                k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j,       % Node 5
                k*(nelx+1)*(nely+1) + i*(nely+1) + j,           % Node 6
                k*(nelx+1)*(nely+1) + (i-1)*(nely+1) + j+1,     % Node 7
                k*(nelx+1)*(nely+1) + i*(nely+1) + j+1          % Node 8
            ];

            % Apply self-weight load to each node (in the z-direction, DOF 3)
            for n = 1:8
                dof_z = 3 * element_nodes(n);  % Z DOF for each node
                F_self_weight(dof_z) = F_self_weight(dof_z) - node_force;  % Negative for downward force
            end
        end
    end
end

% Combine self-weight load with the existing load vector
F = F + F_self_weight;

U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
% INITIALIZE ITERATION
x = repmat(volfrac,[nely,nelx,nelz]);
% Dimensions of the 3D matrix x
[nely, nelx, nelz] = size(x);

% Initialize the new matrix with zeros
passive = zeros(nely, nelx, nelz);

% Define the size of the square zone
square_size_y = floor(nelz / 3);
square_size_z = floor(nelz / 3);

% Calculate the starting indices to center the square zone on each plane
start_y = floor((nely - square_size_y) / 2) + 1;
start_z = floor((nelz - square_size_z) / 2) + 1;

% Define the indices for the square region on the leftmost plane (x = 1)
passive(start_y:start_y+square_size_y-1, 1, start_z:start_z+square_size_z-1) = 1;

% Define the indices for the square region on the rightmost plane (x = nelx)
passive(start_y:start_y+square_size_y-1, nelx, start_z:start_z+square_size_z-1) = 1;

x(find(passive)) = 1;
xPhys = x; 
loop = 0; 
change = 1;
% START ITERATION
while change > tolx && loop < maxloop
    loop = loop+1;
    % FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx,nelz);
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    dc(:) = H*(dc(:)./Hs);  
    dv(:) = H*(dv(:)./Hs);
    % OPTIMALITY CRITERIA UPDATE
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        x(find(passive)) = 1;
        xPhys(:) = (H*xnew(:))./Hs;
        if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    % PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
    % PLOT DENSITIES
    if displayflag, clf; display_3D(xPhys); end %#ok<UNRCH>
end
clf; 
% xPhys = xPhys >= 0.75;
display_3D(xPhys);
max(max(U))
%Top3dSTL_v3
% end


% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];

K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end
% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end

