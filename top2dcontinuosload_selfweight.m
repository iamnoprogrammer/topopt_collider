%This code is developed from the 99 line code written by Segmund, 2001. (https://rdcu.be/dwLxi)
%This code produces the 2D optimal topology as output for particular volume
%fraction. Starting domain is rectangular. Loads are continuous loads in
%top and bottom boundaries and self weight for each element. Material is
%aluminum. The entire geomtry can be tilted at any angle.
% INITIALIZE
clear all;
nelx=100;%elements in x direction/horizontal
nely=30;%elements in y direction/vertical
volfrac=0.5;%volume fraction to be filled
penal=3;%penalty factor for topology optimisation
rmin=1.3;%filter radius for topology optimisation
x(1:nely,1:nelx) = volfrac;%initial topology
loop = 0;
change = 1.;

% START ITERATION
while change > 0.1 && loop < 100
    % LOOP-PARAMETERS
    loop = loop + 1;
    xold = x;
    % FE-ANALYSIS
    [U]=FE(nelx,nely,x,penal);
    strain_energy_density(nelx,nely,U,x,penal);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [c,dc]=sensitivity_analysis(nelx,nely,x,penal,U);
    % FILTERING OF SENSITIVITIES
    [dc]   = check(nelx,nely,rmin,x,dc);
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [x]    = OC(nelx,nely,x,volfrac,dc);
    % PRINT RESULTS
%     change = max(max(abs(x-xold)));
     %[change]=plot_results(x,xold,loop,c,nelx,nely);
    figure(1)
colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end
[change]=plot_results(x,xold,loop,c,nelx,nely);
 figure(1)
colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
    lmid = 0.5*(l2+l1);
    xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
    if sum(sum(xnew)) - volfrac*nelx*nely > 0
        l1 = lmid;
    else
        l2 = lmid;
    end
end
end

%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
    for j = 1:nely
        sum=0.0;
        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
            for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
                fac = rmin-sqrt((i-k)^2+(j-l)^2);
                sum = sum+max(0,fac);
                dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
            end
        end
        dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
    end
end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U]=FE(nelx,nely,x,penal)
[KE] = lk;
[K]=assemble(KE,nelx,nely,x,penal);
[F,freedofs,fixeddofs]=load_boundary(nelx,nely,x);
[U]=solve_FE(K,F,freedofs,fixeddofs,nelx,nely);
end
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
E = 70000000000;%material's Young's modulus in Pa
nu = 0.3;%material's Poisson ratio
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
    -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
%%%%%%%%%% ASSEMBLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K]=assemble(KE,nelx,nely,x,penal)
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));  %use zeros or sparse

for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
    end
end
end
%%%%%%%%%% APPLY LOADS AND BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%
function [F,freedofs,fixeddofs]=load_boundary(nelx,nely,x)
rho=2700;%material density in kg/m^3
elemass=3*0.9*1.7*rho/nelx/nely;%material based element mass normalised for the number of elements in the domain and calculated using FCC domain size
F = sparse(2*(nely+1)*(nelx+1),1);
exmass=70;%mass of the tubes in top and bottom each in kg
fload=exmass*9.82/(nelx+1-2);%continuous load per element in top and bottom boundaries
angle=90;%degrees % angle of tilt of the geomtry from horizontal
% F(2*(nelx/2)*(nely+1)+2*(floor(nely)+1),1) = -fload;
% F(2*(nelx/2)*(nely+1)+2*(floor(nely)+1),1) = -1;
 F(2*(nely+1):2*(nely+1):2*(nelx)*(nely+1),1)=-fload*sin(angle*3.14/180);%continuous load bottom vertical component
 F(2*(nely+1)+2:2*(nely+1):2*(nelx)*(nely+1)+2,1)=-fload*sin(90*3.14/180);%continuous load top vertical component
 F(2*(nely+1)-1:2*(nely+1):2*(nelx)*(nely+1)-1,1)=-fload*cos(angle*3.14/180);%continuous load bottom horizontal component
 F(2*(nely+1)+1:2*(nely+1):2*(nelx)*(nely+1)+1,1)=-fload*cos(angle*3.14/180);%continuous load top horizontal component
counter=2*(nely+1)+4;
%self weight as average of quarter weight of surrounding cells in each node
for cx = 2:nelx-1
    for cy = 2:nely-1
        F(counter,1)=F(counter,1)+0.25*(x(cy,cx-1)+x(cy-1,cx-1)+x(cy,cx)+x(cy-1,cx))*9.82*elemass;%self weight of surrounding elements
      counter=counter+2;
    end
end

fixeddofs   = [1:2*(nely+1) 2*nelx*(nely+1)+1:2*(nelx+1)*(nely+1)];%fixed nodes in left and right boundaries
alldofs     = 1:2*(nely+1)*(nelx+1);
freedofs    = setdiff(alldofs,fixeddofs);
end
%%%%%%%%%%% SOLVE FINITE ELEMENT ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%
function [U]=solve_FE(K,F,freedofs,fixeddofs,nelx,nely)
U = sparse(2*(nely+1)*(nelx+1),1);

U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs,:)= 0;
end
%%%%%%%%%%% DISPLAY OUTPUT AND RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%
function [change]=plot_results(x,xold,loop,c,nelx,nely)
% PRINT RESULTS
change = max(max(abs(x-xold)));
disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
    ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
    ' ch.: ' sprintf('%6.3f',change )])
% PLOT DENSITIES
% figure(1)
% colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end
%%%%%%%%%% SENSITIVITY ANALYIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,dc]=sensitivity_analysis(nelx,nely,x,penal,U)
[KE] = lk;
dc=zeros(nelx,nely);
c = 0.;
for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
        c = c + x(ely,elx)^penal*Ue'*KE*Ue;
        dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
    end
end
end
%%%%%%%%%%% CALCULATION OF STRAIN ENERGY DENSITY %%%%%%%%%%%%%%%%
function strain_energy_density(nelx,nely,U,x,penal)
W=zeros(nelx,nely);
S=W;
B=1/4*[-1 0 1 0 1 0 -1 0; 0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
D=[1 .5 0;.5 1 0;0 0 0.25];  %for plain stress

%you might recycle the part from the assemble routine by copying it here
for elx=1:nelx
    for ely=1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        W(ely,elx)=0.5*U(edof)'*B'*D*B*U(edof)*penal*x(ely,elx)^(penal-1);
        S_el=D*B*U(edof)*penal*x(ely,elx)^(penal-1);
        S(ely,elx)=sqrt(S_el(1)^2+S_el(2)^2-S_el(1)*S_el(2)+3*S_el(3)^2);
    end
end
end
