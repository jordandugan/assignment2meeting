%set the geometry of the region
L = 30;
W = 20;

%set the distance between meshpoints and determine the number points used
meshspace = 0.5;
nx = floor(L/meshspace + 1);
ny = floor(W/meshspace + 1);

G = sparse(nx*ny); %% Discretized Laplacian Operator
B = zeros(1,nx*ny); %% Mostly zeros but a few ones for BCs
for i = 1:nx
    for j = 1:ny
        n = j +(i-1)*ny; %Map 2d Geometry to 1D Vector
        
        %Left side V=1 BC
        if i == 1
            G(n,n) = 1;
            B(n) = 1;
        
        %Right side V=0 BC
        elseif i == nx 
            G(n,n) = 1;
            
        %Bottom Side Absorbing boundary condition
        elseif j == 1 
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nyp = j+1 + (i-1)*ny;
            
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nyp) = 1; 
            
            
        %Top side absorbing boundary condition
        elseif j == ny
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
        
        
        %Inner Nodes
        else
            nxm = j + (i-2)*ny;
            nxp = j + i*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;
            
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1; 
        
        end
    end
end

%Solve Matrix equation to find V
V = G\B';

%Map solution back to 2D space
Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j +(i-1)*ny;
        Vmap(i,j) = V(n);
    end
end

%Plot potential
x = 0:meshspace:L;
figure(1)
plot(x,Vmap(:,floor(ny/2)));
title('Potential vs x-position')
xlabel('x')
ylabel('Potential')
figure(2)
[X, Y] = meshgrid(0:meshspace:W,0:meshspace:L);
surf(X,Y,Vmap)
colorbar
hold on
imagesc([0 20],[0 30],Vmap)
        