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
        n = j +(i-1)*ny; %Map 2D geometry to 1D vector
        
        % V=1 @ x=0 BC
        if i == 1
            G(n,n) = 1;
            B(n) = 1;
            
        % V=1 @x=L BC
        elseif i == nx 
            G(n,n) = 1;
            B(n) = 1;
        
        %V=0 @ y=0 BC
        elseif j == 1 
           G(n,n) = 1;
        
        %V=0 @ y=W BC
        elseif j == ny
            G(n,n) = 1;
        
        %Matrix elemenets for inner nodes
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

%Plot 
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
