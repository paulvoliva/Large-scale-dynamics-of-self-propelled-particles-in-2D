function data=Model_mod_numerical_tvar(p)

%This function takes in the parameters of the class p and then runs
%throught the following implementation, returning a data with structure
%elements data.x and data.y for the x and y coordinates of the SPPs,
%data.thetas to store the SPPs' orientation, data.X1 and data.X2 for the x
%and y coordinates of the obstacles as well as their fixed anchor points
%data.Y1 and data.Y2. This is for each data entry which corresponds to each
%frame of the simulation.

%% Parameters

% if no parameters are provided, use these 
if nargin<1
    p.dt = 1; %time step
    p.N = 10; % number of particles
    p.M = 10; % number of tethers
    p.k = 1; %stiffness constant
    p.L = 1; % domain size
    p.tmax = 1; % number of timesteps
    p.noise = 0.3;
    p.noiseOb = 0.2;
    p.r = 0.2; % interaction radius of SPPs to themselves
    p.r1 = 0.2; % interaction radius of objects to SPPs
    p.vl = 1; %constant velocity
    p.eta = 1;   %friction coefficients
    p.zeta = 1;
    p.A = pi; %force mass
    p.AR = 0.02;
    p.nu = 1;
    p.res = 100;
    p.dx = p.L/p.res;
    p.xVals = linspace(-p.L/2+p.dx/2,p.L/2-p.dx/2,p.res);
    p.yVals = linspace(-p.L/2+p.dx/2,p.L/2-p.dx/2,p.res);
end

%% Initialize
%
% initialize results structure
data(p.tmax+1).x=[];
data(p.tmax+1).y=[];
data(p.tmax+1).thetas=[];
data(p.tmax+1).X1=[];   %Object positions
data(p.tmax+1).X2=[];
data(p.tmax+1).Y1=[];   %Tether positions
data(p.tmax+1).Y2=[];
data(p.tmax+1).dt=[];

% random initial positions
x=1*(p.L*rand(p.M,1)-p.L/2);
y=1*(p.L*rand(p.M,1)-p.L/2);
X1=1*(p.L*rand(p.M,1)-p.L/2);
X2=1*(p.L*rand(p.M,1)-p.L/2);
Y1=X1;
Y2=X2;
%}
% random initial orientations
thetas=2*pi*rand(p.N,1);

%
% save initial conditions
data(1).x=x;
data(1).y=y;
data(1).thetas=thetas;
data(1).X1=X1;
data(1).X2=X2;
data(1).Y1=X1;
data(1).Y2=X2;

i = 1;
totaldt = 0;

%% Run

% look through timesteps
while i-1<p.tmax/0.05
    
    tic
    
    % update
    [x,y,thetas,X1,X2,dt]=update(x,y,thetas,X1,X2,Y1,Y2,p);
    
    totaldt = totaldt+dt;
    
    % store
    if totaldt > 0.05
        data(i).x=x;
        data(i).y=y;
        data(i).thetas=thetas;
        data(i).X1=X1;
        data(i).X2=X2;
        data(i).Y1=Y1;
        data(i).Y2=Y2;
        data(i).dt=dt;
        i = i+1;
        totaldt = 0;
    end
    
    % some reporting
    fprintf('Time step %s took %s sec\n', num2str(i), num2str(toc))
    
end

end

%% Ancillary function

function [x,y,thetas,X1,X2,dt]=update(x,y,thetas,X1,X2,Y1,Y2,p)
%{
update the positions and directions
%}
    X1m=mod(X1+p.L/2,p.L)-p.L/2;
    X2m=mod(X2+p.L/2,p.L)-p.L/2;
    
    denss = density(x,y,p);    %swimmer density
    denso = density(X1m,X2m,p);  %object density
    [NPh,X_,Y_] = NPhi(p.r1,p);    %Phi evaluated on a p.dx by p.dx grid
    [NPs,X__,Y__] = NPhi(p.r,p);   %Psi evaluated on a p.dx by p.dx grid
    [anx,any] = convrp(X1m,X2m,denss,NPh,X_,Y_,p);    %object-swimmer interaction for object
    [anx1,any1] = convrp(x,y,denso,NPh,X_,Y_,p);    %object-swimmer interaction for swimmer
    [anx2,any2] = convrp(x,y,denss,NPs,X__,Y__,p);  %swimmer-swimmer interaction
    
    %calculate all the forces
    fx = -1/p.zeta*(p.A*anx1+p.AR*anx2);
    fy = -1/p.zeta*(p.A*any1+p.AR*any2);
    Fx = 1/p.eta*(p.k*(X1-Y1)+p.A*anx);
    Fy = 1/p.eta*(p.k*(X2-Y2)+p.A*any);
    
    %determine the timestep to be taken
    m = max(abs([fx fy Fx Fy]),[],'all');
    
    dt = p.r/(2*m);
    
    %we limit the highest possible timestep
    if dt > 0.075
        dt = 0.075;
    end

    % find neighbors of all particles
    [nb, indexExt]=neighbours(x,y,x,y,p.r,p);

    % determine average angle
    avrg_angle=avrg_anglev(thetas,indexExt,nb,p);
    
    %aux vectors
    B_X=cos(thetas)+p.nu*dt/2*(cos(avrg_angle)-cos(thetas));
    B_Y=sin(thetas)+p.nu*dt/2*(sin(avrg_angle)-sin(thetas));
    
    %the signed angle between B and omega
    angle = -atan2d(B_X.*sin(thetas)-B_Y.*cos(thetas),B_X.*cos(thetas)+B_Y.*sin(thetas))*pi/180;
    
    % add noise
    avrg_angle_noise=thetas+2*angle+sqrt(2*p.noise*dt).*randn(p.N,1);
    
    % update positions & direction
    x=x+dt*(p.vl*cos(avrg_angle_noise)+fx);
    y=y+dt*(p.vl*sin(avrg_angle_noise)+fy);
    thetas = avrg_angle_noise;
    X1=X1-dt*Fx+sqrt(2*dt*p.noiseOb)*randn(p.M,1);
    X2=X2-dt*Fy+sqrt(2*dt*p.noiseOb)*randn(p.M,1);
    
    % account for periodicity
    x=mod(x+p.L/2,p.L)-p.L/2;
    y=mod(y+p.L/2,p.L)-p.L/2;
    
end

function [xExt, yExt, indexExt]=do_periodic_ext(x, y,p)
%{
make copies of particles to account for periodocity
%}
    border=p.r;
    
    % who is near a boundary
    isNearBoundary_L = x<=-p.L/2+border; isNearBoundary_R = x>=p.L/2-border;
    isNearBoundary_D = y<=-p.L/2+border; isNearBoundary_U = y>=p.L/2-border;
    
    % who is near an edge
    isNearEdge_LU=isNearBoundary_L & isNearBoundary_U; isNearEdge_LD=isNearBoundary_L & isNearBoundary_D;
    isNearEdge_RU=isNearBoundary_R & isNearBoundary_U; isNearEdge_RD=isNearBoundary_R & isNearBoundary_D;
   
    % extend
    xExt=[x; x(isNearBoundary_L)+p.L; x(isNearBoundary_R)-p.L; x(isNearBoundary_D); x(isNearBoundary_U);...
            x(isNearEdge_LU)+p.L; x(isNearEdge_LD)+p.L; x(isNearEdge_RU)-p.L; x(isNearEdge_RD)-p.L];
    yExt=[y; y(isNearBoundary_L); y(isNearBoundary_R); y(isNearBoundary_D)+p.L; y(isNearBoundary_U)-p.L;...
            y(isNearEdge_LU)-p.L; y(isNearEdge_LD)+p.L; y(isNearEdge_RU)-p.L; y(isNearEdge_RD)+p.L];
        
    % copy also index infos to identify correct neighbors later
    index=(1:length(x))';    
    indexExt=[index; index(isNearBoundary_L); index(isNearBoundary_R); index(isNearBoundary_D); index(isNearBoundary_U);...
            index(isNearEdge_LU); index(isNearEdge_LD); index(isNearEdge_RU); index(isNearEdge_RD)];
        
end

function [nb, indexExt] = neighbours(x,y,X1,X2,r,p)
%{
Find neighbours within radius p.r
%}

    % extend to account for periodocity
    [xExt, yExt, indexExt]=do_periodic_ext(x, y,p);    
    
    % find neighbors
    nb = rangesearch([xExt yExt],[X1 X2],r);
 
end

%the function phi evaluated on a grid
function [R,X,Y] = NPhi(r,p)
    [X,Y] = meshgrid(p.xVals);
    R = phi(X,Y,r,p);
end

function R = phi(x,y,r,p)
    absv = sqrt(x.^2+y.^2);
    axis equal tight
    xlim([-1 1]*p.L/2); ylim([-1 1]*p.L/2)
    con = 3/(2*r^3*pi);
    I = absv < r;
    R = con*(r-absv).^2.*I;
end

%the convolution of the functions
function [anx,any] = convrp(x,y,dens,NPh,X0,Y0,p)
    [convx,convy] = gradient(convol(NPh,dens,p),p.dx);
    %We extend the grid with the following values
    lhx = interp2(X0,Y0,convx,-p.L/2+p.dx,p.yVals,'linear',0);
    rhx = interp2(X0,Y0,convx,p.L/2-p.dx,p.yVals,'linear',0);
    uhx = interp2(X0,Y0,convx,p.xVals,p.L/2-p.dx,'linear',0);
    bhx = interp2(X0,Y0,convx,p.xVals,-p.L/2+p.dx,'linear',0);
    lhy = interp2(X0,Y0,convy,-p.L/2+p.dx,p.yVals,'linear',0);
    rhy = interp2(X0,Y0,convy,p.L/2-p.dx,p.yVals,'linear',0);
    uhy = interp2(X0,Y0,convy,p.xVals,p.L/2-p.dx,'linear',0);
    bhy = interp2(X0,Y0,convy,p.xVals,-p.L/2+p.dx,'linear',0);
    [X,Y] = meshgrid(linspace(-p.L/2,p.L/2,p.res+2));
    convx = [(bhx(1)+rhx(1))/2 bhx (bhx(end)+lhx(1))/2;rhx convx lhx;(uhx(1)+rhx(1))/2 uhx (uhx(end)+lhx(1))/2];
    convy = [(bhy(1)+rhy(1))/2 bhy (bhy(end)+lhy(1))/2;rhy convy lhy;(uhy(1)+rhy(1))/2 uhy (uhy(end)+lhy(1))/2];
    anx = interp2(X,Y,convx,x,y,'linear',0);
    any = interp2(X,Y,convy,x,y,'linear',0);
end

%function to determine the weighted density
function rho = density(x,y,p)
    Xedges=[p.xVals p.L/2+p.dx/2]-p.dx/2;
    Yedges=[p.yVals p.L/2+p.dx/2]-p.dx/2;
    N = length(x);
    [r,~,~] = histcounts2(x,y,Xedges,Yedges);
    rho = r'/(N*p.dx^2);
end

function Y = convol(z,v,p)
    a = fft2(ifftshift(z));
    b = fft2(ifftshift(v));
    Y = p.dx^2*fftshift(ifft2(a.*b, 'symmetric'));
end

function avrg_angle = avrg_anglev(thetas,indexExt,nb, p)
%{
Find mean direction
%}

% intialize
avrg_angle=zeros(p.N,1);

    % loop through particles
    for k=1:p.N

        % neighbor indices
        iNb=indexExt(nb{k});
        
        % find mean direction by adding vectors
        xs=sum(cos(thetas(iNb)));
        ys=sum(sin(thetas(iNb)));
        
        % determine corresponding direction
        theta = atan2(ys,xs);
        
        % store
        avrg_angle(k) = theta;
        
    end
end
