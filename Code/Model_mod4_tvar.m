function data=Model_mod4_tvar(p)

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
    p.dt = 0.1; %time step
    p.N = 5000; % number of particles
    p.M = 5000; % number of tethers
    p.k = 1; %stiffness constant
    p.L = 10; % domain size
    p.tmax = 10; % number of timesteps
    p.noise = 0.3;
    p.noiseOb = 0.2;
    p.r = 0.5; % interaction radius of SPPs to themselves
    p.r1 = 0.2; % interaction radius of objects to SPPs
    p.vl = 1; %constant velocity
    p.eta = 1;   %friction coefficients
    p.zeta = 1;
    p.A = 1; %force mass
    p.AR = 0.02;
    p.nu = 1;
end

%% Initialize

% initialize results structure
data(p.tmax+1).x=[];
data(p.tmax+1).y=[];
data(p.tmax+1).thetas=[];
data(p.tmax+1).X1=[];   %Object positions
data(p.tmax+1).X2=[];
data(p.tmax+1).Y1=[];   %Tether positions
data(p.tmax+1).Y2=[];

% random initial positions
x=p.L*rand(p.N,1)-p.L/2;
y=p.L*rand(p.N,1)-p.L/2;
X1=p.L*rand(p.M,1)-p.L/2;
X2=p.L*rand(p.M,1)-p.L/2;
Y1=X1;
Y2=X2;

% random initial orientations
thetas=2*pi*rand(p.N,1);

% save initial conditions
data(1).x=x;
data(1).y=y;
data(1).thetas=thetas;
data(1).X1=X1;
data(1).X2=X2;
data(1).Y1=X1;
data(1).Y2=X2;

%% Run

i = 1;
%variable to keep track of timestep
totaldt = 0;

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

    % find neighbors of all particles
    [nb, indexExt]=neighbours(x,y,x,y,p.r,p);
    
    % find objects that are neighbours
    [NB1, index1] = neighbours(X1m,X2m,x,y,p.r1,p);

    % determine average angle
    avrg_angle=avrg_anglev(thetas,indexExt,nb,p);
    
    %aux vectors
    B_X=cos(thetas)+p.nu*p.dt/2*(cos(avrg_angle)-cos(thetas));
    B_Y=sin(thetas)+p.nu*p.dt/2*(sin(avrg_angle)-sin(thetas));
    
    %the signed angle between B and omega
    angle = -atan2d(B_X.*sin(thetas)-B_Y.*cos(thetas),B_X.*cos(thetas)+B_Y.*sin(thetas))*pi/180;
    
    % add noise
    avrg_angle_noise=thetas+2*angle+sqrt(2*p.noise*p.dt).*randn(p.N,1);
    
    %evaluate all of the sums in the ODEs
    [sumsw1, sumob1] = sumF(x,X1m,y,X2m,index1,NB1,p);
    [sumsw2, sumob2] = sumF(y,X2m,x,X1m,index1,NB1,p);
    [sum1, ~] = sumF(x,x,y,y,indexExt,nb,p);
    [sum2, ~] = sumF(y,y,x,x,indexExt,nb,p);
    
    %determine the overall forces acting on the SPPs and obstacles
    fx = 1/p.zeta*(p.A*sumsw1+p.AR*p.N/p.M*sum1);
    fy = 1/p.zeta*(p.A*sumsw2+p.AR*p.N/p.M*sum2);
    Fx = 1/p.eta*(p.k*(X1-Y1)+p.A*sumob1);
    Fy = 1/p.eta*(p.k*(X2-Y2)+p.A*sumob2);
    
    %determine the timestep to be taken
    m = max(abs([fx fy Fx Fy]),[],'all');

    dt = p.r/(2*m);

    %we limit the highest possible timestep
    if dt > 0.075
        dt = 0.075;
    end
    
    % update positions & direction
    x=x+dt*p.vl*cos(avrg_angle_noise)-dt*fx;
    y=y+dt*p.vl*sin(avrg_angle_noise)-dt*fy;
    thetas = avrg_angle_noise;
    X1=X1-dt*Fx+sqrt(2*dt*p.noiseOb)*randn(p.M,1);
    X2=X2-dt*Fy+sqrt(2*dt*p.noiseOb)*randn(p.M,1);
    
    % account for periodicity
    x=mod(x+p.L/2,p.L)-p.L/2;
    y=mod(y+p.L/2,p.L)-p.L/2;
    
    %fprintf('diff %s\n', sum(sumsw1)+sum(sumob1))
    
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

function [sumsw, sumob] = sumF(x,X1,y,X2,index,NB,p)

N=length(x);
M=length(X1);
sumsw=zeros(N,1);
sumob=zeros(M,1);

    for k=1:N
        iNB=index(NB{k});
        for i = 1:length(iNB)
            if X1(iNB(i)) == x(k) && X2(iNB(i)) == y(k) % to avoid singularities in dphi
            else
                P = dphi(x(k),X1(iNB(i)),y(k),X2(iNB(i)),p);
                sumsw(k) = sumsw(k)+P;
                sumob(iNB(i)) = sumob(iNB(i))-P;
            end
        end
    end
    
    sumsw=1/N*sumsw;
    sumob=1/N*sumob;

end

function r = dphi(x,X1,y,X2,p)

if abs(x-X1) > p.r1
    X1 = X1-sign(X1)*p.L;
elseif abs(y-X2) > p.r1
    X2 = X2-sign(X2)*p.L;
end

absv=sqrt((x-X1)^2+(y-X2)^2);
con=3/((p.r1^3)*pi);

r=con*(p.r1-absv)*(X1-x)/absv;

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
