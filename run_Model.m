
%% Define parameters
%% So far I have assumed that rI = rR

p.dt = 0.01; %time step
p.N = 5000; % number of particles
p.M = 5000; % number of tethers
p.k = 100; %stiffness constant
p.L = 1; % domain size
p.tmax = 30; % number of timesteps
p.noise = 0.1;
p.noiseOb = 0;
p.r = 0.05; % interaction radius of SPPs to themselves
p.r1 = 0.05; % interaction radius of objects to SPPs
p.vl = 1; %constant velocity
p.eta = 1;   %friction coefficients
p.zeta = 1;
p.A = pi; %force mass
p.AR = 0.02*4/p.r1;
p.nu = 10;
p.res = 100;
p.dx = p.L/p.res;
p.xVals = linspace(-p.L/2+p.dx/2,p.L/2-p.dx/2,p.res);
p.yVals = linspace(-p.L/2+p.dx/2,p.L/2-p.dx/2,p.res);

%% Run

%initialise

%data=Model_mod4(p);

data=Model_mod_numerical_tvar(p);
save('Num','data')
data=Model_mod4_tvar(p);
save('IBM','data')
%{
for mu = 0.0125:0.0025:0.0175
    for r = [0.2 0.3]
        p.r1 = r;
        p.AR = mu*4/p.r;
        filename = 'NM4_mu='+string(mu)+'_r='+string(p.r1)+'.mat';
        Var = zeros(20,1);
        Var(1)=p.dt;
        Var(2)=p.N;
        Var(3)=p.M;
        Var(4)=p.k;
        Var(5)=p.L;
        Var(6)=p.tmax;
        Var(7)=p.noise;
        Var(8)=p.noiseOb;
        Var(9)=p.r;
        Var(10)=p.r1;
        Var(11)=p.vl;
        Var(12)=p.eta;
        Var(13)=p.zeta;
        Var(14)=p.A;
        Var(15)=p.AR;
        Var(16)=p.nu;
        Var(17)=p.res;
        Var(18)=p.dx;
        data=Model_mod_numerical_tvar(p);
        save(filename,'data')
        save(filename,'Var','-append')
    end
end
%}

%_numerical
%% Store
%{
Var = zeros(20,1);
Var(1)=p.dt;
Var(2)=p.N;
Var(3)=p.M;
Var(4)=p.k;
Var(5)=p.L;
Var(6)=p.tmax;
Var(7)=p.noise;
Var(8)=p.noiseOb;
Var(9)=p.r;
Var(10)=p.r1;
Var(11)=p.vl;
Var(12)=p.eta;
Var(13)=p.zeta;
Var(14)=p.A;
Var(15)=p.AR;
Var(16)=p.nu;
Var(17)=p.res;
Var(18)=p.dx;

filename = 'r=0.3_' + string(datetime);

save(filename,'data')
save(filename,'Var','-append')

myVideo = VideoWriter('movingclusters');
myVideo.FrameRate = 10;
open(myVideo)
%}

%
%% Plot
%{
clf;
l=0.025;

for i=1:p.tmax/0.1
    
    clf;
    hold on
    
    X1m=mod(data(i).X1+p.L/2,p.L)-p.L/2;
    X2m=mod(data(i).X2+p.L/2,p.L)-p.L/2;
    
    quiver(data(i).x-l/2*cos(data(i).thetas),data(i).y-l/2*sin(data(i).thetas),l*cos(data(i).thetas), l*sin(data(i).thetas), 'AutoScale', 'off')
    scatter(X1m,X2m,p.r1*100, 'r', 'filled')
    axis equal
    xlim([-1 1]*p.L/2)
    ylim([-1 1]*p.L/2)
    
    pause(0.01)
    %{
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    %}
end
%}

savefig(filename)

%close(myVideo)

%}
