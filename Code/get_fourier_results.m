function [pattNr, pattDir, score]=get_fourier_results(xPos, yPos,nr)

%This function takes in the x-coordinates of the SPPs 'xPos', the
%y-coordinates 'yPos' and the number of the peak we want 'nr' and returns
%the pattern number 'pattNr', its direction 'pattDir' and the score it
%assigns to this results' accuracy 'score'.

isPlot=1;

N=length(xPos);

%% Define radius and direction values

nn=200;
rad=linspace(1.5,20,nn);
% phi=linspace(0,pi/2,nn);
phi=linspace(0,pi,nn);

[R,P]=meshgrid(rad,phi);

%% Calculate Fourier coefficient
fourierCoeff=0;

for i_part=1:N
    fourierCoeff=fourierCoeff+exp(-2*pi*1i*(R.*cos(P)*xPos(i_part)+R.*sin(P)*yPos(i_part)));
end

fourierCoeff=fourierCoeff/N;

% determine absolute value
absFourier=abs(fourierCoeff);

%% Find Peaks

peakIndices=FastPeakFind(absFourier);
xi=peakIndices(1:2:end);
yi=peakIndices(2:2:end);


% determine linear indices
ind=sub2ind(size(absFourier),yi,xi);
[vals, IX]=sort(absFourier(ind), 'desc');

pattNr=rad(xi(IX(nr)));
pattDir=phi(yi(IX(nr)));

% determine score of highest peak
% score=(vals(1)-vals(2))/vals(1);
score=vals(1);

if isPlot==1

    clf;
    
    subplot(1,2,1)
    cla; hold on
    title('SPP Positions')
    scatter(xPos, yPos,'*')
    axis equal
    xlim([-.5 .5])
    ylim([-.5 .5])
   
    subplot(1,2,2)
    cla; hold on
    title('Fourier Transform of SPP Density')
    imagesc(rad,phi,absFourier)
    scatter(rad(xi), phi(yi), 10, 'k', 'filled')
    scatter(pattNr, pattDir, 50, 'k', 'filled')
    
    
    xlabel('Pattern Number r')
    ylabel('\theta')
    colorbar
    colormap jet
    axis tight
end



