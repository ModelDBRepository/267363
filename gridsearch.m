% Last updated: Apr 27th, 2022
% Grid search of gTRPC4 x gGIRK

clc; clear all; close all; 

x_n_pts = 21; 
y_n_pts = 21; 

gTRPC4 = linspace(0,5,x_n_pts);
gGIRK = linspace(0,5,y_n_pts);

time = 8.5e+3;
dt = 0.01;

[ggirk,gtrpc4] = ndgrid(gGIRK,gTRPC4);
spikes = zeros(y_n_pts,x_n_pts,time/dt);


for x = 1:numel(gTRPC4)
    for y = 1:numel(gGIRK)
        spikes(y,x,:) = mML_TRPC_GIRK(gTRPC4(x),gGIRK(y));
    end
    fprintf([num2str(x/x_n_pts*100),' %%\n'])
end

FileName=[datestr(now, 'yyyymmdd'),'_20x20-gridsearch_steadystates-gCa=p02_time=8p5_range=0-5.mat'];

save(FileName,'-v7.3')