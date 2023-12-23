% Last updated: Apr 27th, 2022
% Reproduces panels in Figure 5

clear all; close all; clc; 

load('20200406_20x20_gridsearch_data.mat')

dt = 0.01;

PAUSE = zeros(y_n_pts,x_n_pts); 
TAU = zeros(y_n_pts,x_n_pts);
CLASS = zeros(y_n_pts,x_n_pts); 
class_c = zeros(y_n_pts,x_n_pts); %store colors for each class
color_list = ['r','y','g','c','b'];% r = 1; y = 2; g = 3; c = 4; b = 5 
for i = 1:x_n_pts % gTRPC
    for j = 1:y_n_pts % gGIRK     
        spike = spikes(j,i,:);
        spike = spike(500/dt:end); % remove first 500 msec
        
        spike_times = find(spike~=0);
        
        ISI = diff(spike_times);
        [pause pause_ind] = max(ISI); % maximum interspike interval
        restart = pause_ind+1; % index of first spike after pause << double check
        
        ISI = ISI*dt/1000;
        IFR = 1./ISI; % reciprocal of interspike interval
        d_IFF = IFR./IFR(1); %convert to percent change for tau calculation
        
        spike_times = spike_times*dt/1000;% convert to sec
        PAUSE(j,i) = pause*dt/1000; % convert to sec
        
%         x = spike_times(restart:end-1);
%         y = d_IFF(restart:end);
        
        x = spike_times(1:end-1);
        y = IFR; %d_IFF(1:end);

        figure(1)
        
        scatter(x,y,'.k')
        xlabel('Time (s)'); ylabel('IFR (spk/s)')
        %         f = fit(x,y,'exp1'); % fit exponential function get tau
        %         [A, b] = coeffnames(f);
        %         TAU(j,i)  = -1/b;
        
        set(gcf,'position',[795   358   560   194])
        
        
        %% For spike pattern classification:
        [min_IFF min_IFF_ind] = min(IFR); % find minimum IFF
        pre_pause = IFR(12:min_IFF_ind-1);
        d_pre_pause = diff(pre_pause);
        decel = IFR(min_IFF_ind-1) < 25;
        
        if min_IFF >3 %>3.5 % no long pause (continuous)
            CLASS(j,i) = 1;
        else
            if min_IFF_ind-12 < 2 %15 % I (strong GIRK) - long pause without Na inactivation block
                CLASS(j,i) = 2;
            else
                
                if decel == 0 % acceleration
                    if IFR(min_IFF_ind+1) > IFR(min_IFF_ind+2) %compare adjacent indices
                        CLASS(j,i) = 3; % deceleration
                    else
                        CLASS(j,i) = 4; % acceleration
                    end
                else % deceleration: slow spikes just before pause
%                     if min_IFF > 2 % no long pause (continuous)
%                         CLASS(j,i) = 1;
%                     else
                        CLASS(j,i) = 5; % q + i:% V shaped IFR
%                     end
                end % end of decel ==0
            end
        end
        
        %% FOR DE-BUGGING
        fprintf(['gGIRK = ',num2str(gGIRK(j)),'    gTRPC = ',num2str(gTRPC4(i)),'\n'])
        fprintf(['Class: ',num2str(CLASS(j,i)),'\n'])
        fprintf('\n\n')
        
        class_c(j,i) = color_list(CLASS(j,i)); 
    end
end




figure('name','LongestISI')
% RHEO(RHEO>100) = NaN;
[C,h] = contourf(gTRPC4,gGIRK,PAUSE,50);
set(h,'LineColor','none')
colorbar
xlabel('gTRPC4'); ylabel('gGIRK')
cb = colorbar; ylabel(cb,'Longest ISI (s)')
colormap(flipud(jet))
set(gca,'TickDir','out')
pbaspect([1 1 1])

figure
[ggirk,gtrpc4] = ndgrid(gGIRK,gTRPC4);
% V = surf(gtrpc4, ggirk,CLASS,'FaceColor','texturemap');
V = surf(gtrpc4, ggirk,CLASS);
pbaspect([1 1 1]); view([0 90])
xlabel('gTRPC4'); ylabel('gGIRK')
set(V,'edgecolor','none')
set(gca,'TickDir','out')
% V.EdgeColor = 'flat';

figure
[ggirk,gtrpc4] = ndgrid(gGIRK,gTRPC4);
V = surf(gtrpc4, ggirk,CLASS);
V.EdgeColor = 'flat';
V.Marker = '.';
V.FaceColor ='none';

% V.MarkerFaceColor='flat';
pbaspect([1 1 1]); view([0 90])
xlabel('gTRPC4'); ylabel('gGIRK')
set(gca,'TickDir','out')
% V.EdgeColor = 'flat';

