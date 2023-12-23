% Last updated: Apr 27th, 2022
% Modified from Ratte et al. 2018 to include:
% 1) fast h and slow n Na inactivation (from Rho and Prescott 2012)
% 2) GIRK channel (Yim et al. 2015)
% 3) Exp2Syn (from NEURON) used for onset and offset of gGIRK & gTRPC4

function [spike] = mML_TRPC_GIRK(i_gTRPC,i_gGIRK)
graph = 1;
% close all;
%% Constants:
F = 96485; % [C/mol] Faraday's constant
SAV = 3/(10)*1e4; % corrected: surface-to-volume ratio: s/r * scaling factor(um to cm)
% SAV = 3/(10* 0.01); % incorrect: from XPP code for Ratte et al.,2018

n_val = 2; % # of valence e- for Ca2+
C = 2;

%% Initialize variables:
time = 8.5e+3; % [msec]
dt = 0.01;  % time step for forward euler method
loop  = time/dt;   % no. of iterations of euler

% Reversal potentials:
eNa = 50; % [mV] from Ratte et al. 2018
eK =  -90; % [mV] from Ratte et al. 2018
eCa = 100; % [mV] from Ratte et al. 2018
eCAN = 10; % [mV] % increased from 0 (Ratte et al. 2018) for ~ -10 ATPD
eL = -50; % [mV] % increased from -70 (Ratte et al. 2018) for Vrest ~45mV

% Conductance densities:
gL = 2;
gNa = 20; %(Ratte et al. 2018)
gK = 20; % (Ratte et al. 2018)
gfAHP = 10; %20; % [mS/cm2] 50 (Ratte et al.,2018)
gsAHP = 0; % [mS/cm2] 25 (Ratte et al.,2018)
gCa = 0.03; %0.005; % [mS/cm2]

tau_af = 200; % [ms] (Ratte et al.,2018)
tau_as = 2e3; % [ms] (Ratte et al.,2018)
tau_b = 1; % [ms] - (Ratte et al.,2018)
tau_Ca = 2e3; % Ratte et al. 2018 (before = 0.5e3)

beta_m = -1.2; %-15; % -1.2 (Ratte et al., 2018)
gamma_m =  18 ; % 18 (Ratte et al., 2018)

beta_h =  -30; %-27/-30 (Fig. S2 from Rho&Prescott, 2012)
gamma_h = -7; %-8 (Rho&Prescott, 2012)

beta_n = -22; % modified from -30 (Rho&Prescott, 2012)
gamma_n = -5;% 2019-04-08 -5 (Rho&Prescott, 2012)

beta_w = 5; %0 (Ratte et al., 2018)
gamma_w = 5; %10 (Ratte et al., 2018)

% Scaling factors
phi_w = 0.15; % (Rho&Prescott, 2012)
phi_h = 0.03; % Modified from (0.15 from Rho&Prescott) for AP shape
phi_n = 5;


%% TESTING
beta_h =  -27; %-27/-30 (Fig. S2 from Rho&Prescott, 2012)
beta_n = -25; % modified from -30 (Rho&Prescott, 2012)
gamma_n = -5;% 2019-04-08 -5 (Rho&Prescott, 2012)
gamma_h = -8; % (Rho&Prescott, 2012)
phi_n = 10;


% NOTE: tau_rise must be shorter than tau_decay
tau_rise_TRPC =  200;
tau_decay_TRPC =  500;
tau_rise_GIRK =  400;
tau_decay_GIRK = 800;

gCa = 0.02; %0.02;

%% Initialize empty arrays:
t = (1:loop)*dt;
V = zeros(1,loop);

% gating variables
m = zeros(1,loop); % fast Na - activation
h = zeros(1,loop); % fast Na - fast inactivation
n = zeros(1,loop); % fast Na - slow inactivation
w = zeros(1,loop); % slow K
af = zeros(1,loop); % fast AHP
as = zeros(1,loop); % slow AHP
b = zeros(1,loop); % gCa
z = zeros(1,loop); % TRPC
q = zeros(1,loop); % gGIRK

spike = zeros(1,loop);
ref = zeros(1,loop);

Ca_in = zeros(1,loop); % intracellular [Ca]

%% initialize empty vectors
gCAN = zeros(loop,1);
gGIRK = zeros(loop,1);

% currents
INa = zeros(loop,1);
IK = zeros(loop,1);
IfAHP = zeros(loop,1);
IsAHP = zeros(loop,1);
ICat = zeros(loop,1);
ICa = zeros(loop,1);
IGIRK = zeros(loop,1);

%% Initial values
on = 2.5e3; %[msec] agonist application
V(1) = -45;
Ca_in(1) = 0.43; %0.62; % [uM]
gCAN(on/dt:end) = i_gTRPC;
gGIRK(on/dt:end) = i_gGIRK;

% steady state values for gating variables
m(1) = 0;
h(1) = 0.9699;
n(1) = 0.9734;
w(1) = 6.4464e-09;
af(1) = 0.0050;
b(1)=2.1637e-04;
z(1) = 0;
q(1) = 0.5522;

Vth = 0; 
%% Loop through time
for step=1:loop-1
    dV_dt = (0 - gL*(V(step)-eL)...
        - gNa*m_inf(V(step),beta_m,gamma_m)*h(step)*n(step)*(V(step)-eNa)...
        - gK*w(step)*(V(step)-eK)...
        - gfAHP*af(step)*(V(step) - eK)...
        - gsAHP*as(step)*(V(step) - eK)...
        - gCa*b(step)*(V(step) - eCa)...
        - gCAN(step)*z_inf(Ca_in(step))*(V(step) - eCAN)*(exp(-(step*dt-on)/tau_decay_TRPC)-exp(-(step*dt-on)/tau_rise_TRPC))...
        - gGIRK(step)*q(step)*(V(step) - eK)*(exp(-(step*dt-on)/tau_decay_GIRK)-exp(-(step*dt-on)/tau_rise_GIRK))...
        )/C;
    % n (slow inactivation gating variable) was added to prolong Na
    % inactivation
    % Exp2Syn was added to implement EPSP waveforms for changes in GIRK and
    % TRPC4 channel. Exp2Syn is a function in NEURON based on the following
    % equation:
    % G = weight * factor * (exp(-t/tau2) - exp(-t/tau1))
    % tau1 - rise; tau2 - decay; tau2 > tau1
    
    
    
    %% Forward Euler Eq.
    V(step+1) = V(step) + dV_dt*dt;
    
    % fast Na
    m(step) = m_inf(V(step),beta_m,gamma_m);
    % fast inactivation
    dh_dt = phi_h*(h_inf(V(step),beta_h,gamma_h)-h(step))/h_tau(V(step),beta_h,gamma_h);
    h(step+1) = h(step) + dt*dh_dt;
    % slow inactivation
    dn_dt = phi_n*(n_inf(V(step),beta_n,gamma_n)-n(step))/n_tau(V(step),beta_n,gamma_n);
    n(step+1) = n(step) + dt*dn_dt;
    
    % slow K
    dw_dt = phi_w*(w_inf(V(step),beta_w,gamma_w)-w(step))/tau_w(V(step),beta_w,gamma_w);
    w(step+1) = w(step) + dt*dw_dt;
    
    % fast AHP
    daf_dt = (x_inf(V(step)) - af(step)) / tau_af;
    af(step+1) = af(step) + dt*daf_dt;
    
    % slow AHP
    das_dt = (x_inf(V(step)) - as(step)) / tau_as;
    as(step+1) = as(step) + dt*das_dt;
    
    % gCa
    db_dt = (x_inf(V(step)) - b(step)) / tau_b;
    b(step+1) = b(step) + dt*db_dt;
    
    % gCAN (TRPC4)
    z(step) = z_inf(Ca_in(step));
    
    % gGIRK
    dq_dt =  (q_inf(V(step)) - q(step))/q_tau(V(step));
    q(step+1) = q(step) + dt*dq_dt;
    
    % Intracellular Ca
    dCa_dt = -1*SAV * (  gCa*b(step)*(V(step)-eCa)  ) / (n_val*F) - Ca_in(step)/tau_Ca;
    Ca_in(step+1) = Ca_in(step) + dt*dCa_dt;
    
    
    %% Store currents
    INa(step) = gNa*m_inf(V(step),beta_m,gamma_m)*h(step)*n(step)*(V(step)-eNa);
    IK(step) = gK*w(step)*(V(step)-eK);
    IfAHP(step) = gfAHP*af(step)*(V(step) - eK);
    IsAHP(step) = gsAHP*as(step)*(V(step) - eK);
    ICat(step) = gCAN(step)*z_inf(Ca_in(step))*(V(step) - eCAN)*(exp(-(step*dt-on)/tau_decay_TRPC)-exp(-(step*dt-on)/tau_rise_TRPC));
    ICa(step) = gCa*b(step)*(V(step) - eCa);
    IGIRK(step) = gGIRK(step)*q(step)*(V(step)-eK)*(exp(-(step*dt-on)/tau_decay_GIRK)-exp(-(step*dt-on)/tau_rise_GIRK));
    
    spike(step) = (V(step) > Vth).*(~ref(step));
    ref(step+1) = (V(step) > Vth);
end


%% Graph
if graph == 1
    %     close all;
    %     figure('name','V_vs_T')
    %     subplot(3,1,1)
    %     plot(t/1000,V)
    %     ylabel('Voltage [mV]')
    %     % plot(t/1000,h)
    %
    %     subplot(3,1,2)
    %     plot(t/1000,m,'r')%gNa - m
    %     hold on
    %     plot(t/1000,h,'r')%gNa - h
    %     hold on
    %     plot(t/1000,n,'b')%gNa - n
    %     hold on
    %     plot(t/1000,w,'b')%gK
    %     hold on
    %     plot(t/1000,af,'c') % gfahp
    %     hold on
    %     plot(t/1000,b,'m') % gCa
    %     hold on
    %     plot(t/1000,z,'g') % TRPC4
    %     hold on
    %     plot(t/1000,q,'color',[1.0000    0.2706         0]) % GIRK
    %
    %     ylabel('Gating Variables [mS/cm^{2}]')
    
    figure
    plot(t/1000,Ca_in)
    ylabel('[Ca]_{in}')
    xlabel('Time [sec]')
    set(gcf, 'Position',[  288    94   498   841])
    
    %
    % Graph 2 - currents
    %     figure('name','different_currents')
    %     plot(t/1000,INa,'r')
    %     hold on
    %     plot(t/1000,IK,'b')
    %     hold on
    %     plot(t/1000,ICa,'m')
    %     hold on
    %     plot(t/1000,IfAHP,'c')
    %     hold on
    %     plot(t/1000,ICat,'g')
    %     hold on
    %     plot(t/1000,IGIRK,'color',[1.0000    0.2706         0]) % orange
    %     ylabel('Current [\muA/cm^{2}]')
    
    % Graph 3 - final
    figure('name','final')
    % figure(3)
    hold on
    plot(t/1000,V,'r')
    ylabel('Voltage [mV]')
    xlabel('Time [sec]')
    set(gcf,'position',[ 440   378   560   214])
    axis([0.5 8.5 -70 40])
    set(gca,'TickDir','out')
    
    figure('name','IFR')
    %     spike = spike(500/dt:end); % remove first 500 msec
    spike_times = find(spike~=0);
    
    ISI = diff(spike_times);
    ISI = ISI*dt/1000; % convert to sec
    IFR = 1./ISI; % reciprocal of interspike interval
    
    spike_times = spike_times*dt/1000;% convert to sec
    
    x = spike_times(1:end-1);
    y = IFR; %d_IFF(1:end);
    
    scatter(x,y,'.k')
    xlabel('Time (s)'); ylabel('IFR (spk/s)')
    set(gcf,'position',[795   358   560   194])
    set(gca,'TickDir','out')
    xlim([0.5 8.5])
    
    figure('name','sodium gating variables')
    plot(t/1000,m,'r')
    hold on
    plot(t/1000,h,'g') % fast inactivation
    hold on
    plot(t/1000,n,'b') % slow inactivation
    xlabel('Time [sec]')
    legend m h n
    xlim([0.5 8.5])
    
    figure('name','w')
    plot(t/1000,w,'b')
    xlabel('Time [sec]')
    legend w
    
    figure('name','GIRK&TRPC4')
    plot(t/1000,ICat,'r')
    hold on
    plot(t/1000,IGIRK,'b')
    xlabel('Time [sec]')
    
end
end % end of function mML_TRPC_GIRK()

%% fast Na
function [minf] = m_inf(V,beta_m,gamma_m)
% beta_m = -1.2; %-15; % -1.2 (Ratte et al., 2018)
% gamma_m =  18 ;%17;  % 18 (Ratte et al., 2018)
minf = 0.5*(1+tanh((V-beta_m)/gamma_m));
end

%% slow K
function [winf] = w_inf(V, beta_w, gamma_w)
winf = 0.5*(1+tanh((V-beta_w)/gamma_w)); %(Ratte et al., 2018)
end

function [tauw] = tau_w(V, beta_w, gamma_w)
tauw = 1/cosh((V-beta_w)/(2*gamma_w)); %(Ratte et al., 2018)
end

%% fast Na inactivation (modified from Rho&Prescott, 2012)
function [hinf] = h_inf(V, beta_h, gamma_h)
% beta_h =  -27; %-27/-30 (Fig. S2 from Rho&Prescott, 2012)
% gamma_h = -8;% -7; %-8 (Rho&Prescott, 2012)
hinf = 0.5*(1+tanh((V-beta_h)/gamma_h)); %(Fig. S2 from Rho&Prescott, 2012)

end

function [htau] = h_tau(V, beta_h, gamma_h)
% beta_h = -27; %-27/-30 (Rho&Prescott, 2012)
% gamma_h = -8;% -7; %8 (Rho&Prescott, 2012)
htau = 1/cosh((V-beta_h)/(2*gamma_h)); %(Fig. S2 from Rho&Prescott, 2012)
% htau = 2.5;%3.5;
end


%% slow Na inactivation
function [ninf] = n_inf(V,beta_n,gamma_n)
ninf = 1/(1+exp(-(V-beta_n)/gamma_n));
end

function [ntau] = n_tau(V,beta_n,gamma_n)
% ntau = 500/cosh((V-beta_n)/(2*gamma_n));%2000
ntau = 2000; % (Rho&Prescott, 2012)  phi_n = 4
end


%% x = gating variable for both AHP & Ca currents
function[xinf] = x_inf(V)
x_slope = 5; % [mV] % Ratte et al., 2018
x_half = 0; % [mV] % Ratte et al., 2018
xinf = (1+exp(-(V-x_half)/x_slope))^(-1);
end

%% TRPC - nonselective cation channel
function[zinf] = z_inf(Ca)
Ca_half = 0.4; % [uM] Ratte et al., 2018
Ca_slope = 0.2; % [uM] Ratte et al., 2018
zinf = (1+exp(- (Ca - Ca_half)/Ca_slope))^(-1);
end

%% GIRK from Yim et al. 2015
% Changes were made to fit I-V curves,
function[qinf] = q_inf(V)
vhalfl = -70;  %-98.92  ;%(mV)    		: fitted to patch data, Stegen et al. 2012
kl =   20; %10.89 ;   %   (mV)    		: Stegen et al. 2012
% qinf = 1/(1 + exp((V-vhalfl)/kl)); % from Yim et al. 2015
qinf = 1/(1 + exp((V-vhalfl)/kl))+ 0.8/(1 + exp((V-vhalfl)/100)); % modified to fill in (get rid of the hump) in I-V curve
% qinf = 1/(1 + exp((V-vhalfl)/kl))*(1/1.8)+ 0.8/(1 + exp((V-vhalfl)/100))*(0.8/1.8);
% if qinf > 1
%     qinf = 1;
% end
end

function [qtau] = q_tau(V)
vhalft = 67.0828;   % (mV)    		: fitted #100 \muM sens curr 350a,  Stegen et al. 2012
at = 0.00610779;%	 (/ms)   		: Stegen et al. 2012
bt = 0.0817741;%	 (/ms)
qtau = 1/((at*exp(-V/vhalft) + bt*exp(V/vhalft) ))	;
end