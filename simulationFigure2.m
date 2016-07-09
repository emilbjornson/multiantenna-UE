%This Matlab script can be used to generate Figure 2 in the paper:
%
%Xueru Li, Emil Björnson, Shidong Zhou, Jing Wang, “Massive MIMO with
%Multi-Antenna Users: When are Additional User Antennas Beneficial?,”
%Proceedings of International Conference on Telecommunications (ICT),
%Thessaloniki, Greece, May 2016.
%
%Download article: http://arxiv.org/pdf/1603.09052
%
%This is version 1.0 (Last edited: 2016-07-09)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original paper listed above.
%
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.


clear all
close all
clc


%% ====================System Basic Parameters==================

%Number of Monte Carlo realizations
loop=200;

% Number of users per cell
nbrUEs=10;
nbrBSs=1;

%Number of antennas per BS
Mt = 10.^(1:7);

%Number of antennas per UE
N=3;

% Number of pilots per cell
nbrPilots=nbrUEs*N;

%Radius of cells
radius=500;

% Path loss coefficient
p=3.7;

% Shadow fading variance
stdShadow=5; %dB, standard deviation

% length of a coherence block
U=200;

snr=0;   % dB
PowerAverage=10^(snr/10);

% All available pilots
pilotBook=dftmtx(nbrPilots);

rateDetermUL = zeros(1,length(Mt));
rateDetermDL = zeros(1,length(Mt));
rateLimitUL = 0;
rateLimitDL = 0;

% range of power reduction coefficient (this is called alpha in the paper)
tauRange = [0.5 1];

for t = 1:length(tauRange)
    
    tau=tauRange(t);
    
    for l=1:loop
        % Shadow Fading
        shadowFadingdB=stdShadow*randn(nbrUEs*nbrBSs,nbrBSs);
        shadowFadingLinear=10.^(shadowFadingdB/10);
        
        % UE to BS Distance
        distanceToBSs=SystemPlot(nbrBSs, nbrUEs, radius);
        pathLoss=distanceToBSs.^(-p).*shadowFadingLinear;
        
        % Total pilot and payload power
        powerTotal = PowerAverage./pathLoss;
        powerDL = radius^p*0.5*ones(nbrUEs,1);
        
        % correlation matrices
        phase = rand(nbrUEs, 1)*2*pi;
        amp = 0.4;
        a = amp.*exp(1j* phase);
        
        lambda = zeros(N,nbrUEs);
        for k=1:nbrUEs
            R = pathLoss(k)*toeplitz( a(k).^(0:N-1) );
            lambda(:,k) = svd(R);  % column vector
        end
        
        for m=1:length(Mt)
            
            M=Mt(m);
            
            % power reduction
            pilotPower = powerTotal/(M^tau);
            payloadPower = powerTotal/(M^(1-tau));
            payloadPowerDL = powerDL/(M^(1-tau));
            
            % deterministic parameters of channel estimation
            alpha = zeros(N, nbrUEs);
            Beta=zeros(nbrUEs,1);
            for k=1: nbrUEs
                alpha(:,k) = lambda(:,k).*(pilotPower(k)/N) ./( lambda(:,k).*(pilotPower(k)/N) + 1/nbrPilots);
                Beta (k) = sum(lambda(:,k)* payloadPower(k)/N .*(1-alpha(:,k)));
            end
            beta = sum(Beta) + 1;
            
            Rb=lambda(:).*alpha(:).* kron(payloadPower,ones(N,1))*M/N;
            
            %==============Large scale approximation: UL==============
            rou = beta;
            d_old = zeros(nbrPilots,1);
            d_new = 1/rou*ones(nbrPilots,1);
            while(max(abs(d_old-d_new)) > 1e-5)
                d_old = d_new;
                Z = 1/M*sum(Rb./(1+d_old)) ;
                d_new = Rb/(Z+rou);
            end
            delta=d_new;
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            T=1/(1/M*sum(Rb./(1+delta)) +rou);
            rateDetermUL(m)  = rateDetermUL(m) + (1-nbrPilots/U)* sum( log2(1+T*Rb) );
            
            % =====================large scale approximation: DL=================
            theta = sum(alpha,1);
            gamma =  sum(alpha(:) .*  kron(payloadPower/N./theta',ones(N,1) ))/M;
            rateDetermDL (m) =  rateDetermDL (m) + (1-nbrPilots/U)*sum( log2( 1+ lambda(:).*alpha(:).* ...
                kron(payloadPower/N,ones(N,1)).*alpha(:)./kron(theta',ones(N,1)) ./(gamma*lambda(:) + 1/M)));
            
        end
        
        % ===============Limit: UL=========================
        alpha_bar = lambda(:).*kron(powerTotal,ones(N,1))/N*nbrPilots;
        if (tau~=1)
            beta_bar = 1;
        else
            beta_bar = 1+ sum(lambda(:).*kron(powerTotal,ones(N,1))/N);
        end
        Rb_bar = alpha_bar.*lambda(:).*kron(powerTotal,ones(N,1))/N;
        T_bar = 1/beta_bar;
        
        rateLimitUL = rateLimitUL + (1-nbrPilots/U)* sum( log2(1+T_bar*Rb_bar) );
        
        % ===============Limit:DL=========================
        theta_limit = kron(sum(reshape(alpha_bar,N,nbrUEs),1)',ones(N,1));
        if(tau==1)
            gamma = sum(kron(powerTotal/N,ones(N,1)).*alpha_bar./ theta_limit);
            
        else
            gamma = 0;
        end
        rateLimitDL = rateLimitDL + (1-nbrPilots/U)*sum(log2( 1+ lambda(:).* alpha_bar.^2.* ...
            kron(powerTotal/N,ones(N,1))./theta_limit./(1+gamma*lambda(:))));
        
    end
    rateDetermUL = rateDetermUL/loop;
    rateDetermDL = rateDetermDL/loop;
    rateLimitUL= rateLimitUL/loop;
    rateLimitDL=rateLimitDL/loop;
    
    %Plot result for one value of tau
    figure(1)
    semilogx( Mt, ones(1,length(Mt))*rateLimitUL,'b--',Mt, rateDetermUL,'bo-',Mt, ones(1,length(Mt))*rateLimitDL,'r--',Mt,rateDetermDL,'ro-','LineWidth',2,'MarkerSize',6)
    grid on
    hold on
    
end

ylim([0 80]);
legend('Uplink SE limit','Uplink SE','Downlink SE limit','Downlink SE','Location','NorthWest');
xlabel('Number of BS antennas');
ylabel('Achievable sum SE (bit/s/Hz)');
