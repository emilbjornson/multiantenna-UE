%This Matlab script can be used to generate Figure 1 in the paper:
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
loop=100;

% Number of users per cell
nbrUEs=10;
nbrBSs=1;

%Number of antennas per BS
Mt=[10 50 100 200 300 400];

%Number of antennas per UE
NN=[1 3];

%Radius of cells
radius=500;

% Path loss coefficient
p=3.7;

% Shadow fading variance
stdShadow=5; %dB, standard deviation

% pilot and payload power
snr=0;   % dB
PowerAverage=10^(snr/10);

% length of a coherence block
U=200;

% Rr
phase2 = rand(nbrUEs, 1)*2*pi;
amp2 = 0;
a2 = amp2.*exp(1j* phase2);

% Shadow Fading
shadowFadingdB=stdShadow*randn(nbrUEs*nbrBSs,nbrBSs);
shadowFadingLinear=10.^(shadowFadingdB/10);

%     UE to BS Distance
distanceToBSs=SystemPlot(nbrBSs, nbrUEs, radius);
pathLoss=distanceToBSs.^(-p).*shadowFadingLinear;

% total power
pilotPower = PowerAverage./pathLoss;
payloadPower=pilotPower;
payloadPowerDL=0.5*radius^p*ones(nbrUEs,1);


% prepare to save simulation results
rateDetermUL = zeros(length(NN),length(Mt));
rateInstantUL = zeros(length(NN),length(Mt));
rateInstantUL_MMSE = zeros(length(NN),length(Mt));

rateDetermDL = zeros(length(NN),length(Mt));
rateInstantDL_MMSE = zeros(length(NN),length(Mt));
rateInstantDL = zeros(length(NN),length(Mt));


%Go through number of UE antennas
for n=1:length(NN)
    N=NN(n);
    
    disp([num2str(n) ' UE antennas out of ' num2str(length(NN))]);
    
    % Number of pilots per cell
    nbrPilots=nbrUEs*N;
    
    % All available pilots
    pilotBook=dftmtx(nbrPilots);
    
    % correaltion matrices
    % Rt
    phase = rand(nbrUEs,1)*2*pi;
    amp=0.4;
    a = amp.*exp(1j* phase);
    lambda = zeros(N,nbrUEs);
    for k=1:nbrUEs
        Rt =toeplitz( a(k).^(0:N-1) );
        lambda(:,k) = svd(Rt);  % column vector
    end
    
    %Go through number of BS antennas
    for m=1:length(Mt)
        
        disp([num2str(m) ' BS antennas out of ' num2str(length(Mt))]);
        
        M=Mt(m);
        
        Rr = zeros(M*nbrUEs,M);
        Rr_sqrt = zeros(M*nbrUEs,M);
        for k=1:nbrUEs
            Rr(1+(k-1)*M:k*M,:)  =pathLoss(k)* toeplitz(a2(k).^(0:M-1));
            [u,s,v] = svd(Rr(1+(k-1)*M:k*M,:));
            Rr_sqrt(1+(k-1)*M:k*M,:) = u*sqrt(s)*v';
        end
        
        % deterministic parameters of channel estimation
        Phi = zeros(M, N*M*nbrUEs);
        Phi_tilde = zeros(M, N*M*nbrUEs);
        for k=1: nbrUEs
            idx = 1+(k-1)*M:k*M;
            idx2 = 1+(k-1)*M*N:k*M*N;
            temp = kron( diag(lambda(:,k)*pilotPower(k)/N), Rr(idx,:)) *...
                (kron( diag(lambda(:,k)*pilotPower(k)/N), Rr(idx,:)) + 1/nbrPilots*eye(M*N))^(-1)* kron(eye(N), Rr(idx,:)) ;
            Phi(:,idx2) = reshape(temp(find(kron(eye(N),ones(M,M)))),M,M*N);
            Phi_tilde(:,idx2) = repmat(Rr(idx,:),1,N) - Phi(:,idx2);
        end
        Z =Phi_tilde* kron( lambda(:).* kron(payloadPower/N,ones(N,1)), eye(M));
        
        Rb = Phi * kron(diag(lambda(:).*kron(payloadPower/N,ones(N,1))),eye(M) );
        T=zeros(M,nbrUEs*M);
        
        
        %         ============Uplink========
        %         Large scale approximation, MMSE-SIC receiver
        rou = 1/M;
        S = Z/M;
        for k=1:nbrUEs
            d_old = zeros(nbrPilots-N,1);
            d_new = 1/rou*ones(nbrPilots-N,1);
            idx1 =  (k-1)*M*N;
            idx2 =  k*N*M+1;
            Rb_temp = Rb(:, [1:idx1 idx2:M*N*nbrUEs]);
            while(max(abs(d_old-d_new)) > 1e-5)
                d_old = d_new;
                Temp =  1/M* Rb_temp* kron(1./(1+d_old),eye(M)) + S;
                for b=1: nbrPilots-N
                    d_new(b) =1/M * trace(Rb_temp(:,1+(b-1)*M:b*M)*(Temp + rou*eye(M))^(-1));
                end
            end
            delta=d_new;
            T(:,1+(k-1)*M:k*M) = (1/M* Rb_temp* kron(1./(1+delta),eye(M)) + S + rou*eye(M))^(-1);
        end
        
        rate = zeros(nbrPilots,1);
        lambdatemp=lambda(:);
        powertemp = kron(payloadPower/N,ones(N,1));
        pathLosstemp = kron(pathLoss,ones(N,1));
        for b=1:nbrPilots
            idx = floor((b-1)/N)+1;
            rate(b) = log2(1+1/M*trace(Phi(:,1+(b-1)*M:b*M)*T(:, 1+(idx-1)*M : idx*M))*lambdatemp(b)*powertemp(b));
        end
        rateDetermUL(n,m) =  (1-nbrPilots/U)* sum(rate);
        
        
        %Instantaneous rate
        Rate=0; Rate2=0;
        for l=1:loop
            rate=zeros(nbrUEs,1);
            D = diag(lambda(:).*kron(pilotPower/N,ones(N,1)));
            H = zeros(M,N*nbrUEs);
            for k=1:nbrUEs
                H(:,1+(k-1)*N:k*N) = Rr_sqrt(1+(k-1)*M:k*M,:)*sqrt(0.5)*(randn(M, N)+1j* randn(M, N));
            end
            Noise = sqrt(0.5)*(randn(M, nbrPilots)+1j* randn(M, nbrPilots));
            
            Y = H*sqrt(D)*transp(pilotBook) + Noise;
            
            H_hat = zeros(M,N*nbrUEs);
            for k=1:nbrUEs
                idx = 1+(k-1)*N:k*N;
                idx2 = 1+(k-1)*M:k*M;
                temp = ( kron( sqrt(D(idx, idx)),Rr(idx2,:)))*...
                    (kron(D(idx,idx),Rr(idx2,:)) + 1/nbrPilots*eye(M*N))^(-1)*...
                    reshape(H(:,idx)* sqrt(D(idx,idx))+1/nbrPilots*Noise*conj(pilotBook(:,idx)),M*N,1);
                H_hat(:,idx) = reshape(temp,M,N);
            end
            Omega = diag(lambda(:).*kron(payloadPower/N,ones(N,1)));
            
            % MMSE-SIC receiver
            for k=1:nbrUEs
                idx = (k-1)*N+1:k*N;
                Sigma_k = (H_hat*Omega *H_hat' - H_hat(:,idx)*Omega(idx,idx)*H_hat(:,idx)' + Z + eye(M,M))^(-1);
                rate(k) = real(log2(det(eye(N,N) + Omega(idx,idx)*H_hat(:,idx)'*Sigma_k*H_hat(:,idx))));
            end
            Rate = Rate + sum(rate);
            
            % MMSE receiver
            G =  (H_hat*Omega *H_hat'  + Z + eye(M,M))^(-1)*H_hat *sqrt(Omega);
            all = (G'*H_hat).* conj(G'*H_hat);
            omega = diag(Omega);
            
            noiseP =diag(G'*(Z+eye(M))*G);
            signalP = omega.*diag(all);
            interP = all* omega - signalP;
            
            rate2 = real(log2(1+ signalP./(interP + noiseP)));
            Rate2 = Rate2 + sum(rate2);
        end
        rateInstantUL(n,m) =  (1-nbrPilots/U)*Rate/loop;
        rateInstantUL_MMSE(n,m) = (1-nbrPilots/U)*Rate2/loop;
        
        
        %% ============Downlink ==============
        theta=zeros(nbrUEs,1);
        for k=1:nbrUEs
            idx = 1+(k-1)*N*M:k*N*M;
            Temp = Phi(:,idx);
            for i=1:N
                idx2 = 1+(i-1)*M:i*M;
                theta(k) = theta(k) + 1/M* trace(Temp(:,idx2));
            end
        end
        
        % Instantaneous rate
        H_bar = zeros(N, N*nbrUEs);
        Interf = zeros(N, N*nbrUEs);
        Total2=zeros(N,N);
        
        %Instantaneous rate
        for l=1:loop
            D = diag(lambda(:).*kron(pilotPower/N,ones(N,1)));
            H = zeros(M,N*nbrUEs);
            for k=1:nbrUEs
                H(:,1+(k-1)*N:k*N) = Rr_sqrt(1+(k-1)*M:k*M,:)*sqrt(0.5)*(randn(M, N)+1j* randn(M, N));
            end
            Noise = sqrt(0.5)*(randn(M, nbrPilots)+1j* randn(M, nbrPilots));
            Y = H*sqrt(D)*transp(pilotBook) + Noise;
            H_hat = zeros(M,N*nbrUEs);
            
            for k=1:nbrUEs
                idx = 1+(k-1)*N:k*N;
                idx2 = 1+(k-1)*M:k*M;
                temp = ( kron( sqrt(D(idx, idx)),Rr(idx2,:)))*...
                    (kron(D(idx,idx),Rr(idx2,:)) + 1/nbrPilots*eye(M*N))^(-1)*...
                    reshape(H(:,idx)* sqrt(D(idx,idx))+1/nbrPilots*Noise*conj(pilotBook(:,idx)),M*N,1);
                H_hat(:,idx) = reshape(temp,M,N)/sqrt(theta(k));
            end
            W = H_hat;
            
            % real rate of MMSE-SIC
            Total = W * diag(kron(payloadPowerDL/N, ones(N,1))) * transp(conj(W));
            Temp =  zeros(N, N*nbrUEs);
            Temp2=  zeros(N, N*nbrUEs);
            for k = 1:nbrUEs
                idx = (k-1)*N+1:k*N;
                Temp (:, idx) =diag(sqrt(lambda(:,k))) * H(:, idx)' * W(:,idx) * diag(sqrt(payloadPowerDL(k)/N*ones(N,1)));
                Temp2(:,idx) =diag(sqrt(lambda(:,k))) * H(:, idx)' * Total * H(:,idx) *  diag(sqrt(lambda(:,k))) ;
            end
            H_bar = H_bar + Temp;
            Interf = Interf + Temp2;
        end
        H_bar = H_bar/loop;
        Interf = Interf/loop;
        
        rate = zeros(nbrUEs,1);
        for k=1:nbrUEs
            idx = 1+(k-1)*N:k*N;
            rate(k) = log2(det(eye(N) + H_bar(:,idx)*conj(transp(H_bar(:,idx)))*(Interf(:,idx) - H_bar(:,idx)*conj(transp(H_bar(:,idx)))+eye(N) )^(-1) ));
        end
        rateInstantDL(n,m) = (1-nbrPilots/U)* real(sum(rate));
        
        % approximative rate of MMSE-SIC receiver
        Gamma = Phi*kron(kron(payloadPowerDL/N./theta,ones(N,1)),eye(M));
        Gamma_temp = zeros(nbrUEs,1);
        for k=1:nbrUEs
            idx = 1+(k-1)*M:k*M;
            idx2 =  1+(k-1)*M*N:k*M*N;
            Gamma_temp(k) = 1/M*trace(Rr(idx,:)*(Gamma - Phi(:,idx2)*kron(kron(payloadPowerDL(k)/N/theta(k),ones(N,1)),eye(M)) ));
        end
        
        Xi = zeros(nbrPilots,1);
        for b=1:nbrPilots
            idx = 1+(b-1)*M:b*M;
            Xi(b) = 1/M*trace(Phi(:,idx));
        end
        
        rateDL = log2(1+ Xi.^2.* lambda(:).*kron(payloadPowerDL/N./theta,ones(N,1))./ (1/M*kron(Gamma_temp,ones(N,1)).*lambda(:)+ 1/M));
        rateDetermDL(n,m) = (1-nbrPilots/U)* sum(rateDL);
        
        % real rate of MMSE
        G = zeros(N,nbrUEs*N);
        for k=1:nbrUEs
            idx = (k-1)*N+1:k*N;
            G(:,idx) = (Interf(:,idx) + eye(N))^(-1) * H_bar(:,idx);
        end
        
        signalP = diag(transp(conj(G))*H_bar);
        interfP = zeros(nbrPilots,1);
        for k=1:nbrUEs
            idx = 1+(k-1)*N:k*N;
            interfP(idx) = diag(transp(conj(G(:,idx)))* (Interf(:,idx) + eye(N))*G(:,idx));
        end
        
        rateInstantDL_MMSE(n,m) =  (1-nbrPilots/U)*real(sum(log2(1+ signalP.^2./(interfP - signalP.^2))));
    end
end



figure(1); hold on; box on; grid on;
plot(Mt,rateDetermUL(1,:),'b-',Mt,rateInstantUL(1,:),'bo',Mt,rateDetermDL(1,:),'m-',Mt,rateInstantDL(1,:),'mo',Mt,rateInstantUL_MMSE(1,:),'k--',Mt,rateInstantDL_MMSE(1,:),'k--','LineWidth',2,'MarkerSize',6)
plot(Mt,rateDetermUL(2,:),'b-',Mt,rateInstantUL(2,:),'bo',Mt,rateDetermDL(2,:)+2,'m-',Mt,rateInstantDL(2,:),'mo',Mt,rateInstantUL_MMSE(2,:),'k--',Mt,rateInstantDL_MMSE(2,:),'k--','LineWidth',2,'MarkerSize',6)
ylim([0 200]);
legend('Uplink approximation: MMSE-SIC','Uplink Simulation: MMSE-SIC','Downlink approximation: MMSE-SIC','Downlink simulation: MMSE-SIC','Uplink/Downlink simulation: linear MMSE','Location','NorthWest');
xlabel('Number of BS antennas');
ylabel('Achievable sum SE (bit/s/Hz)');
