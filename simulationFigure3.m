%This Matlab script can be used to generate Figure 3 in the paper:
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
loop=3;

%Number of BS antennas
M=200;

NN=[1 3 10]; %1 and 10
stream10=10:20:200;  %rate1 rate10
stream3=12:15:200;    %rate3

%Radius of cells
radius=500;

% Path loss coefficient
p=3.7;

% correlation at UE and BS
at=0.4;
ar=0;

% Shadow fading variance
stdShadow=5; %dB, standard deviation

% pilot and payload power
PowerPerStream=radius^p*0.5;


% length of a coherence block
U=200;

% prepare to save results
rateDetermUL = zeros(length(NN),length(stream3));
rateDetermDL = zeros(length(NN),length(stream3));
rate_final= zeros(length(NN),length(stream3));


for l=1:loop
    
    disp([num2str(l) ' loop out of ' num2str(loop)]);
    
    %Go through different number of user antennas
    for n=1:length(NN)
        
        disp([num2str(n) ' number of UE antennas out of ' num2str(length(NN))]);
        
        N=NN(n);
        
        % Number of users per cell
        if N == 3
        nbrUES=stream3/N;
        elseif N == 1 || N == 10
        nbrUES=stream10/N;
        end
        
        %Number of cells
        nbrBSs=1;
       
        %Go through different number of users
        for m=1:length(nbrUES)
            
            disp([num2str(m) ' UEs out of ' num2str(length(nbrUES))]);
            
            nbrUEs = nbrUES(m);
            % Number of pilots per cell
            nbrPilots=nbrUEs*N;
            
            % All available pilots
            pilotBook=dftmtx(nbrPilots);
            
            % correlation matrices
            % Rt
            phase = rand(nbrUEs, 1)*2*pi;
            a = at.*exp(1j* phase);
            lambda = zeros(N,nbrUEs);
            for k=1:nbrUEs
                Rt =toeplitz( a(k).^(0:N-1) );
                lambda(:,k) = svd(Rt);  % column vector
            end
            % Shadow Fading
            shadowFadingdB=stdShadow*randn(nbrUEs*nbrBSs,nbrBSs);
            shadowFadingLinear=10.^(shadowFadingdB/10);
            
            % UE to BS Distance
            distanceToBSs=SystemPlot(nbrBSs, nbrUEs, radius);
            pathLoss=distanceToBSs.^(-p).*shadowFadingLinear;
            
            % total power
            pilotPower=PowerPerStream*N*ones(nbrUEs,1);
            payloadPower=pilotPower;
            payloadPowerDL=pilotPower;
            
            % Rr
            phase2 = rand(nbrUEs, 1)*2*pi;
            a2 = ar.*exp(1j* phase2);
            
            Rr = zeros(M*nbrUEs,M);
            Rr_sqrt = zeros(M*nbrUEs,M);
            for k=1:nbrUEs
                Rr(1+(k-1)*M:k*M,:)  =  pathLoss(k) * toeplitz(a2(k).^(0:M-1));  %changed
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
                    ((kron( diag(lambda(:,k)*pilotPower(k)/N), Rr(idx,:)) + 1/nbrPilots*eye(M*N))\ kron(eye(N), Rr(idx,:))) ;
                Phi(:,idx2) = reshape(temp(find(kron(eye(N),ones(M,M)))),M,M*N);
                Phi_tilde(:,idx2) = repmat(Rr(idx,:),1,N) - Phi(:,idx2);
            end
            Z =Phi_tilde* kron( lambda(:).* kron(payloadPower/N,ones(N,1)), speye(M));
            
            Rb = Phi * kron(diag(lambda(:).*kron(payloadPower/N,ones(N,1))),speye(M) );
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
                    Temp =  1/M* Rb_temp* kron(1./(1+d_old),speye(M)) + S;
                    for b=1: nbrPilots-N
                        d_new(b) =1/M * trace(Rb_temp(:,1+(b-1)*M:b*M)/(Temp + rou*speye(M)));
                    end
                end
                delta=d_new;
                T(:,1+(k-1)*M:k*M) = inv(1/M* Rb_temp* kron(1./(1+delta),speye(M)) + S + rou*speye(M));
            end
            
            rate = zeros(nbrPilots,1);
            lambdatemp=lambda(:);
            powertemp = kron(payloadPower/N,ones(N,1));
            
            for b=1:nbrPilots
                idx = floor((b-1)/N)+1;
                rate(b) = log2(1+1/M*trace(Phi(:,1+(b-1)*M:b*M)*T(:, 1+(idx-1)*M : idx*M))*lambdatemp(b)*powertemp(b));
            end
            rateDetermUL(n,m) =  (1-nbrPilots/U)* sum(rate);
            
            %% ============Downlink ==============
            theta=zeros(nbrUEs,1);
            for k=1:nbrUEs
                idx = 1+(k-1)*N*M:k*N*M;
                Temp = Phi(:,idx);
                for i=1:N
                    idx2 = 1+(i-1)*M:i*M;
                    theta(k) = theta(k) + trace(Temp(:,idx2));
                end
            end
            
            % approximative rate of MMSE-SIC receiver
            Gamma = Phi*kron(kron(payloadPowerDL/N./theta,ones(N,1)),speye(M));
            Gamma_temp = zeros(nbrUEs,1);
            for k=1:nbrUEs
                idx = 1+(k-1)*M:k*M;
                idx2 =  1+(k-1)*M*N:k*M*N;
                Gamma_temp(k) = 1/M*trace(Rr(idx,:)*Gamma );
            end
            
            Xi = zeros(nbrPilots,1);
            for b=1:nbrPilots
                idx = 1+(b-1)*M:b*M;
                Xi(b) = 1/M*trace(Phi(:,idx));
            end
            
            rateDL = log2(1+ Xi.^2 .* lambda(:).*kron(payloadPowerDL/N./theta,ones(N,1))./ (1/M*kron(Gamma_temp,ones(N,1)).*lambda(:)+ 1/M^2));
            rateDetermDL(n,m) = (1-nbrPilots/U)* sum(rateDL);
        end
    end
    rate_final=rate_final+(rateDetermUL+rateDetermDL)/2;
end
rate_final = rate_final/loop;


figure; hold on; grid on; box on;
plot(stream10,rate_final(1,1:length(stream10)),'b-',stream10,rate_final(3,1:length(stream10)),'r--',stream3,rate_final(2,:),'k-.','LineWidth',2,'MarkerSize',6)
legend('N=1','N=10','N=3','Location','SouthWest');
axis([0 200 0 250]);
xlabel('Number of total streams NK');
ylabel('Achievable sum SE (bit/s/Hz)');
