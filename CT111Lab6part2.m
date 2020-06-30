%% QPSK Modem

% Modulation
r=1/sqrt(2);
M=4;
L=100;
k=2;
Es=1;
s=[-r,-r;-r,r;r,-r;r,r];
% scatter(s(:,1),s(:,2),'g','filled');
data_packet=randi(2,1,k*L)-1;
transmitted=-ones(L,2);
% Converting into L parallel blocks each having k bits
for i=1:L
    temp=data_packet(i*k-k+1:i*k);
    if(temp==[0 0])
        transmitted(i,:)=s(1,:);
    elseif(temp==[0 1])
        transmitted(i,:)=s(2,:);
    elseif(temp==[1 0])
        transmitted(i,:)=s(3,:);
    elseif(temp==[1 1])
        transmitted(i,:)=s(4,:);
    end
end

% AWGN Channel
noise=-ones(L,2);  % Allocating memory for noise beforehand
SNR=0;
SNR_lin=10.^(0.1*SNR);
sigma2=0.5*Es./(SNR_lin);
sigma=sqrt(sigma2);
noise(:,1)=sigma*randn(1,L);
noise(:,2)=sigma*randn(1,L);
% disp(var(noise(:,1)));
% disp(sigma2);

% Demodulator
received=transmitted+noise;     % Adding AWGN to the transmitted symbols

% % Constellation plot of received symbols superimposed on transmitted
% % symbols
% figure(2);
% scatter(transmitted(:,1),transmitted(:,2),'r','filled');
% hold on;
% scatter(received(:,1),received(:,2),'b','filled');

% Minimum distance symbol detection
demodulated=-ones(L,2);
for i=1:L
    temp=received(i,:);
     min=deuclid(temp,s(1,:));
     for j=1:M
         if(deuclid(temp,s(j,:))<=min)
             min=deuclid(temp,s(j,:));
             demodulated(i,:)=s(j,:);
         end
     end
end

% To find the number of symbols received in error
j=sqrt(-1);
diff=demodulated-transmitted;
comp_diff=diff(:,1)+j*diff(:,2);
comp_diff(comp_diff~=0)=1;
symb_err=(sum(comp_diff))/L;  % prob. of symbols getting received in error
disp(symb_err);
%% Monte Carlo Simulations for QPSK
clc; 
Nsim=100;
r=1/sqrt(2);
M=4;
L=100;
k=2;
Es=1;
s=[-r,-r;-r,r;r,-r;r,r];
symb_err_prob=zeros(1,16);
Psymb_err_UUB=zeros(1,16);
for t=0:1:15
    for z=1:Nsim
        % Modulation
        data_packet=randi(2,1,k*L)-1;
        transmitted=-ones(L,2);
        % Converting into L parallel blocks each having k bits
        for i=1:L
            temp=data_packet(i*k-k+1:i*k);
            if(temp==[0 0])
                transmitted(i,:)=s(1,:);
            elseif(temp==[0 1])
                transmitted(i,:)=s(2,:);
            elseif(temp==[1 0])
                transmitted(i,:)=s(3,:);
            elseif(temp==[1 1])
                transmitted(i,:)=s(4,:);
            end
        end

        % AWGN Channel
        noise=-ones(L,2);  % Allocating memory for noise beforehand
        SNR=t;
        SNR_lin=10.^(0.1*SNR);
        sigma2=0.5*Es./(SNR_lin);
        sigma=sqrt(sigma2);
        noise(:,1)=sigma*randn(1,L);
        noise(:,2)=sigma*randn(1,L);

        % Demodulator
        received=transmitted+noise;     % Adding AWGN to the transmitted symbols

        % Minimum distance symbol detection
        demodulated=-ones(L,2);
        for i=1:L
            temp=received(i,:);
             min=deuclid(temp,s(1,:));
             for j=1:M
                 if(deuclid(temp,s(j,:))<=min)
                     min=deuclid(temp,s(j,:));
                     demodulated(i,:)=s(j,:);
                 end
             end
        end

        % To find the number of symbols received in error
        j=sqrt(-1);
        diff=demodulated-transmitted;
        comp_diff=diff(:,1)+j*diff(:,2);
        comp_diff(comp_diff~=0)=1;
        symb_err=(sum(comp_diff))/L;   % prob. of symbols getting received in error
        symb_err_prob(t+1)=symb_err_prob(t+1)+symb_err;
        disp(symb_err);
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;
    Psymb_err_UUB(t+1)=2*qfunc(sqrt(2)/(2*sigma))+qfunc(2/(2*sigma));
end
SNR=0:1:15;
figure(2);
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_UUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for QPSK Modulation');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid
%% Realistic modem simulator
clearvars;
load('srrcFilter');
p=16;
r=1/sqrt(2);
M=4;
L=100;
k=2;
Es=1;
s=[-r,-r;-r,r;r,-r;r,r];
data_packet=randi(2,1,k*L)-1;
transmitted=zeros(p*L,2);
symbol_stream=zeros(L,2);
% Converting into L parallel blocks each having k bits
for i=1:L
    temp=data_packet(i*k-k+1:i*k);
    if(temp==[0 0])
        transmitted(p*(i-1)+1,:)=s(1,:);
        symbol_stream(i,:)=s(1,:);
    elseif(temp==[0 1])
        transmitted(p*(i-1)+1,:)=s(2,:);
       symbol_stream(i,:)=s(2,:);
    elseif(temp==[1 0])
        transmitted(p*(i-1)+1,:)=s(3,:);
        symbol_stream(i,:)=s(3,:);
    elseif(temp==[1 1])
        transmitted(p*(i-1)+1,:)=s(4,:);
       symbol_stream(i,:)=s(4,:);
    end
end
padded_trasnmitted=[];
padded_transmitted(:,1)=conv(transmitted(:,1),srrcImpulseResponse_alpha03_P16);
padded_transmitted(:,2)=conv(transmitted(:,2),srrcImpulseResponse_alpha03_P16);
padded_transmitted=padded_transmitted(49:end-48,:);
eyediagram(padded_transmitted,2*p);

% Adding noise to the pulse-shaped transmitted sequence
noise=-ones(p*L,2);  % Allocating memory for noise beforehand
SNR=16;
SNR_lin=10.^(0.1*SNR);
sigma2=(0.5*Es./(SNR_lin))*p;
sigma=sqrt(sigma2);
noise(:,1)=sigma*randn(p*L,1);
noise(:,2)=sigma*randn(p*L,1);

% Demodulator
padded_received=padded_transmitted+noise;    % Adding AWGN to the transmitted symbols
received=[];
received(:,1)=conv(padded_received(:,1),srrcImpulseResponse_alpha03_P16);
received(:,2)=conv(padded_received(:,2),srrcImpulseResponse_alpha03_P16);
received=received(49:end-48,:);
eyediagram(received,2*p);


% Minimum distance symbol detection
demodulated=-ones(L,2);
for i=1:L
     temp=received(p*(i-1)+1,:);
     min=deuclid(temp,s(1,:));
     for j=1:M
         if(deuclid(temp,s(j,:))<=min)
             min=deuclid(temp,s(j,:));
             demodulated(i,:)=s(j,:);
         end
     end
end

% To find the number of symbols received in error
j=sqrt(-1);
diff=demodulated-symbol_stream;
comp_diff=diff(:,1)+j*diff(:,2);
comp_diff(comp_diff~=0)=1;
symb_err=(sum(comp_diff))/L;  % prob. of symbols getting received in error
disp(symb_err);
%% Realistic modem QPSK Monte Carlo simulator
clc; 
load('srrcFilter');
Nsim=10000;
p=16;
r=1/sqrt(2);
M=4;
L=100;
k=2;
Es=1;
s=[-r,-r;-r,r;r,-r;r,r];
symb_err_prob=zeros(1,16);
Psymb_err_UUB=zeros(1,16);
for t=0:1:15
    for z=1:Nsim
        % Modulation
        data_packet=randi(2,1,k*L)-1;
        transmitted=zeros(p*L,2);
        symbol_stream=zeros(L,2);
        % Converting into L parallel blocks each having k bits
        for i=1:L
            temp=data_packet(i*k-k+1:i*k);
            if(temp==[0 0])
                transmitted(p*(i-1)+1,:)=s(1,:);
                symbol_stream(i,:)=s(1,:);
            elseif(temp==[0 1])
                transmitted(p*(i-1)+1,:)=s(2,:);
               symbol_stream(i,:)=s(2,:);
            elseif(temp==[1 0])
                transmitted(p*(i-1)+1,:)=s(3,:);
                symbol_stream(i,:)=s(3,:);
            elseif(temp==[1 1])
                transmitted(p*(i-1)+1,:)=s(4,:);
               symbol_stream(i,:)=s(4,:);
            end
        end
        
        % Convolution with the filter
        padded_trasnmitted=[];
        padded_transmitted(:,1)=conv(transmitted(:,1),srrcImpulseResponse_alpha03_P16);
        padded_transmitted(:,2)=conv(transmitted(:,2),srrcImpulseResponse_alpha03_P16);

        % AWGN Channel
        noise=-ones(p*L+96,2);  % Allocating memory for noise beforehand
        SNR=t;
        SNR_lin=10.^(0.1*SNR);
        sigma2=(0.5*Es./(SNR_lin))*p;
        sigma=sqrt(sigma2);
        noise(:,1)=sigma*randn(p*L+96,1);
        noise(:,2)=sigma*randn(p*L+96,1);

        % Demodulator
        padded_received=padded_transmitted+noise;    % Adding AWGN to the transmitted symbols
        received=[];
        received(:,1)=conv(padded_received(:,1),srrcImpulseResponse_alpha03_P16);
        received(:,2)=conv(padded_received(:,2),srrcImpulseResponse_alpha03_P16);
        
        % Minimum distance symbol detection
        demodulated=-ones(L,2);
        for i=1:L
            temp=received(p*(i-1)+97,:);
             min=deuclid(temp,s(1,:));
             for j=1:M
                 if(deuclid(temp,s(j,:))<=min)
                     min=deuclid(temp,s(j,:));
                     demodulated(i,:)=s(j,:);
                 end
             end
        end

        % To find the number of symbols received in error
        j=sqrt(-1);
        diff=demodulated-symbol_stream;
        comp_diff=diff(:,1)+j*diff(:,2);
        comp_diff(comp_diff~=0)=1;
        symb_err=(sum(comp_diff))/L;   % prob. of symbols getting received in error
        symb_err_prob(t+1)=symb_err_prob(t+1)+symb_err;
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;
    Psymb_err_UUB(t+1)=2*qfunc(sqrt(2)/(2*(sigma/4)))+qfunc(2/(2*(sigma/4)));
end
SNR=0:1:15;
figure(2);
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_UUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for Realistic modem QPSK simulator');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid

%% Realistic modem MONTE CARLO simulator for 8 PSK
clc; 
load('srrcFilter');
Nsim=10000;
p=16;
r=1/sqrt(2);
M=8;
L=100;
k=3;
Es=1;
s=[0 1;r r;1 0;r -r;0 -1;-r -r;-1 0;-r r];
symb_err_prob=zeros(1,16);
Psymb_err_UUB=zeros(1,16);
for t=0:1:15
    for z=1:Nsim
        % Modulation
        data_packet=randi(2,1,k*L)-1;
        transmitted=zeros(p*L,2);
        symbol_stream=zeros(L,2);
        % Converting into L parallel blocks each having k bits
        for i=1:L
            temp=data_packet(i*k-k+1:i*k);
             if(temp==[0 0 0])
                transmitted(p*(i-1)+1,:)=s(1,:);
                symbol_stream(i,:)=s(1,:);
            elseif(temp==[0 0 1])
                transmitted(p*(i-1)+1,:)=s(2,:);
                symbol_stream(i,:)=s(2,:);
            elseif(temp==[0 1 1])
                transmitted(p*(i-1)+1,:)=s(3,:);
                symbol_stream(i,:)=s(3,:);
            elseif(temp==[0 1 0])
                transmitted(p*(i-1)+1,:)=s(4,:);
                symbol_stream(i,:)=s(4,:);
            elseif(temp==[1 1 0])
                transmitted(p*(i-1)+1,:)=s(5,:);
                symbol_stream(i,:)=s(5,:);
            elseif(temp==[1 1 1])
                transmitted(p*(i-1)+1,:)=s(6,:);
                symbol_stream(i,:)=s(6,:);
            elseif(temp==[1 0 1])
                transmitted(p*(i-1)+1,:)=s(7,:);
                symbol_stream(i,:)=s(7,:);
            elseif(temp==[1 0 0])
                transmitted(p*(i-1)+1,:)=s(8,:);
                symbol_stream(i,:)=s(8,:);
            end
        end
        
        % Convolution with the filter
        padded_trasnmitted=[];
        padded_transmitted(:,1)=conv(transmitted(:,1),srrcImpulseResponse_alpha03_P16);
        padded_transmitted(:,2)=conv(transmitted(:,2),srrcImpulseResponse_alpha03_P16);

        % AWGN Channel
        noise=-ones(p*L+96,2);  % Allocating memory for noise beforehand
        SNR=t;
        SNR_lin=10.^(0.1*SNR);
        sigma2=(0.5*Es./(SNR_lin))*p;
        sigma=sqrt(sigma2);
        noise(:,1)=sigma*randn(p*L+96,1);
        noise(:,2)=sigma*randn(p*L+96,1);

        % Demodulator
        padded_received=padded_transmitted+noise;    % Adding AWGN to the transmitted symbols
        received=[];
        received(:,1)=conv(padded_received(:,1),srrcImpulseResponse_alpha03_P16)/p;
        received(:,2)=conv(padded_received(:,2),srrcImpulseResponse_alpha03_P16)/p;
        
        % Minimum distance symbol detection
        demodulated=-ones(L,2);
        for i=1:L
            temp=received(p*(i-1)+97,:);
             min=deuclid(temp,s(1,:));
             for j=1:M
                 if(deuclid(temp,s(j,:))<=min)
                     min=deuclid(temp,s(j,:));
                     demodulated(i,:)=s(j,:);
                 end
             end
        end

        % To find the number of symbols received in error
        j=sqrt(-1);
        diff=demodulated-symbol_stream;
        comp_diff=diff(:,1)+j*diff(:,2);
        comp_diff(comp_diff~=0)=1;
        symb_err=(sum(comp_diff))/L;   % prob. of symbols getting received in error
        symb_err_prob(t+1)=symb_err_prob(t+1)+symb_err;
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;   
    Psymb_err_UUB(t+1)=2*qfunc(sqrt((2-sqrt(2)))/(2*(sigma/4)));
end
SNR=0:1:15;
figure(2);
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_UUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for Realistic modem 8PSK simulator');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid

%% Realistic modem MONTE CARLO simulator for 16 APSK
clc; 
load('srrcFilter');
Nsim=10000;
p=16;
r=1/sqrt(2);
M=16;
L=100;
k=4;
Es=15/4;
s=[-1,-1;-1,1;1,-1;1,1;-3,-3;-3,3;3,-3;3,3;
   -3,-1;-3,1;3,-1;3,1;-1,-3;-1,3;1,-3;1,3];
symb_err_prob=zeros(1,16);
Psymb_err_UUB=zeros(1,16);
for t=0:1:15
    for z=1:Nsim
        % Modulation
        data_packet=randi(2,1,k*L)-1;
        transmitted=zeros(p*L,2);
        symbol_stream=zeros(L,2);
        % Converting into L parallel blocks each having k bits
        for i=1:L
            temp=data_packet(i*k-k+1:i*k);
            if(temp==[0 1 0 1])
                transmitted(p*(i-1)+1,:)=s(1,:);
                symbol_stream(i,:)=s(1,:);
            elseif(temp==[1 0 0 1])
                transmitted(p*(i-1)+1,:)=s(2,:);
                symbol_stream(i,:)=s(2,:);
            elseif(temp==[0 1 1 0])
                transmitted(p*(i-1)+1,:)=s(3,:);
                symbol_stream(i,:)=s(3,:);
            elseif(temp==[1 0 1 0])
                transmitted(p*(i-1)+1,:)=s(4,:);
                symbol_stream(i,:)=s(4,:);
            elseif(temp==[0 0 0 0])
                transmitted(p*(i-1)+1,:)=s(5,:);
                symbol_stream(i,:)=s(5,:);
            elseif(temp==[1 1 0 0])
                transmitted(p*(i-1)+1,:)=s(6,:);
                symbol_stream(i,:)=s(6,:);
            elseif(temp==[0 0 1 1])
                transmitted(p*(i-1)+1,:)=s(7,:);
                symbol_stream(i,:)=s(7,:);
            elseif(temp==[1 1 1 1])
                transmitted(p*(i-1)+1,:)=s(8,:);
                symbol_stream(i,:)=s(8,:);
            elseif(temp==[0 1 0 0])
                transmitted(p*(i-1)+1,:)=s(9,:);
                symbol_stream(i,:)=s(9,:);
            elseif(temp==[1 0 0 0])
                transmitted(p*(i-1)+1,:)=s(10,:);
                symbol_stream(i,:)=s(10,:);
            elseif(temp==[0 1 1 1])
                transmitted(p*(i-1)+1,:)=s(11,:);
                symbol_stream(i,:)=s(11,:);
            elseif(temp==[1 0 1 1])
                transmitted(p*(i-1)+1,:)=s(12,:);
                symbol_stream(i,:)=s(12,:);
            elseif(temp==[0 0 0 1])
                transmitted(p*(i-1)+1,:)=s(13,:);
                symbol_stream(i,:)=s(13,:);
            elseif(temp==[1 1 0 1])
                transmitted(p*(i-1)+1,:)=s(14,:);
                symbol_stream(i,:)=s(14,:);
            elseif(temp==[0 0 1 0])
                transmitted(p*(i-1)+1,:)=s(15,:);
                symbol_stream(i,:)=s(15,:);
            elseif(temp==[1 1 1 0])
                transmitted(p*(i-1)+1,:)=s(16,:);
                symbol_stream(i,:)=s(16,:);
            end
        end
        
        % Convolution with the filter
        padded_trasnmitted=[];
        padded_transmitted(:,1)=conv(transmitted(:,1),srrcImpulseResponse_alpha03_P16);
        padded_transmitted(:,2)=conv(transmitted(:,2),srrcImpulseResponse_alpha03_P16);

        % AWGN Channel
        noise=-ones(p*L+96,2);  % Allocating memory for noise beforehand
        SNR=t;
        SNR_lin=10.^(0.1*SNR);
        sigma2=(0.5*Es./(SNR_lin))*p;
        sigma=sqrt(sigma2);
        noise(:,1)=sigma*randn(p*L+96,1);
        noise(:,2)=sigma*randn(p*L+96,1);

        % Demodulator
        padded_received=padded_transmitted+noise;    % Adding AWGN to the transmitted symbols
        received=[];
        received(:,1)=conv(padded_received(:,1),srrcImpulseResponse_alpha03_P16)/p;
        received(:,2)=conv(padded_received(:,2),srrcImpulseResponse_alpha03_P16)/p;
        
        % Minimum distance symbol detection
        demodulated=-ones(L,2);
        for i=1:L
            temp=received(p*(i-1)+97,:);
             min=deuclid(temp,s(1,:));
             for j=1:M
                 if(deuclid(temp,s(j,:))<=min)
                     min=deuclid(temp,s(j,:));
                     demodulated(i,:)=s(j,:);
                 end
             end
        end

        % To find the number of symbols received in error
        j=sqrt(-1);
        diff=demodulated-symbol_stream;
        comp_diff=diff(:,1)+j*diff(:,2);
        comp_diff(comp_diff~=0)=1;
        symb_err=(sum(comp_diff))/L;   % prob. of symbols getting received in error
        symb_err_prob(t+1)=symb_err_prob(t+1)+symb_err;
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;   
    Psymb_err_UUB(t+1)=4*qfunc(2/(2*(sigma/sqrt(p))));
end
SNR=0:1:15;
figure(2);
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_UUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for Realistic modem 16-APSK simulator');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid


%% Realistic modem MONTE CARLO simulator for 32 QAM
clc; 
load('srrcFilter');
Nsim=10000;
p=16;
r=1/sqrt(2);
M=32;
L=100;
k=5;
Es=150/32;
s=[0.5 0.5;0.5 1.5;0.5 2.5;1.5 0.5;1.5 1.5;1.5 2.5;2.5 0.5;2.5 1.5;
   -0.5 0.5;-0.5 1.5;-0.5 2.5;-1.5 0.5;-1.5 1.5;-1.5 2.5;-2.5 0.5;-2.5 1.5;
   -0.5 -0.5;-0.5 -1.5;-0.5 -2.5;-1.5 -0.5;-1.5 -1.5;-1.5 -2.5;-2.5 -0.5;-2.5 -1.5;
   0.5 -0.5;0.5 -1.5;0.5 -2.5;1.5 -0.5;1.5 -1.5;1.5 -2.5;2.5 -0.5;2.5 -1.5;];
symb_err_prob=zeros(1,16);
Psymb_err_UUB=zeros(1,16);
for t=0:1:15
    for z=1:Nsim
        % Modulation
        data_packet=randi(2,1,k*L)-1;
        transmitted=zeros(p*L,2);
        symbol_stream=zeros(L,2);
        % Converting into L parallel blocks each having k bits
        for i=1:L
            temp=data_packet(i*k-k+1:i*k);
            if(temp==[0 0 0 1 1])
                transmitted(p*(i-1)+1,:)=s(1,:);
                symbol_stream(i,:)=s(1,:);
            elseif(temp==[1 0 1 1 1])
                transmitted(p*(i-1)+1,:)=s(2,:);
                symbol_stream(i,:)=s(2,:);
            elseif(temp==[0 1 0 0 1])
                transmitted(p*(i-1)+1,:)=s(3,:);
                symbol_stream(i,:)=s(3,:);
            elseif(temp==[0 0 1 1 0])
                transmitted(p*(i-1)+1,:)=s(4,:);
                symbol_stream(i,:)=s(4,:);
            elseif(temp==[0 1 0 1 0])
                transmitted(p*(i-1)+1,:)=s(5,:);
                symbol_stream(i,:)=s(5,:);
            elseif(temp==[0 1 1 0 0])
                transmitted(p*(i-1)+1,:)=s(6,:);
                symbol_stream(i,:)=s(6,:);
            elseif(temp==[0 1 1 1 1])
                transmitted(p*(i-1)+1,:)=s(7,:);
                symbol_stream(i,:)=s(7,:);
            elseif(temp==[0 0 1 0 1])
                transmitted(p*(i-1)+1,:)=s(8,:);
                symbol_stream(i,:)=s(8,:);
            elseif(temp==[0 0 0 0 0])
                transmitted(p*(i-1)+1,:)=s(9,:);
                symbol_stream(i,:)=s(9,:);
            elseif(temp==[1 0 0 1 0])
                transmitted(p*(i-1)+1,:)=s(10,:);
                symbol_stream(i,:)=s(10,:);
            elseif(temp==[1 1 1 1 0])
                transmitted(p*(i-1)+1,:)=s(11,:);
                symbol_stream(i,:)=s(11,:);
            elseif(temp==[1 0 0 0 1])
                transmitted(p*(i-1)+1,:)=s(12,:);
                symbol_stream(i,:)=s(12,:);
            elseif(temp==[1 1 0 1 1])
                transmitted(p*(i-1)+1,:)=s(13,:);
                symbol_stream(i,:)=s(13,:);
            elseif(temp==[1 1 1 0 1])
                transmitted(p*(i-1)+1,:)=s(14,:);
                symbol_stream(i,:)=s(14,:);
            elseif(temp==[1 1 0 0 0])
                transmitted(p*(i-1)+1,:)=s(15,:);
                symbol_stream(i,:)=s(15,:);
            elseif(temp==[1 0 1 0 0])
                transmitted(p*(i-1)+1,:)=s(16,:);
                symbol_stream(i,:)=s(16,:);
            elseif(temp==[0 1 1 0 1])
                transmitted(p*(i-1)+1,:)=s(17,:);
                symbol_stream(i,:)=s(17,:);
            elseif(temp==[1 1 1 1 1])
                transmitted(p*(i-1)+1,:)=s(18,:);
                symbol_stream(i,:)=s(18,:);
            elseif(temp==[1 0 0 1 1])
                transmitted(p*(i-1)+1,:)=s(19,:);
                symbol_stream(i,:)=s(19,:);
            elseif(temp==[0 1 1 1 0])
                transmitted(p*(i-1)+1,:)=s(20,:);
                symbol_stream(i,:)=s(20,:);
            elseif(temp==[1 0 1 1 0])
                transmitted(p*(i-1)+1,:)=s(21,:);
                symbol_stream(i,:)=s(21,:);
            elseif(temp==[0 0 0 1 0])
                transmitted(p*(i-1)+1,:)=s(22,:);
                symbol_stream(i,:)=s(22,:);
            elseif(temp==[0 0 1 1 1])
                transmitted(p*(i-1)+1,:)=s(23,:);
                symbol_stream(i,:)=s(23,:);
            elseif(temp==[0 1 0 1 1])
                transmitted(p*(i-1)+1,:)=s(24,:);
                symbol_stream(i,:)=s(24,:);
            elseif(temp==[1 1 1 0 0])
                transmitted(p*(i-1)+1,:)=s(25,:);
                symbol_stream(i,:)=s(25,:);
            elseif(temp==[1 1 0 1 0])
                transmitted(p*(i-1)+1,:)=s(26,:);
                symbol_stream(i,:)=s(26,:);
            elseif(temp==[0 0 1 0 0])
                transmitted(p*(i-1)+1,:)=s(27,:);
                symbol_stream(i,:)=s(27,:);
            elseif(temp==[1 1 0 0 1])
                transmitted(p*(i-1)+1,:)=s(28,:);
                symbol_stream(i,:)=s(28,:);
            elseif(temp==[1 0 1 0 1])
                transmitted(p*(i-1)+1,:)=s(29,:);
                symbol_stream(i,:)=s(29,:);
            elseif(temp==[0 0 0 0 1])
                transmitted(p*(i-1)+1,:)=s(30,:);
                symbol_stream(i,:)=s(30,:);
            elseif(temp==[1 0 0 0 0])
                transmitted(p*(i-1)+1,:)=s(31,:);
                symbol_stream(i,:)=s(31,:);
            elseif(temp==[0 1 0 0 0])
                transmitted(p*(i-1)+1,:)=s(32,:);
                symbol_stream(i,:)=s(32,:);
            end
        end
        
             
        
        % Convolution with the filter
        padded_trasnmitted=[];
        padded_transmitted(:,1)=conv(transmitted(:,1),srrcImpulseResponse_alpha03_P16);
        padded_transmitted(:,2)=conv(transmitted(:,2),srrcImpulseResponse_alpha03_P16);

        % AWGN Channel
        noise=-ones(p*L+96,2);  % Allocating memory for noise beforehand
        SNR=t;
        SNR_lin=10.^(0.1*SNR);
        sigma2=(0.5*Es./(SNR_lin))*p;
        sigma=sqrt(sigma2);
        noise(:,1)=sigma*randn(p*L+96,1);
        noise(:,2)=sigma*randn(p*L+96,1);

        % Demodulator
        padded_received=padded_transmitted+noise;    % Adding AWGN to the transmitted symbols
        received=[];
        received(:,1)=conv(padded_received(:,1),srrcImpulseResponse_alpha03_P16)/p;
        received(:,2)=conv(padded_received(:,2),srrcImpulseResponse_alpha03_P16)/p;
        
        % Minimum distance symbol detection
        demodulated=-ones(L,2);
        for i=1:L
            temp=received(p*(i-1)+97,:);
             min=deuclid(temp,s(1,:));
             for j=1:M
                 if(deuclid(temp,s(j,:))<=min)
                     min=deuclid(temp,s(j,:));
                     demodulated(i,:)=s(j,:);
                 end
             end
        end

        % To find the number of symbols received in error
        j=sqrt(-1);
        diff=demodulated-symbol_stream;
        comp_diff=diff(:,1)+j*diff(:,2);
        comp_diff(comp_diff~=0)=1;
        symb_err=(sum(comp_diff))/L;   % prob. of symbols getting received in error
        symb_err_prob(t+1)=symb_err_prob(t+1)+symb_err;
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;   
    Psymb_err_UUB(t+1)=4*qfunc(sqrt(1)/(2*(sigma/sqrt(p))));
end
SNR=0:1:15;
figure(2);
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_UUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for Realistic modem 32-QAM simulator');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid