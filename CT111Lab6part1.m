% CT111 Lab 6
% Name: Adit Shah
% ID: 201901454

%% QPSK Modem

% Modulation
r=1/sqrt(2);
M=4;
L=100;
k=2;
Es=2.5;
s=[-r,-r;-r,r;r,-r;r,r];
scatter(s(:,1),s(:,2),'r','filled');
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
SNR=13;
SNR_lin=10.^(0.1*SNR);
sigma2=0.5*Es./(SNR_lin);
sigma=sqrt(sigma2);
noise(:,1)=sigma*randn(1,L);
noise(:,2)=sigma*randn(1,L);
% disp(var(noise(:,1)));
% disp(sigma2);

% Demodulator
received=transmitted+noise;     % Adding AWGN to the transmitted symbols

% Constellation plot of received symbols superimposed on transmitted
% symbols
figure(2);
scatter(transmitted(:,1),transmitted(:,2),'r','filled');
hold on;
scatter(received(:,1),received(:,2),'b','filled');
hold on;
xlim([-3 3])
ylim([-3 3])
xL = xlim;
yL = ylim;
line([0 0], yL);  %y-axis
line(xL, [0 0]);  %x-axis
title('Received Symbols at Es/No=13dB');
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

%% Monte Carlo Simulations for QPSK
clc; close all;
Nsim=10000;
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
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;
    Psymb_err_UUB(t+1)=2*qfunc(sqrt(2)/(2*sigma))+qfunc(2/(2*sigma));
end
SNR=0:1:15;
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_UUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for QPSK Modulation');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid


%% 8-PSK Modem
% Modulation
r=1/sqrt(2);
M=8;
L=100;
k=3;
Es=1;
s=[0 1;r r;1 0;r -r;0 -1;-r -r;-1 0;-r r];
scatter(s(:,1),s(:,2),'r','filled');
data_packet=randi(2,1,k*L)-1;
transmitted=-ones(L,2);
% Converting into L parallel blocks each having k bits(gray code)
for i=1:L
    temp=data_packet(i*k-k+1:i*k);
    if(temp==[0 0 0])
        transmitted(i,:)=s(1,:);
    elseif(temp==[0 0 1])
        transmitted(i,:)=s(2,:);
    elseif(temp==[0 1 1])
        transmitted(i,:)=s(3,:);
    elseif(temp==[0 1 0])
        transmitted(i,:)=s(4,:);
    elseif(temp==[1 1 0])
        transmitted(i,:)=s(5,:);
    elseif(temp==[1 1 1])
        transmitted(i,:)=s(6,:);
    elseif(temp==[1 0 1])
        transmitted(i,:)=s(7,:);
    elseif(temp==[1 0 0])
        transmitted(i,:)=s(8,:);
    end
end

% AWGN Channel
noise=-ones(L,2);  % Allocating memory for noise beforehand
SNR=13;
SNR_lin=10.^(0.1*SNR);
sigma2=0.5*Es./(SNR_lin);
sigma=sqrt(sigma2);
noise(:,1)=sigma*randn(1,L);
noise(:,2)=sigma*randn(1,L);

% Demodulator
received=transmitted+noise;     % Adding AWGN to the transmitted symbols

% Constellation plot of received symbols superimposed on transmitted
% symbols
figure(2);
scatter(transmitted(:,1),transmitted(:,2),'r','filled');
hold on;
scatter(received(:,1),received(:,2),'b','filled');
title('Received Symbols at Es/No=13dB');
xlim([-2.5 2.5])
ylim([-2.5 2.5])
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

%% Monte Carlo Simulations for 8-PSK
clc; close all;
Nsim=10000;
r=1/sqrt(2);
M=8;
L=100;
k=3;
Es=1;
s=[0 1;r r;1 0;r -r;0 -1;-r -r;-1 0;-r r];
symb_err_prob=zeros(1,16);
Psymb_err_IUB=zeros(1,16);
for t=0:1:15
    for z=1:Nsim
        % Modulation
        data_packet=randi(2,1,k*L)-1;
        transmitted=-ones(L,2);
        % Converting into L parallel blocks each having k bits
        for i=1:L
            temp=data_packet(i*k-k+1:i*k);
            if(temp==[0 0 0])
                transmitted(i,:)=s(1,:);
            elseif(temp==[0 0 1])
                transmitted(i,:)=s(2,:);
            elseif(temp==[0 1 1])
                transmitted(i,:)=s(3,:);
            elseif(temp==[0 1 0])
                transmitted(i,:)=s(4,:);
            elseif(temp==[1 1 0])
                transmitted(i,:)=s(5,:);
            elseif(temp==[1 1 1])
                transmitted(i,:)=s(6,:);
            elseif(temp==[1 0 1])
                transmitted(i,:)=s(7,:);
            elseif(temp==[1 0 0])
                transmitted(i,:)=s(8,:);
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
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;
    Psymb_err_IUB(t+1)=2*qfunc(sqrt((2-sqrt(2)))/(2*sigma));
end
SNR=0:1:15;
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_IUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for 8-PSK Modulation');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid;

%% 16-APSK Modem

% Modulation
M=16;
r=1/sqrt(2);
x=cosd(22.5);
y=sind(22.5);
L=100;
k=4;
Es=15/4;
s=[-1,-1;-1,1;1,-1;1,1;-3,-3;-3,3;3,-3;3,3;
   -3,-1;-3,1;3,-1;3,1;-1,-3;-1,3;1,-3;1,3];
scatter(s(:,1),s(:,2),'r','filled');
data_packet=randi(2,1,k*L)-1;
transmitted=-ones(L,2);
% Converting into L parallel blocks each having k bits(gray code)
for i=1:L
    temp=data_packet(i*k-k+1:i*k);
    if(temp==[0 1 0 1])
        transmitted(i,:)=s(1,:);
    elseif(temp==[1 0 0 1])
        transmitted(i,:)=s(2,:);
    elseif(temp==[0 1 1 0])
        transmitted(i,:)=s(3,:);
    elseif(temp==[1 0 1 0])
        transmitted(i,:)=s(4,:);
    elseif(temp==[0 0 0 0])
        transmitted(i,:)=s(5,:);
    elseif(temp==[1 1 0 0])
        transmitted(i,:)=s(6,:);
    elseif(temp==[0 0 1 1])
        transmitted(i,:)=s(7,:);
    elseif(temp==[1 1 1 1])
        transmitted(i,:)=s(8,:);
    elseif(temp==[0 1 0 0])
        transmitted(i,:)=s(9,:);
    elseif(temp==[1 0 0 0])
        transmitted(i,:)=s(10,:);
    elseif(temp==[0 1 1 1])
        transmitted(i,:)=s(11,:);
    elseif(temp==[1 0 1 1])
        transmitted(i,:)=s(12,:);
    elseif(temp==[0 0 0 1])
        transmitted(i,:)=s(13,:);
    elseif(temp==[1 1 0 1])
        transmitted(i,:)=s(14,:);
    elseif(temp==[0 0 1 0])
        transmitted(i,:)=s(15,:);
    elseif(temp==[1 1 1 0])
        transmitted(i,:)=s(16,:);
    end
end

% AWGN Channel
noise=-ones(L,2);  % Allocating memory for noise beforehand
SNR=13;
SNR_lin=10.^(0.1*SNR);
sigma2=0.5*Es./(SNR_lin);
sigma=sqrt(sigma2);
noise(:,1)=sigma*randn(1,L);
noise(:,2)=sigma*randn(1,L);

% Demodulator
received=transmitted+noise;     % Adding AWGN to the transmitted symbols

% Constellation plot of received symbols superimposed on transmitted
% symbols
figure(2);
scatter(transmitted(:,1),transmitted(:,2),'r','filled');
hold on;
scatter(received(:,1),received(:,2),'b','filled');
title('Received Symbols at Es/No=13dB');
y = 0;
line([-4,4],[y,y])
y = 2;
line([-4,4],[y,y])
y = -2;
line([-4,4],[y,y])
y = 4;
line([-4,4],[y,y])
y = -4;
line([-4,4],[y,y])
x = 0;
line([x,x],[-4,4])
x = 2;
line([x,x],[-4,4])
x = 4;
line([x,x],[-4,4])
x = -2;
line([x,x],[-4,4])
x = -4;
line([x,x],[-4,4])
xlim([-4 4 ])
ylim([-4 4])
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

%% Monte Carlo Simulations for 16-APSK
clc; close all;
Nsim=10000;
M=16;
L=100;
k=4;
Es=15/4;
s=[-1,-1;-1,1;1,-1;1,1;-3,-3;-3,3;3,-3;3,3;
   -3,-1;-3,1;3,-1;3,1;-1,-3;-1,3;1,-3;1,3];
symb_err_prob=zeros(1,16);
Psymb_err_IUB=zeros(1,16);
for t=0:1:15
    for z=1:Nsim
        % Modulation
        data_packet=randi(2,1,k*L)-1;
        transmitted=-ones(L,2);
        % Converting into L parallel blocks each having k bits
        for i=1:L
            temp=data_packet(i*k-k+1:i*k);
           if(temp==[0 1 0 1])
                transmitted(i,:)=s(1,:);
            elseif(temp==[1 0 0 1])
                transmitted(i,:)=s(2,:);
            elseif(temp==[0 1 1 0])
                transmitted(i,:)=s(3,:);
            elseif(temp==[1 0 1 0])
                transmitted(i,:)=s(4,:);
            elseif(temp==[0 0 0 0])
                transmitted(i,:)=s(5,:);
            elseif(temp==[1 1 0 0])
                transmitted(i,:)=s(6,:);
            elseif(temp==[0 0 1 1])
                transmitted(i,:)=s(7,:);
            elseif(temp==[1 1 1 1])
                transmitted(i,:)=s(8,:);
            elseif(temp==[0 1 0 0])
                transmitted(i,:)=s(9,:);
            elseif(temp==[1 0 0 0])
                transmitted(i,:)=s(10,:);
            elseif(temp==[0 1 1 1])
                transmitted(i,:)=s(11,:);
            elseif(temp==[1 0 1 1])
                transmitted(i,:)=s(12,:);
            elseif(temp==[0 0 0 1])
                transmitted(i,:)=s(13,:);
            elseif(temp==[1 1 0 1])
                transmitted(i,:)=s(14,:);
            elseif(temp==[0 0 1 0])
                transmitted(i,:)=s(15,:);
            elseif(temp==[1 1 1 0])
                transmitted(i,:)=s(16,:);
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
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;
    Psymb_err_IUB(t+1)=4*qfunc(2/(2*sigma));
end
SNR=0:1:15;
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_IUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for 16-APSK Modulation');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid;
%% 32-QAM Modem

% Modulation
M=32;
L=100;
k=5;
Es=150/32;
s=[0.5 0.5;0.5 1.5;0.5 2.5;1.5 0.5;1.5 1.5;1.5 2.5;2.5 0.5;2.5 1.5;
   -0.5 0.5;-0.5 1.5;-0.5 2.5;-1.5 0.5;-1.5 1.5;-1.5 2.5;-2.5 0.5;-2.5 1.5;
   -0.5 -0.5;-0.5 -1.5;-0.5 -2.5;-1.5 -0.5;-1.5 -1.5;-1.5 -2.5;-2.5 -0.5;-2.5 -1.5;
   0.5 -0.5;0.5 -1.5;0.5 -2.5;1.5 -0.5;1.5 -1.5;1.5 -2.5;2.5 -0.5;2.5 -1.5;];
scatter(s(:,1),s(:,2),'r','filled');
data_packet=randi(2,1,k*L)-1;
transmitted=-ones(L,2);
% Converting into L parallel blocks each having k bits(gray code)
for i=1:L
    temp=data_packet(i*k-k+1:i*k);
    if(temp==[0 0 0 1 1])
        transmitted(i,:)=s(1,:);
    elseif(temp==[1 0 1 1 1])
        transmitted(i,:)=s(2,:);
    elseif(temp==[0 1 0 0 1])
        transmitted(i,:)=s(3,:);
    elseif(temp==[0 0 1 1 0])
        transmitted(i,:)=s(4,:);
    elseif(temp==[0 1 0 1 0])
        transmitted(i,:)=s(5,:);
    elseif(temp==[0 1 1 0 0])
        transmitted(i,:)=s(6,:);
    elseif(temp==[0 1 1 1 1])
        transmitted(i,:)=s(7,:);
    elseif(temp==[0 0 1 0 1])
        transmitted(i,:)=s(8,:);
    elseif(temp==[0 0 0 0 0])
        transmitted(i,:)=s(9,:);
    elseif(temp==[1 0 0 1 0])
        transmitted(i,:)=s(10,:);
    elseif(temp==[1 1 1 1 0])
        transmitted(i,:)=s(11,:);
    elseif(temp==[1 0 0 0 1])
        transmitted(i,:)=s(12,:);
    elseif(temp==[1 1 0 1 1])
        transmitted(i,:)=s(13,:);
    elseif(temp==[1 1 1 0 1])
        transmitted(i,:)=s(14,:);
    elseif(temp==[1 1 0 0 0])
        transmitted(i,:)=s(15,:);
    elseif(temp==[1 0 1 0 0])
        transmitted(i,:)=s(16,:);
    elseif(temp==[0 1 1 0 1])
        transmitted(i,:)=s(17,:);
    elseif(temp==[1 1 1 1 1])
        transmitted(i,:)=s(18,:);
    elseif(temp==[1 0 0 1 1])
        transmitted(i,:)=s(19,:);
    elseif(temp==[0 1 1 1 0])
        transmitted(i,:)=s(20,:);
    elseif(temp==[1 0 1 1 0])
        transmitted(i,:)=s(21,:);
    elseif(temp==[0 0 0 1 0])
        transmitted(i,:)=s(22,:);
    elseif(temp==[0 0 1 1 1])
        transmitted(i,:)=s(23,:);
    elseif(temp==[0 1 0 1 1])
        transmitted(i,:)=s(24,:);
    elseif(temp==[1 1 1 0 0])
        transmitted(i,:)=s(25,:);
    elseif(temp==[1 1 0 1 0])
        transmitted(i,:)=s(26,:);
    elseif(temp==[0 0 1 0 0])
        transmitted(i,:)=s(27,:);
    elseif(temp==[1 1 0 0 1])
        transmitted(i,:)=s(28,:);
    elseif(temp==[1 0 1 0 1])
        transmitted(i,:)=s(29,:);
    elseif(temp==[0 0 0 0 1])
        transmitted(i,:)=s(30,:);
    elseif(temp==[1 0 0 0 0])
        transmitted(i,:)=s(31,:);
    elseif(temp==[0 1 0 0 0])
        transmitted(i,:)=s(32,:);
    end
end

% AWGN Channel
noise=-ones(L,2);  % Allocating memory for noise beforehand
SNR=13;
SNR_lin=10.^(0.1*SNR);
sigma2=0.5*Es./(SNR_lin);
sigma=sqrt(sigma2);
noise(:,1)=sigma*randn(1,L);
noise(:,2)=sigma*randn(1,L);

% Demodulator
received=transmitted+noise;     % Adding AWGN to the transmitted symbols

% Constellation plot of received symbols superimposed on transmitted
% symbols
figure(2);
scatter(transmitted(:,1),transmitted(:,2),'r','filled');
hold on;
scatter(received(:,1),received(:,2),'b','filled');
title('Received Symbols at Es/No=13dB');
xlim([-6 6 ])
ylim([-6 6])

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
%% Monte Carlo Simulations for 32 QAM
clc; close all;
Nsim=10000;
M=32;
L=100;
k=5;
Es=150/32;
s=[0.5 0.5;0.5 1.5;0.5 2.5;1.5 0.5;1.5 1.5;1.5 2.5;2.5 0.5;2.5 1.5;
   -0.5 0.5;-0.5 1.5;-0.5 2.5;-1.5 0.5;-1.5 1.5;-1.5 2.5;-2.5 0.5;-2.5 1.5;
   -0.5 -0.5;-0.5 -1.5;-0.5 -2.5;-1.5 -0.5;-1.5 -1.5;-1.5 -2.5;-2.5 -0.5;-2.5 -1.5;
   0.5 -0.5;0.5 -1.5;0.5 -2.5;1.5 -0.5;1.5 -1.5;1.5 -2.5;2.5 -0.5;2.5 -1.5;];
symb_err_prob=zeros(1,16);
Psymb_err_IUB=zeros(1,16);
for t=0:1:15
    for z=1:Nsim
        % Modulation
        data_packet=randi(2,1,k*L)-1;
        transmitted=-ones(L,2);
        % Converting into L parallel blocks each having k bits
        for i=1:L
             temp=data_packet(i*k-k+1:i*k);
             if(temp==[0 0 0 1 1])
                transmitted(i,:)=s(1,:);
            elseif(temp==[1 0 1 1 1])
                transmitted(i,:)=s(2,:);
            elseif(temp==[0 1 0 0 1])
                transmitted(i,:)=s(3,:);
            elseif(temp==[0 0 1 1 0])
                transmitted(i,:)=s(4,:);
            elseif(temp==[0 1 0 1 0])
                transmitted(i,:)=s(5,:);
            elseif(temp==[0 1 1 0 0])
                transmitted(i,:)=s(6,:);
            elseif(temp==[0 1 1 1 1])
                transmitted(i,:)=s(7,:);
            elseif(temp==[0 0 1 0 1])
                transmitted(i,:)=s(8,:);
            elseif(temp==[0 0 0 0 0])
                transmitted(i,:)=s(9,:);
            elseif(temp==[1 0 0 1 0])
                transmitted(i,:)=s(10,:);
            elseif(temp==[1 1 1 1 0])
                transmitted(i,:)=s(11,:);
            elseif(temp==[1 0 0 0 1])
                transmitted(i,:)=s(12,:);
            elseif(temp==[1 1 0 1 1])
                transmitted(i,:)=s(13,:);
            elseif(temp==[1 1 1 0 1])
                transmitted(i,:)=s(14,:);
            elseif(temp==[1 1 0 0 0])
                transmitted(i,:)=s(15,:);
            elseif(temp==[1 0 1 0 0])
                transmitted(i,:)=s(16,:);
            elseif(temp==[0 1 1 0 1])
                transmitted(i,:)=s(17,:);
            elseif(temp==[1 1 1 1 1])
                transmitted(i,:)=s(18,:);
            elseif(temp==[1 0 0 1 1])
                transmitted(i,:)=s(19,:);
            elseif(temp==[0 1 1 1 0])
                transmitted(i,:)=s(20,:);
            elseif(temp==[1 0 1 1 0])
                transmitted(i,:)=s(21,:);
            elseif(temp==[0 0 0 1 0])
                transmitted(i,:)=s(22,:);
            elseif(temp==[0 0 1 1 1])
                transmitted(i,:)=s(23,:);
            elseif(temp==[0 1 0 1 1])
                transmitted(i,:)=s(24,:);
            elseif(temp==[1 1 1 0 0])
                transmitted(i,:)=s(25,:);
            elseif(temp==[1 1 0 1 0])
                transmitted(i,:)=s(26,:);
            elseif(temp==[0 0 1 0 0])
                transmitted(i,:)=s(27,:);
            elseif(temp==[1 1 0 0 1])
                transmitted(i,:)=s(28,:);
            elseif(temp==[1 0 1 0 1])
                transmitted(i,:)=s(29,:);
            elseif(temp==[0 0 0 0 1])
                transmitted(i,:)=s(30,:);
            elseif(temp==[1 0 0 0 0])
                transmitted(i,:)=s(31,:);
            elseif(temp==[0 1 0 0 0])
                transmitted(i,:)=s(32,:);
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
    end
    symb_err_prob(t+1)=symb_err_prob(t+1)/Nsim;
    Psymb_err_IUB(t+1)=4*qfunc(sqrt(1)/(2*sigma));
end
SNR=0:1:15;
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_IUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for 32 QAM Modulation');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid;

