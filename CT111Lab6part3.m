%% Question 6: Radio frequency simulation with frequency upconversion and downconversion

Rs = 50e3;                              % Symbol rate in symbols/sec (arbitrarily set)
Ts = 1/Rs;                              % Symbol duration in seconds
Lpacket = 100;                          % Packet length in symbols
P_sampPerSym = 16;                      % SRRC oversample rate in samples/symbol
Lsamp = Lpacket*P_sampPerSym;           % Number of samples in the packet
Fsamp = Rs*P_sampPerSym;
Tsamp = 1/Fsamp;
t0 = (0:1696-1)*Tsamp; t0 = t0(:);     % Packet duration in seconds
fc = 3*Rs;                              % Carrier frequency in Hz (150 kHz in this case)
cosineSignal =cos(2*pi*fc*t0);          % Electromagnetic signals
sineSignal =sin(2*pi*fc*t0);            % Electromagnetic signals
Fs=P_sampPerSym*Rs;

% The inphase and the quadrature branches' output after the pulse shaping 
%at the transmitter are denoted as s_tx_inphase and s_tx_quadphase, respectively. 
% The following is the baseband complex-envelope of the transmitted signal.

load('srrcFilter');
p=16;
r=1/sqrt(2);
M=4;
L=100;
k=2;
Es=1;
s=[-r,-r;-r,r;r,-r;r,r];
% symb_err_prob=zeros(1,16);
% Psymb_err_UUB=zeros(1,16);
% for t=0:1:15
%     for z=1:Nsim
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
        s_tx_inphase=padded_transmitted(:,1);
        s_tx_quadphase=padded_transmitted(:,2);
        s_tx = s_tx_inphase + 1i*s_tx_quadphase;
        
        
    
    % Here we upconvert the baseband modulated pulse-shaped signal to radio frequency (RF)
    % cosSignal = cos(2*pi*fc*n/Fs) and sinSignal = sin(2*pi*fc*n/Fs) need
    % to be defined as described in the Lab Manual
    
    % Implement the quadrature notation of the transmitted signal as defined in
    % the lecture slides; note s_tx_RF is not complex-valued anymore
    % Verify mathematically that sqrt(2) multiplier is needed to preserve the symbol
    % energy Es
    
    s_tx_RF = sqrt(2)*(s_tx_inphase.*cosineSignal - s_tx_quadphase.*sineSignal);
    
    % Observe the power spectral density of the transmitted signal
    
    [P_tx,f_tx] = pwelch(s_tx_RF,[],[],[],Fs,'twosided');
    plot(f_tx-Fs/2,10*log10(fftshift(P_tx)),'linewidth',2);
    xlabel('Frequency in Hertz'); ylabel('dB'); 
    title('Power Spectral Density of the Transmitted Signal');
    
    % Add the AWGN at RF. Note that the AWGN is real-valued at the RF
    % noiseSTD has the same value as in the baseband simulator with pulse-shaping
    SNR=12;
    SNR_lin=10.^(0.1*SNR);
    sigma2=(0.5*Es./(SNR_lin))*p;
    noiseSTD=sqrt(sigma2);
    nsig_RF = noiseSTD*randn(1696,1);        % Nsamp=1600;
    r_RF = s_tx_RF + nsig_RF ;
    
    % Observe the power spectral density of the received signal in the presence of the AWGN
    
    P_rx = pwelch(r_RF,[],[],[],Fs,'twosided');
    figure; plot(f_tx-Fs/2,10*log10(fftshift(P_rx)),'linewidth',2);
    xlabel('Frequency in Hertz'); ylabel('dB'); 
    title('Power Spectral Density of the Received Signal in the AWGN');
       
    % Downconvert the received signal at RF to obtain the baseband
    % received signal
    
    r_inphase = r_RF.*cosineSignal;
    r_quadphase = r_RF.*sineSignal;
    padded_received=[];
    padded_received(:,1)=r_inphase;
    padded_received(:,2)=r_quadphase;
     % Demodulator
        received=[];
        received(:,1)=conv(padded_received(:,1),srrcImpulseResponse_alpha03_P16)/p;
        received(:,2)=conv(padded_received(:,2),srrcImpulseResponse_alpha03_P16)/p;
        received(:,2)=-received(:,2);
       % received=received(49:end-48,:);
        
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
        %disp(demodulated');
        %disp(symbol_stream');
        diff=demodulated-symbol_stream;
        comp_diff=diff(:,1)+j*diff(:,2);
        comp_diff(comp_diff~=0)=1;
        symb_err=(sum(comp_diff))/L;   % prob. of symbols getting received in error
        disp(symb_err);

        
 %% Monte Carlo simulations
Nsim=10000;
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
        s_tx_inphase=padded_transmitted(:,1);
        s_tx_quadphase=padded_transmitted(:,2);
        s_tx = s_tx_inphase + 1i*s_tx_quadphase;
        
        
    
    % Here we upconvert the baseband modulated pulse-shaped signal to radio frequency (RF)
    % cosSignal = cos(2*pi*fc*n/Fs) and sinSignal = sin(2*pi*fc*n/Fs) need
    % to be defined as described in the Lab Manual
    
    % Implement the quadrature notation of the transmitted signal as defined in
    % the lecture slides; note s_tx_RF is not complex-valued anymore
    % Verify mathematically that sqrt(2) multiplier is needed to preserve the symbol
    % energy Es
    
    s_tx_RF = sqrt(2)*(s_tx_inphase.*cosineSignal - s_tx_quadphase.*sineSignal);
%     
%     % Observe the power spectral density of the transmitted signal
%     
%     [P_tx,f_tx] = pwelch(s_tx_RF,[],[],[],Fs,'twosided');
%     plot(f_tx-Fs/2,10*log10(fftshift(P_tx)),'linewidth',2);
%     xlabel('Frequency in Hertz'); ylabel('dB'); 
%     title('Power Spectral Density of the Transmitted Signal');
    
    % Add the AWGN at RF. Note that the AWGN is real-valued at the RF
    % noiseSTD has the same value as in the baseband simulator with pulse-shaping
    SNR=t;
    SNR_lin=10.^(0.1*SNR);
    sigma2=(0.5*Es./(SNR_lin))*p;
    sigma=sqrt(sigma2);
    noiseSTD=sigma;
    nsig_RF = noiseSTD*randn(1696,1);        % Nsamp=1600;
    r_RF = s_tx_RF + nsig_RF;
    
%     % Observe the power spectral density of the received signal in the presence of the AWGN
%     
%     P_rx = pwelch(r_RF,[],[],[],Fs,'twosided');
%     figure; plot(f_tx-Fs/2,10*log10(fftshift(P_rx)),'linewidth',2);
%     xlabel('Frequency in Hertz'); ylabel('dB'); 
%     title('Power Spectral Density of the Received Signal in the AWGN');
       
    % Downconvert the received signal at RF to obtain the baseband
    % received signal
    
    r_inphase = r_RF.*cosineSignal;
    r_quadphase = r_RF.*sineSignal;
    padded_received=[r_inphase,r_quadphase];
     % Demodulator
        received=[];
        received(:,1)=(conv(padded_received(:,1),srrcImpulseResponse_alpha03_P16))/p;
        received(:,2)=(conv(padded_received(:,2),srrcImpulseResponse_alpha03_P16))/p;
        %received=received(49:end-48,:);
        received(:,2)=-received(:,2);
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
figure(3);
semilogy(SNR,symb_err_prob,'bo:','linewidth',2,'markerfacecolor','b');
hold on;
semilogy(SNR,Psymb_err_UUB,'r','linewidth',2);
legend('Simulatation','Theory');
title('Symbol Error Probability Evaluation for Radio Frequency simulation(QPSK)');
ylabel('Probability of Symbol Error');
xlabel('E_s/N_0 in dB'); grid

   

    








