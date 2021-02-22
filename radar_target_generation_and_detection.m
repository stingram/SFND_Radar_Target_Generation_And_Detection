close all;
clear all;
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e8;
dres = 1;
rmax = 200;
Rtarget = 50;
Vtarget = 70;

 %% FMCW Waveform Generation

B = c/(2*dres);
Tchirp = 5.5*2*rmax/c; % 7e-6;
slope = B/Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq

                                                          
%The number of chirps in one sequence.
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  % for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    %For each time stamp update the Range of the Target for constant velocity. 
    r_t(i) = Rtarget+(Vtarget*t(i));
    td(i) = 2*r_t(i)/c;
    
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+slope*t(i)^2/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-td(i))+slope*(t(i)-td(i))^2/2));
    
    %Now by mixing the Transmit and Receive generate the beat signal
    Mix(i) = Tx(i)*Rx(i);
    
end

%% RANGE MEASUREMENT

%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix=reshape(Mix,[Nr,Nd]);

%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
sig_fft = fft(Mix,Nr);

% Take the absolute value of FFT output
sig_fft = abs(sig_fft./Nr);

% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
sig_fft = sig_fft(1:Nr/2);

%plotting the range
figure ('Name','Range from First FFT')

% plot FFT output 
%plot(r_t,sig_fft)
plot(sig_fft); grid on;
xlabel('Measured Range (meters)')
ylabel('FFT Magnitude')
axis ([0 200 0 0.5]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
xlabel('Measured Range (meters)')
ylabel('Measured Range Rate (meters/second)')
zlabel('2D FFT Power Response')

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

%Select the number of Training Cells in both the dimensions.
Tr = 12;
Td = 12;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 5;
Gd = 3;

% offset the threshold by SNR value in dB
offset = 6;

%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(size(RDM));
RDM_final = zeros(size(RDM));

for m=(Tr+Gr+1):(Nr/2-(Tr+Gr))
    for n=(Td+Td+1):(Nd-(Td+Gd))
        
        % get CUT
        signal = RDM(m,n);
        
        % Get Big Grid
        bg = db2pow(RDM(m-(Tr+Gr):m+(Tr+Gr),n-(Td+Gd):n+(Td+Gd)));
      
        % Get Small Grid
        sg = db2pow(RDM(m-(Gr):m+(Gr),n-(Gd):n+(Gd)));
        
        % Subtract small grid from big grid
        res = sum(sum(bg)) - sum(sum(sg));
        avg = pow2db(res/(numel(bg)-numel(sg)));
        noise_level(m,n) = avg + offset;

        % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
        % CFAR
        if signal > noise_level(m,n)
            RDM_final(m,n)=1;
        end
    end
end

% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 for m=(Tr+Gr+1):(Nr/2-(Tr+Gr))
    for n=(Td+Gd+1):(Nd-(Td+Gd))
        if ((m<(Tr+Gr+1) || m>Nr-(Tr+Gr)) && (n<(Td+Gd+1) || n>Nd-(Td+Gd)))
            RDM_final(m,n) = 0;
        end
    end
 end

%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RDM_final);
colorbar;
xlabel('Measured Range (meters)')
ylabel('Measured Range Rate (meters/second)')
zlabel('Final Range-Doppler-Map Response')


 
 