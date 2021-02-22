# SFND Radar Target Generation and Detection

PUT AN IMAGE HERE

## Implementation Details
---
#### 1. FMCW Waveform Design
Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.

```Matlab
%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 3e8;
dres = 1
rmax = 200
Rtarget = 50;
Vtarget = 70;

 %% FMCW Waveform Generation

B = c/(2*dres);
Tchirp = 5.5*2*rmax/c; % 7e-6;
slope = B/Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq
```

#### 2. Simulation Loop
Simulate Target movement and calculate the beat or mixed signal for every timestamp.

```Matlab
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
```

#### 3. Range FFT (1st FFT)

Implement the Range FFT on the Beat or Mixed Signal and plot the result.

```Matlab
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
subplot(2,1,1)

% plot FFT output 
%plot(r_t,sig_fft)
plot(sig_fft);
axis ([0 200 0 0.5]);
```
![result](https://github.com/godloveliang/SFND-Radar-Target-Generation-and-Detection-/blob/master/img/range_FFT.png)

#### 4. Doppler FFT (2st FFT)

```Matlab
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
```
![result](https://github.com/godloveliang/SFND-Radar-Target-Generation-and-Detection-/blob/master/img/2D_FFT.PNG)

#### 5. 2D CFAR
Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map.

Determine the number of Training cells for each dimension. Similarly, pick the number of guard cells.

```Matlab
%Select the number of Training Cells in both the dimensions.
Tr = 12;
Td = 12;

%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 5;
Gd = 3;

% offset the threshold by SNR value in dB
offset = 6;
```

Create a vector to store the average noise level for each iteration on training cells，and a matrix for the final Range Doppler Map

```Matlab
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);
RDM_final = zeros(size(RDM));
```

Slide the cell under test (CUT) across the grid, leaving a margin for Training and Guard cells around the edges.

```Matlab
for m=(Tr+Gr+1):(Nr/2-(Tr+Gr))
    for n=(Td+Td+1):(Nd-(Td+Gd))
        ...
    end
end
```

For each CUT, convert the signal from logarithmic to linear using the `db2pow` function, next sum the signal level within all the training cells, and get the mean value by dividing by the number of cells in the grid. Finally, convert the value back to db using pow2db.

```Matlab
        % get CUT
        signal = RDM(m,n);
        
        % Get Big Grid
        bg = db2pow(RDM(m-(Tr+Gr):m+(Tr+Gr),n-(Td+Gd):n+(Td+Gd)));
      
        % Get Small Grid
        sg = db2pow(RDM(m-(Gr):m+(Gr),n-(Gd):n+(Gd)));
        
        % Subtract small grid from big grid
        res = sum(sum(bg)) - sum(sum(sg));
        avg = pow2db(res/(numel(bg)-numel(sg)));
```

Add the offset to the noise level to determine the threshold.

```Matlab
        noise_level(m,n) = avg + offset;
```

Next, compare the signal in CUT against this new noise level. If the signal is greater than the threshold, assign it a value of 1.

```Matlab
        % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
        % CFAR
        if signal > noise_level(m,n)
            RDM_final(m,n)=1;
        end
```

The process above will generate a thresholded block, which is smaller than the Range Doppler Map as the CUT cannot be located at the edges of matrix. Hence,few cells will not be thresholded. To keep the map size same set those values to 0. 
```Matlab
 for m=(Tr+Gr+1):(Nr/2-(Tr+Gr))
    for n=(Td+Gd+1):(Nd-(Td+Gd))
        if ((m<(Tr+Gr+1) || m>Nr-(Tr+Gr)) && (n<(Td+Gd+1) || n>Nd-(Td+Gd)))
            RDM_final(m,n) = 0;
        end
    end
 end
```


Number of training cells, number of guard cells and offset were selected by matching the image shared in walkthrough.

![result](https://github.com/godloveliang/SFND-Radar-Target-Generation-and-Detection-/blob/master/img/CA_CFAR.PNG)