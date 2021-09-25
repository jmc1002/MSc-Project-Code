%% FFT Respiratory Impedance Calculation
%Split into 1 sec blocks and take average
%Using longhand power spectrum calculations

flow = flow_subject2_no_d9_d81;
pressure = pressure_subject2_no_d9_d81;
Fs=128;

%% 1.1 Select Input File

%Filename including extension
question1 = 'Enter filename (including extension eg .csv/.txt):';
filename = input(question1,'s'); %Ask user the question and store the response in filename

%Filetype
question2 = 'Enter file type (eg csv/txt):';
filetype = input(question2,'s');

%Sampling Frequency
question3 = 'Enter sampling frequency in Hz:';
Fs = input(question3);


%% 1.2 Read in file

%Read in csv file
if filetype=='csv'
    
    %Read in data from csv file
    data = readtable(filename);

    %Split the data into time, pressure and flow
    %Note:switch pressure and flow for simulated data (f,p) vs real data(p,f)
    time=data(1:1:end,1);
    flow=data(1:1:end,2);
    pressure=data(1:1:end,3);

    %Change from table to double
    time = table2array(time);
    pressure = table2array(pressure);
    flow = table2array(flow);

    %Change orientation of data
    time = time';
    pressure = pressure';
    flow = flow';

end

%Read in txt file
if filetype=='txt'
    
    %Open file for reading
    fileID = fopen(filename,'r');

    %Read in the file contents to data
    %This gives a single column of data with pressure in odd rows and flow in even rows
    formatSpec = '%f %f';
    data = fscanf(fileID,formatSpec);

    %Split the data into pressure and flow using odd/even columns
    pressure=data(1:2:end,:);
    flow=data(2:2:end,:);

    %Change orientation of data
    pressure = pressure';
    flow = flow';
    
end

%% 2.1 Split data into 0.5 sec blocks

%Total number of samples
LP = length(pressure);
LF = length(flow);

%Number of samples per 0.5s
ns_half_sec = Fs/2;

%Pressure
blockP=[];
for i=1:ns_half_sec:LP %from 1 to the last data point in steps of 0.5s
    
    stopBlockP = i+(ns_half_sec-1); %Set stop point for data extract as the start point (i) plus 1 less than the total no points in 0.5s
    extractBlockP = pressure(i:stopBlockP); %Pull out the data from P between the start and stop points
    blockP = [blockP; extractBlockP]; %Store in block
    
end

%Flow
blockF=[];
for i=1:ns_half_sec:LF %from 1 to the last data point in steps of 0.5s
    
    stopBlockF = i+(ns_half_sec-1); %Set stop point for data extract as the start point (i) plus 1 less than the total no points in 0.5s
    extractBlockF = flow(i:stopBlockF); %Pull out the data from F between the start and stop points
    blockF = [blockF; extractBlockF]; %Store in block
    
end

%% 2.2 Combine 0.5 sec blocks to form 1 sec block with a 0.5 sec overlap with the previous
%Combine 2 adjacent blocks starting at each 0.5 sec point
%Should be same number of points per 1 sec block as Fs

blockNoP = size(blockP,1); %Gives number of 0.5 sec blocks in Pressure
blockNoF = size(blockF,1); %Gives number of 0.5 sec blocks in Flow

%Pressure
blockOverlapP = [];
for i=1:1:(blockNoP-1) %from 1 to total no of 0.5 sec blocks in steps of 1 to start at each 0.5 sec point. Stop 1 away from the total number

    %Identify the number of the 2 consecutive blocks from the starting point i
    j=i;
    k=i+1;
    
    %Takes data from specified row (brackets are row:column)
    a=blockP(j,:); 
    b=blockP(k,:);
    
    %Combine the 2 blocks in a single row
    extractOverlapP = [a,b]; 

    %Store in a single location, one new combined block per row
    blockOverlapP = [blockOverlapP; extractOverlapP]; %Store in block
    
end


%Flow
blockOverlapF = [];
for i=1:1:(blockNoF-1) %from 1 to total no of 0.5 sec blocks in steps of 1 to start at each 0.5 sec point. Stop 1 away from the total number

    %Identify the number of the 2 consecutive blocks from the starting point i
    j=i;
    k=i+1;
    
    %Takes data from specified row (brackets are row:column)
    a=blockF(j,:); 
    b=blockF(k,:);
    
    %Combine the 2 blocks in a single row
    extractOverlapF = [a,b]; 

    %Store in a single location, one new combined block per row
    blockOverlapF = [blockOverlapF; extractOverlapF]; %Store in block
    
end

%% 3.1 FFT of Pressure and Flow (1 sec blocks)

window = chebwin(128);
window = window';

%Pressure
fftBlockP = [];
for i=1:1:(blockNoP-1) %For all 1 sec blocks

    %Extract each row in turn
    extractBlockP=blockOverlapP(i,:);

    %Apply window
    extractBlockPwin = extractBlockP.*window;
    
    %Calculate FFT of that row - window
    fftP = fft(extractBlockPwin);
        
    %Calculate FFT of that row - no window
    %fftP = fft(extractBlockP);
    
    %Store in a single location, one new combined block per row
    fftBlockP = [fftBlockP; fftP]; %Store in block
    
end

%Flow
fftBlockF = [];
for i=1:1:(blockNoF-1) %For all 1 sec blocks

    %Extract each row in turn
    extractBlockF=blockOverlapF(i,:);
    
    %Apply window
    extractBlockFwin = extractBlockF.*window;
    
    %Calculate FFT of that row - window
    fftF = fft(extractBlockFwin);
        
    %Calculate FFT of that row - no window
    %fftF = fft(extractBlockF);
        
    %Store in a single location, one new combined block per row
    fftBlockF = [fftBlockF; fftF]; %Store in block
    
end


%% 4.1 Power Spectra of each 1 sec block

%Pressure
GppBlock = [];
for i=1:1:(blockNoP-1) %For all 1 sec blocks

    %Extract each row in turn
    extractfftP=fftBlockP(i,:);

    %Calculate Gpp of that row
    Spi = extractfftP;
    SpConji = conj(extractfftP);
    Gppi = Spi.*SpConji; %.* to perform elementwise multiplication
    
    %Store in a single location, one new combined block per row
    GppBlock = [GppBlock; Gppi]; %Store in block
    
end


%Flow
GffBlock = [];
for i=1:1:(blockNoF-1) %For all 1 sec blocks

    %Extract each row in turn
    extractfftF=fftBlockF(i,:);

    %Calculate Gpp of that row
    Sfi = extractfftF;
    SfConji = conj(extractfftF);
    Gffi = Sfi.*SfConji; %.* to perform elementwise multiplication
    
    %Store in a single location, one new combined block per row
    GffBlock = [GffBlock; Gffi]; %Store in block
    
end

%Cross Spectra
GfpBlock = [];
for i=1:1:(blockNoF-1) %For all 1 sec blocks

    %Extract each row in turn
    extractfftP=fftBlockP(i,:);
    extractfftF=fftBlockF(i,:);

    %Calculate Gpp of that row
    Sfi = extractfftF;
    SpConji = conj(extractfftP);
    Gfpi = Sfi.*SpConji;
    
    %Store in a single location, one new combined block per row
    GfpBlock = [GfpBlock; Gfpi]; %Store in block
    
end

%% 4.2 Average Gpp, Gff and Gfp of each column (not block)
%Effectively average at each freq point

StopAvCol = length(GppBlock);

GppAvBlock = [];
for i=1:1:StopAvCol %For all frequency points

    %Extract each column in turn
    extractGpp=GppBlock(:,i); %Extract row
    
    %Calculate average Gpp of that row
    MGppi = mean(extractGpp); %2 denotes its per row
    
    %Store in a single location, one new combined block per row
    GppAvBlock = [GppAvBlock, MGppi]; %Store in block
    
end

GffAvBlock = [];
for i=1:1:StopAvCol %For all frequency points

    %Extract each column in turn
    extractGff=GffBlock(:,i); %Extract row
    
    %Calculate average Gpp of that column
    MGffi = mean(extractGff);
    
    %Store in a single location, one new combined block per row
    GffAvBlock = [GffAvBlock, MGffi]; %Store in block
    
end

GfpAvBlock = [];
for i=1:1:StopAvCol %For all frequency points

    %Extract each row in turn
    extractGfp=GfpBlock(:,i); %Extract row
    
    %Calculate average Gpp of that row
    MGfpi = mean(extractGfp);
    
    %Store in a single location, one new combined block per row
    GfpAvBlock = [GfpAvBlock, MGfpi]; %Store in block
    
end


%% 5.1 ModZ from 1s block Averages

finishZBlock = length(GfpAvBlock);

ModZBlock = [];
for i=1:1:finishZBlock %For all 1 sec blocks

    %Extract the single values
    AvGfpRi = real(GfpAvBlock(i));
    AvGfpIi = imag(GfpAvBlock(i));
    AvGppi = GppAvBlock(i);

    Denom = (AvGfpRi^2)+(AvGfpIi^2);
    ModZBlockCalc = (AvGppi/(sqrt(Denom)));
    
    %Store in a single location
    ModZBlock = [ModZBlock; ModZBlockCalc]; 
    
end

%% 6.1 Phase Angle from Power and Cross Power Spectrum

finishPhA = size(GfpBlock,2);

PhA = [];
for i=1:1:finishPhA %For all 1 sec blocks

    %Extract the single values
    GfpAvi = GfpAvBlock(i);
    
    %Split into real and imaginary parts
    GfpAvReali = real(GfpAvi);
    GfpAvImagi = imag(GfpAvi);

    %Calculate PhA for that point
    PhACalc = -atand(GfpAvImagi/GfpAvReali);
    
    %Store in a single location
    PhA = [PhA; PhACalc]; 
    
end

%% 7.1 Coherence function of each 1 sec block

cfBlock = [];
for i=1:1:(blockNoP-1) %For all 1 sec blocks

    %Extract the 1 sec block info
    GppBlocki = GppBlock(i,:);
    GffBlocki = GffBlock(i,:);
    
    cfBlocki = mscohere(GppBlocki,GffBlocki);
    cfBlocki=cfBlocki';
    cfBlock = [cfBlock;cfBlocki];
    
end    

%Mean of coherence function (average each column)
StopAvcf = size(cfBlock,2);

cfAv = [];
for i=1:1:StopAvcf %For all frequency points

    %Extract each column in turn
    extractcf=cfBlock(:,i);
    
    %Calculate average Gpp of that row
    Mcf = mean(extractcf);
    
    %Store in a single location, one new combined block per row
    cfAv = [cfAv, Mcf]; %Store in block
    
end
    
%% 8.1 Resistance and Reactance

Res = ModZBlock .* (cosd(PhA));
X = ModZBlock .* (sind(PhA));

%% 9.1 RLC Parameter Estimation

N = 40; %Known that there are 40 frequency points

Ri = [];
Li1 = [];
Li2 = [];
Ei1 = [];
Ei2 = [];
a = [];
b = [];

% Calculate values for each freq
for i=1:N
    
    ReZi = Res(i);
    ImZi = X(i);

    Riint = ReZi; 
    
    Li1int = ImZi/(2*pi*i);
    Li2int = 2*pi*i*ImZi;
    
    Ei1int = ImZi/(2*pi*i);
    Ei2int = 2*pi*i*ImZi;
    
    aint = 1/((2*pi*i)^2);
    
    bint = (2*pi*i)^2;
    
    Ri = [Ri, Riint];
    Li1 = [Li1, Li1int];
    Li2 = [Li2, Li2int];
    Ei1 = [Ei1, Ei1int];
    Ei2 = [Ei2, Ei2int];
    a = [a, aint];
    b = [b, bint];
    
end

RiSum = sum(Ri);
Li1Sum = sum(Li1);
Li2Sum = sum(Li2);
Ei1Sum = sum(Ei1);
Ei2Sum = sum(Ei2);
aSum = sum(a);
bSum = sum(b);

c = (N^2)-(aSum*bSum);

% Calculate RLC 
R = (1/N)*RiSum

L = ((N/c)*Li1Sum) - ((aSum/c)*Li2Sum)

E = ((bSum/c)*Ei1Sum) - ((N/c)*Ei2Sum);

C = 1/E

% Calculate resonant freq based on RLC values
FRes = 1/(2*pi*sqrt(L*C))


%% 10.1 All Plots

tiledlayout(2,3) %Plot using same x axis but different y axis

%Time Calculation
%Fs = 128;
LP = length(pressure);
t = ((0:LP-1)*(1/128)); %Time is each point x the time period (1/Fs)

%Pressure signal
ax1 = nexttile;
plot(t,pressure)
title('Pressure')
xlabel('Time (s)')
ylabel('cm H20')
legend(filename)

%Flow signal
ax2 = nexttile;
plot (t,flow)
title('Flow')
xlabel('Time (s)')
ylabel('Litres/s')
legend(filename)

%%
%Resistance and Reactance
ax3 = nexttile;
f = 0:1:(Fs-1);
plot(f,Res,'x')
hold on
plot(f,X,'x')
xaxis = [0 40]; %Only plot out to 40Hz
xlim(xaxis)
%yaxis = [0 1];
%ylim(yaxis)
ax = gca;
ax.XAxisLocation = 'origin' %Draw the origin
title('Resistance and Reactance')
xlabel('Frequency (Hz)')
ylabel('R and X')
legend('Resistance', 'Reactance')
%%
%ModZ
ax4 = nexttile;
plot(f,ModZBlock,'x')
%Only plot out to 40Hz
xaxis = [1 40];
xlim(xaxis)
yaxis = [0 6];
ylim(yaxis)
title('Modulus of Z')
xlabel('Frequency (Hz)')
ylabel('ZR')
%legend(filename)
%%
%Phase Angle
ax5 = nexttile;
plot(f,PhA,'x')
xaxis = [0 40]; %Only plot out to 40Hz
xlim(xaxis)
%yaxis = [-1 2.5];
%ylim(yaxis)
ax = gca;
ax.XAxisLocation = 'origin' %Display the origin
title('Phase Angle')
xlabel('Frequency (Hz)')
ylabel('PhA')
%legend(filename)
%%
%Coherence Function
ax6 = nexttile;
f1=0:1:Fs;
plot(f1,cfAv,'x')
xaxis = [0 40]; %Only plot out to 40Hz
xlim(xaxis)
yaxis = [0 1]; %Only plot between 0 and 1
ylim(yaxis)
title('Coherence Function')
xlabel('Frequency (Hz)')
ylabel('CF')
%legend(filename)

