%% Dynamic focusing
clc,clear

% ladda data
load("PreRF_ImageA.mat");
data = preBeamformed.Signal;
Fs = preBeamformed.SampleFreq;
pitch = preBeamformed.Pitch;
c = preBeamformed.SoundVel;
deadzone = preBeamformed.DeadZone;
deadzone_sample = round((deadzone/c)*Fs);
elementWidth = preBeamformed.ElementWidth;

% Beräkna sampel djup
% samples 0->2048 motsvarar amplituder som funktion av tid, så sampel 100
% motsvarar ett djup, sampel 300 ett annat osv...
depths = (1:2048)*c/(Fs)+deadzone;

% tom beamformed image, där varje focused linje ska lagras
beamformedImage = zeros(2048,128);

%matris för debugging
new_data = zeros(2048,64);

for line = 1:1:128
    %hämta data för en linje, 2048x64
    line_data = preBeamformed.Signal(:,:,line);

    %HP filter, bara kommentera ut

    %tom vektor för en fokuserad linje
    focused_line = zeros(2048,1);

    for element = 1:1:64
        %iterera varje sampel, som motsvarar ett viss djup
        for sample = deadzone_sample:1:2048

            %hämta djupet för en sampel, motsvarar djupet för mitten elementet
            %därav fokus djupet
            depth = depths(sample);

            %tiden det tar för mitten elementet att få en signal
            time_middle = 2*depth/c;

            %beräkna avståndet från elementet som ska delayas till mitten
            %elementet
            dx = pitch*abs(32-element-1);%+elementWidth*abs(32-element);

            %beräkna avståndet från elementet som ska delayas till fokus
            %punkten
            d = sqrt(dx^2+depth^2);

            %beräkna tiden för ekot från elementet som ska delayas till
            %fokus punkten
            time = 2*d/c;

            %delay som ska sättas blir skillnaden mellan tiden för mitten
            %elementet och tiden för elementet som ska delayas
            delay = time-time_middle;

            %överför till sample delay
            sample_delay = round(delay*Fs);

            %hitta sampeln som är den egentliga som ska användas
            fixed_sample = sample+sample_delay;
            
            %if satsen kan göras mer noggrannt här, hitills finns några
            %fel:
            % * Om delayen som behövs går out-of-array-bounds, så kommer
            %   den inte att delayas
            %   Potentiell lösning: välj sampeln på index 2048 om OOB.
            if fixed_sample > 0 && fixed_sample <= 2048
                %sparar den riktiga sampeln som ska användas, som hittats i
                %hänsyn till delay, och summerar med alla andra lagrade
                %samples
                focused_line(sample) = focused_line(sample) + line_data(fixed_sample, element);
                
                %lagrar vektorn för debugging
                new_data(sample,element) = line_data(fixed_sample, element);
            end

        end
    end
    %efter att ha summerat alla 64 element så har en linje formats, lägg
    %till i indexerad plats i beamform image
    beamformedImage(:,line) = focused_line;
end
beamformedImage = highpass(beamformedImage,0.5e5,Fs);
Image=abs(hilbert(beamformedImage)); %where the fpostbeamformed"-variable is already filtered 
figure; 
imagesc(Image);
colormap(gray)

%%
prebeamformed = new_data;

figure
subplot(511)
plot(prebeamformed(:,1)), hold on, xlim([810 900])
subplot(512)
plot(prebeamformed(:,31)), hold on, xlim([810 900])
subplot(513)
plot(prebeamformed(:,32)), hold on, xlim([810 900])
subplot(514)
plot(prebeamformed(:,45)), hold on, xlim([810 900])
subplot(515)
plot(prebeamformed(:,64)), hold on, xlim([810 900])
%%
prebeamformed = preBeamformed.Signal;
%prebeamformed = highpass(prebeamformed,0.5e6,Fs);

ThreeToTwo=squeeze(sum(prebeamformed,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)

%%
clc,clear
load("PostRF_Phantom.mat");
data = PostRF.Signal;

Fs = 50e6;
fltrd = highpass(data,0.5e3,Fs);
Image=abs(hilbert(fltrd)); %where the fpostbeamformed"-variable is already filtered 
figure; 
imagesc(Image);
colormap(gray)

%%
Fs = 50e6;
ThreeToTwo=squeeze(sum(data,2));
Image=abs(hilbert(data)); %where the fpostbeamformed"-variable is already filtered 
figure; 
imagesc(Image);
colormap(gray)

%%
clc,clear
load("PreRF_ImageA.mat");
Fs = 50e6;

prebeamformed = preBeamformed.Signal;
prebeamformed = highpass(prebeamformed,0.5e7,Fs);

ThreeToTwo=squeeze(sum(prebeamformed,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray), title("Innan")

%%

prebeamformed = preBeamformed.Signal;

figure
subplot(511)
plot(prebeamformed(:,1,1)), hold on, xlim([810 900])
subplot(512)
plot(prebeamformed(:,15,1)), hold on, xlim([810 900])
subplot(513)
plot(prebeamformed(:,30,1)), hold on, xlim([810 900])
subplot(514)
plot(prebeamformed(:,45,1)), hold on, xlim([810 900])
subplot(515)
plot(prebeamformed(:,60,1)), hold on, xlim([810 900])

prebeamformed = preBeamformed.Signal;
%prebeamformed = highpass(prebeamformed,0.5e6,Fs);

ThreeToTwo=squeeze(sum(prebeamformed,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)

%%
clc,clear
load("PreRF_ImageA.mat");
data = preBeamformed.Signal;

goal_height = size(data,1,1);
Fs = 50e6;

data = preBeamformed.Signal(:,:,:);
for n = 1:1:128
    idx_middle = find(data(:,32,n) == max(data(:,32,n)));
    idx_middle = idx_middle(1);
    for t = 1:1:64
        idx_to_be_delayed = find(data(:,t,n) == max(data(:,t,n)));
        idx_to_be_delayed = idx_to_be_delayed(1);
        delay = idx_to_be_delayed - idx_middle;

        new_column = vertcat(zeros(delay,1), data(:,t,n));
        new_column = new_column(1:goal_height);
        data(:,t,n) = new_column;
    end
end
ThreeToTwo=squeeze(sum(data,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)
%%
data = preBeamformed.Signal;
ThreeToTwo=squeeze(sum(data,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)
%%
data = preBeamformed.Signal;
figure
subplot(511)
plot(data(:,1,1)), hold on, xlim([810 900])
subplot(512)
plot(data(:,15,1)), hold on, xlim([810 900])
subplot(513)
plot(data(:,30,1)), hold on, xlim([810 900])
subplot(514)
plot(data(:,45,1)), hold on, xlim([810 900])
subplot(515)
plot(data(:,60,1)), hold on, xlim([810 900])

%%
data = preBeamformed.Signal;
figure
subplot(511)
plot(data(:,1,1))
subplot(512)
plot(data(:,15,1))
subplot(513)
plot(data(:,30,1))
subplot(514)
plot(data(:,45,1))
subplot(515)
plot(data(:,60,1))