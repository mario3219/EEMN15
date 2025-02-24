%% Dynamic focusing
clc,clear

% ladda data
load("PreRF_ImageC.mat");
Fs = preBeamformed.SampleFreq; % 50MHz
pitch = preBeamformed.Pitch; % Distance between elements
c = 1540; %preBeamformed.SoundVel; % Velocity
deadzone = preBeamformed.DeadZone; % Deadzone
deadzone_sample = round((deadzone/c)*Fs);
elementWidth = preBeamformed.ElementWidth;
channels = preBeamformed.Channels;



% Beräkna sampel djup
% samples 0->2048 motsvarar amplituder som funktion av tid, så sampel 100
% motsvarar ett djup, sampel 300 ett annat osv...
depths = (1:2048)*c/(Fs)+deadzone; %meter

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

    for element = 1:1:channels
        %iterera varje sampel, som motsvarar ett viss djup
        for sample = 1:1:2048

            %hämta djupet för en sampel, motsvarar djupet för mitten elementet
            %därav fokus djupet
            depth = depths(sample);

            %tiden det tar för mitten elementet att få en signal
            time_middle = 2.05*depth/c;
            
            %time_middle = 2.07.*depth/c; Ändrad från 2 till 2.05, bättre
            %lateral upplösning. 



            %beräkna avståndet från elementet som ska delayas till mitten
            %elementet, pitch=avstånd från mitt till mitt av element, 0.5:
            %avstånd från element till absolut mitt (jämnt antal element)
            
            dx = pitch*abs(channels/2-element);

           % dx = pitch*abs(channels/2-element+0.5); Utan 0.5 blir lateral
           % resolution bättre

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
beamformedImage = highpass(beamformedImage,4e6,Fs);
Image=abs(hilbert(beamformedImage)); %where the fpostbeamformed"-variable is already filtered 
figure; 
imagesc(Image);
colormap(gray)

%% Inspect individual lines
prebeamformed = new_data;

figure
subplot(511)
plot(prebeamformed(:,1))
subplot(512)
plot(prebeamformed(:,15))
subplot(513)
plot(prebeamformed(:,32))
subplot(514)
plot(prebeamformed(:,45))
subplot(515)
plot(prebeamformed(:,64))

%% PostRF
clc,clear
load("PostRF_Phantom.mat");
data = PostRF.Signal;

Image=abs(hilbert(data)); %where the fpostbeamformed"-variable is already filtered 
figure; 
imagesc(Image);
colormap(gray)

%% Inspect image before beamforming
clc,clear
load("PreRF_ImageC.mat");
Fs = preBeamformed.SampleFreq;

prebeamformed = preBeamformed.Signal;

ThreeToTwo=squeeze(sum(prebeamformed,2));
ThreeToTwo = highpass(ThreeToTwo,4e6,Fs);
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)