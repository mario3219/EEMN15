%% Dynamic focusing
clc,clear
load("PreRF_ImageA.mat");
data = preBeamformed.Signal;
Fs = preBeamformed.SampleFreq;
pitch = preBeamformed.Pitch;
c = preBeamformed.SoundVel;

depths = (0:2048-1).*c/(2*Fs)+3*10^-3;

beamformedImage = zeros(2048,128);

for line = 1:1:128
    line_data = preBeamformed.Signal(:,:,line);
    %data = highpass(data,0.1e6,Fs);
    focused_line = zeros(2048,1);
    for element = 1:1:64
        for sample = 1:1:2048
            depth = depths(sample);
            time_middle = 2*depth/c;
            dx = pitch*abs(32-element);
            d = sqrt(dx^2+depth^2);
            time = 2*d/c;
            delay = time-time_middle;
            sample_delay = round(delay*Fs);
            fixed_sample = sample+sample_delay;
            if fixed_sample > 0 && fixed_sample <= 2048
                focused_line(sample) = focused_line(sample) + line_data(fixed_sample, element);
            end
        end
    end
    beamformedImage(:,line) = focused_line;
end

Image=abs(hilbert(beamformedImage)); %where the fpostbeamformed"-variable is already filtered 
figure; 
imagesc(Image);
colormap(gray)

%%
clc,clear
load("PostRF_Phantom.mat");
data = PostRF.Signal(:,:,1);

Fs = 50e6;
fltrd = highpass(data,0.5e6,Fs);
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
%prebeamformed = highpass(prebeamformed,0.5e6,Fs);

ThreeToTwo=squeeze(sum(prebeamformed,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)

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