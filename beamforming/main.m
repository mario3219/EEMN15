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
Image=abs(hilbert(fltrd)); %where the fpostbeamformed"-variable is already filtered 
figure; 
imagesc(Image);
colormap(gray)

%%
clc,clear
load("PreRF_ImageA.mat");
%%
Fs = 50e6;

prebeamformed = preBeamformed.Signal;
%prebeamformed = highpass(prebeamformed,0.5e6,Fs);

ThreeToTwo=squeeze(sum(prebeamformed,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)

%%
close all

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
clc,clear,close all
load("PreRF_ImageA.mat");
data = preBeamformed.Signal;

goal_height = size(data,1,1);

%Fs = 50e6;
%F = 0.5e6;
%data = highpass(data,0.5e6,Fs);

before = data(:,:,30);
idx1 = find(data(:,15,15) == max(data(:,15,15)));

new_data = [];

% for n = 1:1:128
%     data = preBeamformed.Signal(:,:,:);
%     idx_middle = find(data(:,32,n) == max(data(:,32,n)));
%     idx_middle = idx_middle(1);
%     for t = 1:1:64
%         idx_to_be_delayed = find(data(:,t,n) == max(data(:,t,n)));
%         idx_to_be_delayed = idx_to_be_delayed(1);
%         delay = idx_to_be_delayed - idx_middle;
% 
%         new_column = vertcat(zeros(delay,1), data(:,t,n));
%         new_column = new_column(1:goal_height);
%         data(:,t,n) = new_column;
%     end
% end

for n = 1:1:128
    data = preBeamformed.Signal(:,:,:); % Extracting the data
    idx_middle = find(data(:,32,n) == max(data(:,32,n)), 1, 'first'); % Get first max index
    
    for t = 1:1:64
        idx_to_be_delayed = find(data(:,t,n) == max(data(:,t,n)), 1, 'first'); % Get first max index
        delay = idx_to_be_delayed - idx_middle; % Calculate delay
        
        if delay > 0
            new_column = [zeros(delay,1); data(1:end-delay,t,n)]; % Shift downward
        elseif delay < 0
            new_column = [data(-delay+1:end,t,n); zeros(-delay,1)]; % Shift upward
        else
            new_column = data(:,t,n); % No shift needed
        end
        
        data(:,t,n) = new_column; % Assign modified column
    end
end

after = data(:,:,30);
idx2 = find(data(:,15,15) == max(data(:,15,15)));

%%
ThreeToTwo=squeeze(sum(data,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)
%%
close all

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