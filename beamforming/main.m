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
clc,clear
load("PreRF_ImageA.mat");

Fs = 50e6;

prebeamformed = preBeamformed.Signal;
%prebeamformed = highpass(prebeamformed,0.5e6,Fs);

ThreeToTwo=squeeze(sum(prebeamformed,2));
Image2=abs(hilbert(ThreeToTwo));
figure; imagesc(Image2); colormap(gray)
