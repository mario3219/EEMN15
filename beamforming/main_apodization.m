clc,clear

% Load data
load("PreRF_ImageC.mat");
Fs = preBeamformed.SampleFreq;
pitch = preBeamformed.Pitch;
c = preBeamformed.SoundVel;
deadzone = preBeamformed.DeadZone;
deadzone_sample = round((deadzone/c)*Fs);
elementWidth = preBeamformed.ElementWidth;
channels = preBeamformed.Channels;

% Compute sample depths
depths = (1:2048)*c/(Fs)+deadzone;

% Initialize beamformed image
beamformedImage = zeros(2048,128);

% Debugging matrix
new_data = zeros(2048,64);

% Define apodization window
apodization = hanning(channels);  % Hanning apodization

for line = 1:128
    % Get data for the current line (2048x64)
    line_data = preBeamformed.Signal(:,:,line);

    % Empty vector for the focused line
    focused_line = zeros(2048,1);

    for element = 1:channels
        for sample = 1:2048
            depth = depths(sample);
            time_middle = 2.05*depth/c;
            dx = pitch * abs(channels/2 - element);
            d = sqrt(dx^2 + depth^2);
            time = 2*d/c;
            delay = time - time_middle;
            sample_delay = round(delay * Fs);
            fixed_sample = sample + sample_delay;
            
            if fixed_sample > 0 && fixed_sample <= 2048
                % Apply apodization weight
                weighted_signal = apodization(element) * line_data(fixed_sample, element);
                focused_line(sample) = focused_line(sample) + weighted_signal;

                % Debugging
                new_data(sample, element) = weighted_signal;
            end
        end
    end
    beamformedImage(:,line) = focused_line;
end

% Apply high-pass filtering
beamformedImage = highpass(beamformedImage, 4e6, Fs);
Image = abs(hilbert(beamformedImage));

% Display the beamformed image
figure;
imagesc(Image);
colormap(gray);
