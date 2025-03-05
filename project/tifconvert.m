clc; clear;

process("1.97_MHz_1_stack_MMStack_Default.ome.tif");
process("1.97_MHz_2_stack_MMStack_Default.ome.tif");
process("3.84_MHz_1_stack_MMStack_Default.ome.tif");
process("3.84_MHz_2_stack_MMStack_Default.ome.tif");
process("streaming_1_stack_MMStack_Default.ome.tif");
process("streaming_2_stack_MMStack_Default.ome.tif");

function process(file)
    v = VideoWriter(strcat(file, ".mp4"), 'MPEG-4');
    info = imfinfo(file);
    numFrames = numel(info);
    open(v);
    for i = 1:numFrames
        image = imread(file, i);
        
        if ~isa(image, 'uint8')
            image = im2double(image); % Normalize to [0,1]
        end
        
        writeVideo(v, image);
    end
    close(v);
end