<<<<<<< HEAD
function [XY] = ImageBlank(fileprefix, channelBF,channelN, channelR)

framerangemin = 1;
framerangemax = 191;
XY = num2str(2, '%02d');

for k = framerangemin:framerangemax
    framenum = num2str(k,'%03d');
    numcentroids = 0;    
    frame_n = imread(strcat(fileprefix,'c',num2str(channelN),'t',framenum,'.tif'));
    blank_n = imread(strcat('blankt1','xy',XY,'c',num2str(channelN),'.tif'));
    difference_n = min(min(frame_n)) - min(min(blank_n));
    
       
    frame_bf = imread(strcat(fileprefix,'c',num2str(channelBF),'t',framenum,'.tif'));
    blank_bf = imread(strcat('blankt1','xy',XY,'c',num2str(channelBF),'.tif'));
    difference_bf = min(min(frame_bf)) - min(min(blank_bf));
    
    frame_r = imread(strcat(fileprefix,'c',num2str(channelR),'t',framenum,'.tif'));
    blank_r = imread(strcat('blankt1','xy',XY,'c',num2str(channelR),'.tif'));
    difference_r = min(min(frame_r)) - min(min(blank_r));

    if k == 1
        imwrite(frame_n - blank_n - difference_n, strcat('Normalized_',fileprefix,'c',num2str(channelN),'.tif'))
        imwrite(frame_bf - blank_bf - difference_bf, strcat('Normalized_',fileprefix,'c',num2str(channelBF),'.tif'))
        imwrite(frame_r - blank_r - difference_r, strcat('Normalized_',fileprefix,'c',num2str(channelR),'.tif'))
    end
    if k > 1
        imwrite(frame_n - blank_n - difference_n, strcat('Normalized_',fileprefix,'c',num2str(channelN),'.tif'),'WriteMode','append')
        imwrite(frame_bf - blank_bf - difference_bf, strcat('Normalized_',fileprefix,'c',num2str(channelBF),'.tif'),'WriteMode','append')
        imwrite(frame_r - blank_r - difference_r, strcat('Normalized_',fileprefix,'c',num2str(channelR),'.tif'),'WriteMode','append')
    end
end
=======
function [XY] = ImageBlank(fileprefix, channelBF,channelN, channelR)

framerangemin = 1;
framerangemax = 191;
XY = num2str(2, '%02d');

for k = framerangemin:framerangemax
    framenum = num2str(k,'%03d');
    numcentroids = 0;    
    frame_n = imread(strcat(fileprefix,'c',num2str(channelN),'t',framenum,'.tif'));
    blank_n = imread(strcat('blankt1','xy',XY,'c',num2str(channelN),'.tif'));
    difference_n = min(min(frame_n)) - min(min(blank_n));
    
       
    frame_bf = imread(strcat(fileprefix,'c',num2str(channelBF),'t',framenum,'.tif'));
    blank_bf = imread(strcat('blankt1','xy',XY,'c',num2str(channelBF),'.tif'));
    difference_bf = min(min(frame_bf)) - min(min(blank_bf));
    
    frame_r = imread(strcat(fileprefix,'c',num2str(channelR),'t',framenum,'.tif'));
    blank_r = imread(strcat('blankt1','xy',XY,'c',num2str(channelR),'.tif'));
    difference_r = min(min(frame_r)) - min(min(blank_r));

    if k == 1
        imwrite(frame_n - blank_n - difference_n, strcat('Normalized_',fileprefix,'c',num2str(channelN),'.tif'))
        imwrite(frame_bf - blank_bf - difference_bf, strcat('Normalized_',fileprefix,'c',num2str(channelBF),'.tif'))
        imwrite(frame_r - blank_r - difference_r, strcat('Normalized_',fileprefix,'c',num2str(channelR),'.tif'))
    end
    if k > 1
        imwrite(frame_n - blank_n - difference_n, strcat('Normalized_',fileprefix,'c',num2str(channelN),'.tif'),'WriteMode','append')
        imwrite(frame_bf - blank_bf - difference_bf, strcat('Normalized_',fileprefix,'c',num2str(channelBF),'.tif'),'WriteMode','append')
        imwrite(frame_r - blank_r - difference_r, strcat('Normalized_',fileprefix,'c',num2str(channelR),'.tif'),'WriteMode','append')
    end
end
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
