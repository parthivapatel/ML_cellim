<<<<<<< HEAD
function [allData] = ImageAnalysis(xypos,varargin)
close all

numvarargs = length(varargin);
if numvarargs > 6
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 5 optional inputs');
end

% set defaults for optional inputs (mincentroids, framerangemin,
% framerangemax, outlierminarea, outliermaxarea, rng seed, receptor std threshhold)
optargs = {20 1 75 1000 7000 1 .5};
optargs(1:numvarargs) = varargin;
% Place optional args in memorable variable names
[mincentroids, framerangemin, framerangemax, outliermin, outliermax, seed, receptorthresh] = optargs{:};
OutlierTolerancePercentage = [outliermin, outliermax];

results = cell(1, framerangemax-framerangemin+1);
finaltracks = [];
%centroiddata = cell(1, 70);

% Iterate over all frames
% bfsize
count = 1;
for i = xypos        
    parfor k = framerangemin:framerangemax

        framenum = num2str(k,'%02d');
        numcentroids = 0;

        %im2double
        frame_n = (imread(strcat('2018_02_24_180223_ML_10ngdur_XY',  num2str(i,'%02d'),'_C2.tif'), k));%./Dapi_Blank;
        frame_bf = (imread(strcat('2018_02_24_180223_ML_10ngdur_XY', num2str(i,'%02d'),'_C1.tif'), k));%./BF_Blank;
        frame_r = (imread(strcat('2018_02_24_180223_ML_10ngdur_XY', num2str(i,'%02d'),'_C3.tif'),k));%./YFP_Blank;

        testnucleus = frame_n;
%        testnucleus = imcomplement(testnucleus);

        se = strel('disk',10);
        testnucleus = imdilate(testnucleus, se);
        testnucleus = imerode(testnucleus, se);
        testnucleus = imdilate(testnucleus, se);
        p = testnucleus';
        b = p(:)';  
        Y = sort(b);

        testnucleus(testnucleus<= Y(1000000)) = 0;
        testnucleus = bpass(testnucleus, 2, 51);
        testnucleus = imdilate(testnucleus, se);
        testnucleus = imerode(testnucleus, se);
        testnucleus = imerode(testnucleus, se);

        bw = im2bw(testnucleus);
        DL = bwlabel(bw);
%         figure();imshow(DL)
%         figure;imshow(frame_n, 'InitialMag', 'fit')
    % 
    %     % Make a truecolor all-green image.
    %         green = cat(3, zeros(size(frame_bf)), ones(size(frame_bf)), zeros(size(frame_bf)));
    %         hold on
    %         h = imshow(green);
    %         hold off
    %         set(h, 'AlphaData', DL)
    %         % Show watershed ridgelines
    %         bgm = DL == 0; 
    %         hold on
    %         imshow(bgm), title('Watershed ridge lines (bgm)')
    %         hold off
    %         Lrgb = label2rgb(DL);
    %         figure; imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
    %         hold on
    %         CircleObjectsInImage(DL, 'r');
    %         hold off
    %       %frame_n_overlay = imread(strcat(filenameN,'.tif'),k); %rgb2gray frame_n

        % build properties out of watershed labels, and if the number of
        % centroids is lower than the minimum, repeat. (could change this
        % to a tolerance using standard deviation....)
        objectProperties = regionprops(DL,frame_n, ...
                'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity');

        numcentroids = length(objectProperties);

        % Brightfield properties
        objectProperties_bf = regionprops(DL,frame_bf, ...
                    'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity','MajorAxisLength','MinorAxisLength','Eccentricity', 'EquivDiameter', 'EulerNumber','Orientation','BoundingBox');

        % Reporter properties
        frame_r_overlay = frame_r; % Red channel
        objectProperties_reporter = regionprops(DL,frame_r_overlay, ...
                    'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity' );

        % create a list of outliers based on given tolerances for nuclear area
        medianMaxIntensity = median( [ objectProperties.Area ] );
        outliers = [ objectProperties.Area ] > OutlierTolerancePercentage(2) | [ objectProperties.Area ] < OutlierTolerancePercentage(1);

        % remove outliers from watershed and objectproperites        
        if( sum( outliers ) > 0 )
            DL = RemoveObjectsFromImage( DL, objectProperties( outliers ) );
            objectProperties( outliers ) = [];  % remove outliers
            objectProperties_reporter( outliers ) = [];
            objectProperties_bf( outliers ) = [];
    %        objectProperties_receptor( outliers ) = [];
        end

        
        % grayscale texture features
         texture_featurelist = [];
         for ii = 1:length(objectProperties_bf)
   
            subImage = imcrop(frame_bf, objectProperties_bf(ii).BoundingBox);
             imshow(subImage,[])
            
            Hausdimslope = [hausDim(subImage)];
            otsuthresh = otsurec(subImage, 5);
            sftafeats = sfta(subImage, 2 );

            % different angles: 0,45,90,135
            glcms_offsets = [0 1; -1 1;-1 0;-1 -1];
            glcms = graycomatrix(subImage,'Offset',glcms_offsets,'GrayLimits',[]);
            glcms_statlist = [];

            % iterate through different gray-level co-occurrence matrix at
            % different angles
            for iii = 1:length(glcms_offsets)
                glcms_stats = GLCMFeatures(glcms(:,:,iii));
                glcms_stats = struct2cell(glcms_stats);
                glcms_stats = cell2mat(glcms_stats);
                glcms_statlist = [glcms_statlist, glcms_stats.'];

            end
            % glcmsstats(19 features others)x4,
            % 1hausdimslope,5otsuthresholds,9sftafeats
            texture_featurelist = [texture_featurelist; glcms_statlist, Hausdimslope, otsuthresh.', sftafeats];      
        end
    %    figure;imshow(totalcellbw)   

    
        XY = [];
        for ii = 1:length(objectProperties)
            XY =[XY;i];
        end

        % list of centroids
        centroids = [];
        for ii = 1:length(objectProperties)
            centroids =[centroids; objectProperties(ii).Centroid(1), objectProperties(ii).Centroid(2)];
        end

        % Bf scalar ojbectproperties
        etcinfo = [];
        for ii = 1:length(objectProperties_bf)
            etcinfo =[etcinfo; objectProperties_bf(ii).MajorAxisLength,objectProperties_bf(ii).MinorAxisLength,objectProperties_bf(ii).Eccentricity,...
                objectProperties_bf(ii).EquivDiameter, objectProperties_bf(ii).EulerNumber, objectProperties_bf(ii).Area];
        end

        % add framelabel
        numberOfCentroids = size(centroids, 1);
        centroids = [ centroids, k * ones( numberOfCentroids, 1 ) ];

        % distribution details
        bf_dist = [];
        for ii = 1:length(objectProperties_reporter)
            bf_dist =[bf_dist; mean(double(objectProperties_bf(ii).PixelValues)), std(double(objectProperties_bf(ii).PixelValues)),...
                sum(double(objectProperties_bf(ii).PixelValues))];
        end

        reporter_dist = [];
        for ii = 1:length(objectProperties_reporter)
            reporter_dist =[reporter_dist; mean(double(objectProperties_reporter(ii).PixelValues)), std(double(objectProperties_reporter(ii).PixelValues)), sum(double(objectProperties_reporter(ii).PixelValues))];
        end

        n_dist = [];
        for ii = 1:length(objectProperties_reporter)
            n_dist =[n_dist; mean(double(objectProperties(ii).PixelValues)), std(double(objectProperties(ii).PixelValues)), sum(double(objectProperties(ii).PixelValues))];
        end

        % find angle of ellipse orientation for each nucleus
        orientations = [];
        for ii = 1:length(objectProperties)
            orientations =[orientations; objectProperties_bf(ii).Orientation];
        end

        % Image sift for each centroid (assume each nucleus is same scale and orientation for analysis purposes)
        VL_data = [];
        for ii = 1:length(centroids)
            %orientations(ii) for orientation using ellipse angle instead of 0
            [~,d] = vl_sift(single(frame_bf),'frames',[centroids(ii,1:2).';1;0]);
            VL_data = [VL_data; d.'];
        end

        %receptor_dist,  texture_featurelist_r, texture_featurelist, double(VL_data)
        framedata = [centroids(:,1:2) reporter_dist n_dist bf_dist etcinfo texture_featurelist double(VL_data) XY centroids(:,3)];

        % centroiddata{k} = centroids;
        results{k} = framedata;
        count = count +1;
        
      end
allData = vertcat( results{:} );

end

end

function CircleObjectsInImage( LabelImage, BorderColor )
    boundaries = bwboundaries( LabelImage );	
    numberOfBoundaries = size( boundaries );
    for k = 1 : numberOfBoundaries
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'Color', BorderColor, 'LineWidth', 2);
    end
end

function OutputBinaryImage = RemoveObjectsFromImage( InputBinaryImage, ObjectProperties )
    OutputBinaryImage = InputBinaryImage;
    eliminatedPixels = vertcat( ObjectProperties.PixelList );
    allObjectIndexes = sub2ind( size( InputBinaryImage ), ...
        eliminatedPixels(:, 2), eliminatedPixels(:,1) );
    OutputBinaryImage( allObjectIndexes ) = 0;
end

function out = Gaussian2DFitFunction( Parameters, Coordinates )
    yCenter = Parameters(1);
    xCenter = Parameters(2);
    amplitude = Parameters(3);
    sigma = Parameters(4);
    offset = Parameters(5);
    
    out = amplitude * ...
        exp( -(( Coordinates(:, 1) - yCenter ).^2 + ( Coordinates(:, 2) - xCenter ).^2 ) ...
        ./ (2 * sigma .^ 2 )) + offset;    
end

=======
function [allData] = ImageAnalysis(xypos,varargin)
close all

numvarargs = length(varargin);
if numvarargs > 6
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 5 optional inputs');
end

% set defaults for optional inputs (mincentroids, framerangemin,
% framerangemax, outlierminarea, outliermaxarea, rng seed, receptor std threshhold)
optargs = {20 1 75 1000 7000 1 .5};
optargs(1:numvarargs) = varargin;
% Place optional args in memorable variable names
[mincentroids, framerangemin, framerangemax, outliermin, outliermax, seed, receptorthresh] = optargs{:};
OutlierTolerancePercentage = [outliermin, outliermax];

results = cell(1, framerangemax-framerangemin+1);
finaltracks = [];
%centroiddata = cell(1, 70);

% Iterate over all frames
% bfsize
count = 1;
for i = xypos        
    parfor k = framerangemin:framerangemax

        framenum = num2str(k,'%02d');
        numcentroids = 0;

        %im2double
        frame_n = (imread(strcat('2018_02_24_180223_ML_10ngdur_XY',  num2str(i,'%02d'),'_C2.tif'), k));%./Dapi_Blank;
        frame_bf = (imread(strcat('2018_02_24_180223_ML_10ngdur_XY', num2str(i,'%02d'),'_C1.tif'), k));%./BF_Blank;
        frame_r = (imread(strcat('2018_02_24_180223_ML_10ngdur_XY', num2str(i,'%02d'),'_C3.tif'),k));%./YFP_Blank;

        testnucleus = frame_n;
%        testnucleus = imcomplement(testnucleus);

        se = strel('disk',10);
        testnucleus = imdilate(testnucleus, se);
        testnucleus = imerode(testnucleus, se);
        testnucleus = imdilate(testnucleus, se);
        p = testnucleus';
        b = p(:)';  
        Y = sort(b);

        testnucleus(testnucleus<= Y(1000000)) = 0;
        testnucleus = bpass(testnucleus, 2, 51);
        testnucleus = imdilate(testnucleus, se);
        testnucleus = imerode(testnucleus, se);
        testnucleus = imerode(testnucleus, se);

        bw = im2bw(testnucleus);
        DL = bwlabel(bw);
%         figure();imshow(DL)
%         figure;imshow(frame_n, 'InitialMag', 'fit')
    % 
    %     % Make a truecolor all-green image.
    %         green = cat(3, zeros(size(frame_bf)), ones(size(frame_bf)), zeros(size(frame_bf)));
    %         hold on
    %         h = imshow(green);
    %         hold off
    %         set(h, 'AlphaData', DL)
    %         % Show watershed ridgelines
    %         bgm = DL == 0; 
    %         hold on
    %         imshow(bgm), title('Watershed ridge lines (bgm)')
    %         hold off
    %         Lrgb = label2rgb(DL);
    %         figure; imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
    %         hold on
    %         CircleObjectsInImage(DL, 'r');
    %         hold off
    %       %frame_n_overlay = imread(strcat(filenameN,'.tif'),k); %rgb2gray frame_n

        % build properties out of watershed labels, and if the number of
        % centroids is lower than the minimum, repeat. (could change this
        % to a tolerance using standard deviation....)
        objectProperties = regionprops(DL,frame_n, ...
                'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity');

        numcentroids = length(objectProperties);

        % Brightfield properties
        objectProperties_bf = regionprops(DL,frame_bf, ...
                    'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity','MajorAxisLength','MinorAxisLength','Eccentricity', 'EquivDiameter', 'EulerNumber','Orientation','BoundingBox');

        % Reporter properties
        frame_r_overlay = frame_r; % Red channel
        objectProperties_reporter = regionprops(DL,frame_r_overlay, ...
                    'Area', 'Centroid', 'PixelList', 'PixelValues', 'MaxIntensity' );

        % create a list of outliers based on given tolerances for nuclear area
        medianMaxIntensity = median( [ objectProperties.Area ] );
        outliers = [ objectProperties.Area ] > OutlierTolerancePercentage(2) | [ objectProperties.Area ] < OutlierTolerancePercentage(1);

        % remove outliers from watershed and objectproperites        
        if( sum( outliers ) > 0 )
            DL = RemoveObjectsFromImage( DL, objectProperties( outliers ) );
            objectProperties( outliers ) = [];  % remove outliers
            objectProperties_reporter( outliers ) = [];
            objectProperties_bf( outliers ) = [];
    %        objectProperties_receptor( outliers ) = [];
        end

        
        % grayscale texture features
         texture_featurelist = [];
         for ii = 1:length(objectProperties_bf)
   
            subImage = imcrop(frame_bf, objectProperties_bf(ii).BoundingBox);
             imshow(subImage,[])
            
            Hausdimslope = [hausDim(subImage)];
            otsuthresh = otsurec(subImage, 5);
            sftafeats = sfta(subImage, 2 );

            % different angles: 0,45,90,135
            glcms_offsets = [0 1; -1 1;-1 0;-1 -1];
            glcms = graycomatrix(subImage,'Offset',glcms_offsets,'GrayLimits',[]);
            glcms_statlist = [];

            % iterate through different gray-level co-occurrence matrix at
            % different angles
            for iii = 1:length(glcms_offsets)
                glcms_stats = GLCMFeatures(glcms(:,:,iii));
                glcms_stats = struct2cell(glcms_stats);
                glcms_stats = cell2mat(glcms_stats);
                glcms_statlist = [glcms_statlist, glcms_stats.'];

            end
            % glcmsstats(19 features others)x4,
            % 1hausdimslope,5otsuthresholds,9sftafeats
            texture_featurelist = [texture_featurelist; glcms_statlist, Hausdimslope, otsuthresh.', sftafeats];      
        end
    %    figure;imshow(totalcellbw)   

    
        XY = [];
        for ii = 1:length(objectProperties)
            XY =[XY;i];
        end

        % list of centroids
        centroids = [];
        for ii = 1:length(objectProperties)
            centroids =[centroids; objectProperties(ii).Centroid(1), objectProperties(ii).Centroid(2)];
        end

        % Bf scalar ojbectproperties
        etcinfo = [];
        for ii = 1:length(objectProperties_bf)
            etcinfo =[etcinfo; objectProperties_bf(ii).MajorAxisLength,objectProperties_bf(ii).MinorAxisLength,objectProperties_bf(ii).Eccentricity,...
                objectProperties_bf(ii).EquivDiameter, objectProperties_bf(ii).EulerNumber, objectProperties_bf(ii).Area];
        end

        % add framelabel
        numberOfCentroids = size(centroids, 1);
        centroids = [ centroids, k * ones( numberOfCentroids, 1 ) ];

        % distribution details
        bf_dist = [];
        for ii = 1:length(objectProperties_reporter)
            bf_dist =[bf_dist; mean(double(objectProperties_bf(ii).PixelValues)), std(double(objectProperties_bf(ii).PixelValues)),...
                sum(double(objectProperties_bf(ii).PixelValues))];
        end

        reporter_dist = [];
        for ii = 1:length(objectProperties_reporter)
            reporter_dist =[reporter_dist; mean(double(objectProperties_reporter(ii).PixelValues)), std(double(objectProperties_reporter(ii).PixelValues)), sum(double(objectProperties_reporter(ii).PixelValues))];
        end

        n_dist = [];
        for ii = 1:length(objectProperties_reporter)
            n_dist =[n_dist; mean(double(objectProperties(ii).PixelValues)), std(double(objectProperties(ii).PixelValues)), sum(double(objectProperties(ii).PixelValues))];
        end

        % find angle of ellipse orientation for each nucleus
        orientations = [];
        for ii = 1:length(objectProperties)
            orientations =[orientations; objectProperties_bf(ii).Orientation];
        end

        % Image sift for each centroid (assume each nucleus is same scale and orientation for analysis purposes)
        VL_data = [];
        for ii = 1:length(centroids)
            %orientations(ii) for orientation using ellipse angle instead of 0
            [~,d] = vl_sift(single(frame_bf),'frames',[centroids(ii,1:2).';1;0]);
            VL_data = [VL_data; d.'];
        end

        %receptor_dist,  texture_featurelist_r, texture_featurelist, double(VL_data)
        framedata = [centroids(:,1:2) reporter_dist n_dist bf_dist etcinfo texture_featurelist double(VL_data) XY centroids(:,3)];

        % centroiddata{k} = centroids;
        results{k} = framedata;
        count = count +1;
        
      end
allData = vertcat( results{:} );

end

end

function CircleObjectsInImage( LabelImage, BorderColor )
    boundaries = bwboundaries( LabelImage );	
    numberOfBoundaries = size( boundaries );
    for k = 1 : numberOfBoundaries
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'Color', BorderColor, 'LineWidth', 2);
    end
end

function OutputBinaryImage = RemoveObjectsFromImage( InputBinaryImage, ObjectProperties )
    OutputBinaryImage = InputBinaryImage;
    eliminatedPixels = vertcat( ObjectProperties.PixelList );
    allObjectIndexes = sub2ind( size( InputBinaryImage ), ...
        eliminatedPixels(:, 2), eliminatedPixels(:,1) );
    OutputBinaryImage( allObjectIndexes ) = 0;
end

function out = Gaussian2DFitFunction( Parameters, Coordinates )
    yCenter = Parameters(1);
    xCenter = Parameters(2);
    amplitude = Parameters(3);
    sigma = Parameters(4);
    offset = Parameters(5);
    
    out = amplitude * ...
        exp( -(( Coordinates(:, 1) - yCenter ).^2 + ( Coordinates(:, 2) - xCenter ).^2 ) ...
        ./ (2 * sigma .^ 2 )) + offset;    
end

>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
