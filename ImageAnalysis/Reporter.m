<<<<<<< HEAD
function [reportertracks] = Reporter(tracks)
close all
timescale = 5;

%tracks(end)
reporterdist = [];
reporterIndicatorList = [];
TimeshiftList = [];
reportertracks = [0,0,0,0];

reporterinfo = {};
reporterinfolist = [];
reporteravg = [];
for k = 1:tracks(end) 
    particleNData = tracks(tracks(:,end) == k, : );
    
    Trajectory=particleNData(:,1:2);
    TauValues=particleNData(:,end-1);
    ReporterMeans = particleNData(:,3);
    ReporterMax = max(ReporterMeans);
    ReporterSum = particleNData(:,5);
    
%    plot(movmean(TauValues,5),(movmean(ReporterMeans,5)-min(ReporterMeans))/(max(ReporterMeans)-min(ReporterMeans)))
    close all
    plot(movmean(TauValues,5),movmean(ReporterMeans,5))
    active = input('active?');
    Timeshift = 0;
    maxreporterindex = 0;
    reporterIndicator = 0;
    if ReporterMax >=400
        maxreporterindex = find(ReporterMeans == ReporterMax);
        Timeshift = 70- TauValues(maxreporterindex);
        reporterIndicator = 1;
    end
    
    TimeshiftList = [TimeshiftList; Timeshift];
    reporterIndicatorList = [reporterIndicatorList;reporterIndicator];
    reporterinfo = [ones(length(movmean(TauValues,5)),1).*k, movmean(TauValues,5),movmean(ReporterMeans,5)];
    reporterinfolist = [reporterinfolist; reporterinfo];
    TauValueShift = TauValues +(Timeshift);
    ReporterMin = median(ReporterMeans);
    
    reporterdist = [reporterdist; (ReporterMax - ReporterMin)/ReporterMin];
    
    numtimesteps = length(ReporterMeans);
    
    timebin = 0;
    if Timeshift >1 &&  Timeshift < 15
        timebin = 1;
    end
    if Timeshift >=15 &&  Timeshift < 30
        timebin = 2;
    end
    if Timeshift >=30 &&  Timeshift < 45
        timebin = 3;
    end
    if Timeshift >=45 &&  Timeshift < 60
        timebin = 4;
    end
    if Timeshift >=60
        timebin = 5;
    end
    
    
    dummyreporterindicator = ones(numtimesteps,1)*reporterIndicator;
    dummytimeshift = ones(numtimesteps,1)*Timeshift;
    dummytimebin = ones(numtimesteps,1)*timebin;
    activelist = ones(numtimesteps,1)*active;
    
    
    results = [dummyreporterindicator dummytimeshift dummytimebin, activelist];

    reportertracks = vertcat(reportertracks, results);
    
    %Msd = MeanSquaredDisplacement( Trajectory, TauValues );
    %figure(1); loglog(TauValues*timescale,Msd)
    %figure();plot(TauValues,ReporterMeans,'-')
    %


end
reporteravg = [];
maxtime = max(reporterinfo(:,2));

timeline = [1:.2:maxtime];
for i = timeline
%    timepoints = find(reporterinfo(2,:)< i+.2 & reporterinfo(2,:) >= i);
    tempavg = reporterinfolist(reporterinfolist(:,2)>i,:);
    tempavg = tempavg(tempavg(:,2)<=i+1,:);
    if length(tempavg(:,1)) > 4
        reporteravg = [reporteravg; [i,mean(tempavg(:,3))]];
    end
end
figure()
reporteravg(any(isnan(reporteravg),2),:) = [];


plot(reporteravg(:,1).*5,reporteravg(:,2))
% reportertracks = [tracks reportertracks(2:end,:)];
% reportertracks(imag(reportertracks) ~= 0) = NaN;
% reportertracks = reportertracks(~any(isnan(reportertracks), 2), :);
% reportertracks = reportertracks(reportertracks(:,3) < 160, :);
% 
% reportertracks = reportertracks(reportertracks(:,end-4) < 6, :);
%reportertracks = reportertracks(reportertracks(:,end-4) < 60, :);
% reportertracks = reportertracks(reportertracks(:,end-3) < 30, :);
% reportertracks = reportertracks(reportertracks(:,end-3) > 5, :);


%{
final_reportertracks = [0, 0, 0];
for k = 1:reportertracks(end,end-3) 
    particleNData = reportertracks(reportertracks(:,end-3) == k, : );
    Trajectory=particleNData(:,1:2);
    timedata = particleNData(:,end-4);
    timesize = size( particleNData);
    %find 3 closest neighbors
    distance1 = ones(timesize(1),1)*10000;
    index1 = ones(timesize(1),1)*k;
    reporter1 = zeros(timesize(1),1);
    distance2 = ones(timesize(1),1)*10000;
    index2 = ones(timesize(1),1)*k;
    reporter2 = zeros(timesize(1),1);
    distance3 = ones(timesize(1),1)*10000;
    index3 = ones(timesize(1),1)*k;
    reporter3 = zeros(timesize(1),1);
    for ii = 1:reportertracks(end)
        if ii ~= k
            XTrajectory = reportertracks(reportertracks(:,end-3) == ii, 1:2 );
            particleXtime = reportertracks(reportertracks(:,end-3) == ii, end-4 );
            Xreporter = reportertracks(reportertracks(:,end-3) == ii, end-2 );
            particleXTrajectory = ones(timesize(1),2)*100000;
            particleXReporter = zeros(timesize(1),1);
            for t = 1:length(particleXtime)
                timeXindex = find(timedata == particleXtime(t));
                
                if ~isempty(timeXindex);
                    particleXTrajectory(timeXindex,1) = XTrajectory(t,1);
                    particleXTrajectory(timeXindex,2) = XTrajectory(t,2);
                    particleXReporter(timeXindex) = Xreporter(t);
                end
            end
            distance = sqrt((Trajectory(:,1)-particleXTrajectory(:,1)).^2 + (Trajectory(:,2)-particleXTrajectory(:,2)).^2);
            reporter = particleXReporter./distance;
            for iii = 1:length(distance)
                if distance(iii) < distance3(iii)
                    distance3(iii) = distance(iii); 
                    index3(iii) = ii;
                    reporter3(iii) = reporter(iii);
                    if distance(iii) < distance2(iii)
                        distance3(iii) = distance2(iii);
                        index3(iii) = index2(iii);
                        reporter3(iii) = reporter2(iii);
                        distance2(iii) = distance(iii);
                        index2(iii) = ii;
                        reporter2(iii) = reporter(iii);
                        if distance(iii) < distance1(iii)
                            distance2(iii) = distance1(iii);
                            index2(iii) = index1(iii);
                            reporter2(iii) = reporter1(iii);
                            distance1(iii) = distance(iii);
                            index1(iii) = ii;
                            reporter1(iii) = reporter(iii);
                        end
                    end
                end
            end
        end
    end
    distanceresults = [reporter1 reporter2 reporter3];
    final_reportertracks = vertcat(final_reportertracks, distanceresults);
end

reportertracks = [reportertracks(:,1:end-5) final_reportertracks(2:end,:) reportertracks(:,end-4:end)];

reportersize = size(reportertracks);
%}
%{
derivresults = zeros((reportersize(2)-4),1).';

% derivitives
for k = 1:reportertracks(end,end-3) 
    particleNData = reportertracks(reportertracks(:,end-3) == k, : );
    timedata = particleNData(:,end-4);
    if length(timedata) >= 1
        particledevdata = zeros((reportersize(2)-4),1).';
    end
    if length(timedata)>1
        for ii = 2:length(timedata)
            particledevdata = [particledevdata;(particleNData(ii,1:end-4)-particleNData(ii-1,1:end-4))/(timedata(ii)-timedata(ii-1))];
        end
    end
    if length(timedata) >= 1
        derivresults = vertcat(derivresults, particledevdata);   
    end
end
reportertracks = [reportertracks(:,1:end-5) derivresults(2:end,1:8) derivresults(2:end,137:end) reportertracks(:,end-4:end)];
%}

% singlecelllist=[];
% for k = 1:reportertracks(end,end-3) 
%     particleNData = reportertracks(reportertracks(:,end-3) == k, : );
%     if ~isempty(particleNData)
%         singlecell = particleNData(end,:);
%         singlecelllist = [singlecelllist;singlecell];
%     end
% end
% reportertracks = singlecelllist;

end
function [ Msd] = MeanSquaredDisplacement( Trajectory, TauValues )
    trajectoryLength = size( Trajectory, 1 );
    Msd = zeros( 1, length( TauValues ) );
    
    for kk = 1:length( TauValues )
        interval = TauValues(kk);
        startIndex = 1:(trajectoryLength - interval);
        endIndex = startIndex + interval;
        
        Msd(kk) = mean( sum( (Trajectory( endIndex, : ) - Trajectory(  startIndex, : ) ).^2, 2 ) );
    end
=======
function [reportertracks] = Reporter(tracks)
close all
timescale = 5;

%tracks(end)
reporterdist = [];
reporterIndicatorList = [];
TimeshiftList = [];
reportertracks = [0,0,0,0];

reporterinfo = {};
reporterinfolist = [];
reporteravg = [];
for k = 1:tracks(end) 
    particleNData = tracks(tracks(:,end) == k, : );
    
    Trajectory=particleNData(:,1:2);
    TauValues=particleNData(:,end-1);
    ReporterMeans = particleNData(:,3);
    ReporterMax = max(ReporterMeans);
    ReporterSum = particleNData(:,5);
    
%    plot(movmean(TauValues,5),(movmean(ReporterMeans,5)-min(ReporterMeans))/(max(ReporterMeans)-min(ReporterMeans)))
    close all
    plot(movmean(TauValues,5),movmean(ReporterMeans,5))
    active = input('active?');
    Timeshift = 0;
    maxreporterindex = 0;
    reporterIndicator = 0;
    if ReporterMax >=400
        maxreporterindex = find(ReporterMeans == ReporterMax);
        Timeshift = 70- TauValues(maxreporterindex);
        reporterIndicator = 1;
    end
    
    TimeshiftList = [TimeshiftList; Timeshift];
    reporterIndicatorList = [reporterIndicatorList;reporterIndicator];
    reporterinfo = [ones(length(movmean(TauValues,5)),1).*k, movmean(TauValues,5),movmean(ReporterMeans,5)];
    reporterinfolist = [reporterinfolist; reporterinfo];
    TauValueShift = TauValues +(Timeshift);
    ReporterMin = median(ReporterMeans);
    
    reporterdist = [reporterdist; (ReporterMax - ReporterMin)/ReporterMin];
    
    numtimesteps = length(ReporterMeans);
    
    timebin = 0;
    if Timeshift >1 &&  Timeshift < 15
        timebin = 1;
    end
    if Timeshift >=15 &&  Timeshift < 30
        timebin = 2;
    end
    if Timeshift >=30 &&  Timeshift < 45
        timebin = 3;
    end
    if Timeshift >=45 &&  Timeshift < 60
        timebin = 4;
    end
    if Timeshift >=60
        timebin = 5;
    end
    
    
    dummyreporterindicator = ones(numtimesteps,1)*reporterIndicator;
    dummytimeshift = ones(numtimesteps,1)*Timeshift;
    dummytimebin = ones(numtimesteps,1)*timebin;
    activelist = ones(numtimesteps,1)*active;
    
    
    results = [dummyreporterindicator dummytimeshift dummytimebin, activelist];

    reportertracks = vertcat(reportertracks, results);
    
    %Msd = MeanSquaredDisplacement( Trajectory, TauValues );
    %figure(1); loglog(TauValues*timescale,Msd)
    %figure();plot(TauValues,ReporterMeans,'-')
    %


end
reporteravg = [];
maxtime = max(reporterinfo(:,2));

timeline = [1:.2:maxtime];
for i = timeline
%    timepoints = find(reporterinfo(2,:)< i+.2 & reporterinfo(2,:) >= i);
    tempavg = reporterinfolist(reporterinfolist(:,2)>i,:);
    tempavg = tempavg(tempavg(:,2)<=i+1,:);
    if length(tempavg(:,1)) > 4
        reporteravg = [reporteravg; [i,mean(tempavg(:,3))]];
    end
end
figure()
reporteravg(any(isnan(reporteravg),2),:) = [];


plot(reporteravg(:,1).*5,reporteravg(:,2))
% reportertracks = [tracks reportertracks(2:end,:)];
% reportertracks(imag(reportertracks) ~= 0) = NaN;
% reportertracks = reportertracks(~any(isnan(reportertracks), 2), :);
% reportertracks = reportertracks(reportertracks(:,3) < 160, :);
% 
% reportertracks = reportertracks(reportertracks(:,end-4) < 6, :);
%reportertracks = reportertracks(reportertracks(:,end-4) < 60, :);
% reportertracks = reportertracks(reportertracks(:,end-3) < 30, :);
% reportertracks = reportertracks(reportertracks(:,end-3) > 5, :);


%{
final_reportertracks = [0, 0, 0];
for k = 1:reportertracks(end,end-3) 
    particleNData = reportertracks(reportertracks(:,end-3) == k, : );
    Trajectory=particleNData(:,1:2);
    timedata = particleNData(:,end-4);
    timesize = size( particleNData);
    %find 3 closest neighbors
    distance1 = ones(timesize(1),1)*10000;
    index1 = ones(timesize(1),1)*k;
    reporter1 = zeros(timesize(1),1);
    distance2 = ones(timesize(1),1)*10000;
    index2 = ones(timesize(1),1)*k;
    reporter2 = zeros(timesize(1),1);
    distance3 = ones(timesize(1),1)*10000;
    index3 = ones(timesize(1),1)*k;
    reporter3 = zeros(timesize(1),1);
    for ii = 1:reportertracks(end)
        if ii ~= k
            XTrajectory = reportertracks(reportertracks(:,end-3) == ii, 1:2 );
            particleXtime = reportertracks(reportertracks(:,end-3) == ii, end-4 );
            Xreporter = reportertracks(reportertracks(:,end-3) == ii, end-2 );
            particleXTrajectory = ones(timesize(1),2)*100000;
            particleXReporter = zeros(timesize(1),1);
            for t = 1:length(particleXtime)
                timeXindex = find(timedata == particleXtime(t));
                
                if ~isempty(timeXindex);
                    particleXTrajectory(timeXindex,1) = XTrajectory(t,1);
                    particleXTrajectory(timeXindex,2) = XTrajectory(t,2);
                    particleXReporter(timeXindex) = Xreporter(t);
                end
            end
            distance = sqrt((Trajectory(:,1)-particleXTrajectory(:,1)).^2 + (Trajectory(:,2)-particleXTrajectory(:,2)).^2);
            reporter = particleXReporter./distance;
            for iii = 1:length(distance)
                if distance(iii) < distance3(iii)
                    distance3(iii) = distance(iii); 
                    index3(iii) = ii;
                    reporter3(iii) = reporter(iii);
                    if distance(iii) < distance2(iii)
                        distance3(iii) = distance2(iii);
                        index3(iii) = index2(iii);
                        reporter3(iii) = reporter2(iii);
                        distance2(iii) = distance(iii);
                        index2(iii) = ii;
                        reporter2(iii) = reporter(iii);
                        if distance(iii) < distance1(iii)
                            distance2(iii) = distance1(iii);
                            index2(iii) = index1(iii);
                            reporter2(iii) = reporter1(iii);
                            distance1(iii) = distance(iii);
                            index1(iii) = ii;
                            reporter1(iii) = reporter(iii);
                        end
                    end
                end
            end
        end
    end
    distanceresults = [reporter1 reporter2 reporter3];
    final_reportertracks = vertcat(final_reportertracks, distanceresults);
end

reportertracks = [reportertracks(:,1:end-5) final_reportertracks(2:end,:) reportertracks(:,end-4:end)];

reportersize = size(reportertracks);
%}
%{
derivresults = zeros((reportersize(2)-4),1).';

% derivitives
for k = 1:reportertracks(end,end-3) 
    particleNData = reportertracks(reportertracks(:,end-3) == k, : );
    timedata = particleNData(:,end-4);
    if length(timedata) >= 1
        particledevdata = zeros((reportersize(2)-4),1).';
    end
    if length(timedata)>1
        for ii = 2:length(timedata)
            particledevdata = [particledevdata;(particleNData(ii,1:end-4)-particleNData(ii-1,1:end-4))/(timedata(ii)-timedata(ii-1))];
        end
    end
    if length(timedata) >= 1
        derivresults = vertcat(derivresults, particledevdata);   
    end
end
reportertracks = [reportertracks(:,1:end-5) derivresults(2:end,1:8) derivresults(2:end,137:end) reportertracks(:,end-4:end)];
%}

% singlecelllist=[];
% for k = 1:reportertracks(end,end-3) 
%     particleNData = reportertracks(reportertracks(:,end-3) == k, : );
%     if ~isempty(particleNData)
%         singlecell = particleNData(end,:);
%         singlecelllist = [singlecelllist;singlecell];
%     end
% end
% reportertracks = singlecelllist;

end
function [ Msd] = MeanSquaredDisplacement( Trajectory, TauValues )
    trajectoryLength = size( Trajectory, 1 );
    Msd = zeros( 1, length( TauValues ) );
    
    for kk = 1:length( TauValues )
        interval = TauValues(kk);
        startIndex = 1:(trajectoryLength - interval);
        endIndex = startIndex + interval;
        
        Msd(kk) = mean( sum( (Trajectory( endIndex, : ) - Trajectory(  startIndex, : ) ).^2, 2 ) );
    end
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
end