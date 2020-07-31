<<<<<<< HEAD
function[traceclustermat] = graphs(tracks)

maxcell = max(tracks(:,239));
cellvals = [];
fn = [];



hold on
sumvals = zeros(18);
countvals = zeros(18);
integrated = [];
maxlist = [];
peaktime = [];
finallist = [];
act = [];
t20 = [];
hold on
traceclustermat = NaN(max(tracks(:,239)),max(tracks(:,238)));
indexlist = [];

for i = 1:maxcell
    tempindex = tracks(find(tracks(:,239)==i),:);
    %/max(tempindex(:,3))
    tempmovmean = movmean(tempindex(:,3) - min(tempindex(:,3)),3);
    for j = 1:size(tempindex)
        traceclustermat(i,j) = tempmovmean(j);
    end
    
     if tempindex(1,243) == 1
        plot(tempindex(:,238)*5-100,tempmovmean,'Color', [1 0 0 .15])
     end
     if tempindex(1,243) == 0
        plot(tempindex(:,238)*5-100,tempmovmean,'Color', [0 0 1 .15])
     end
%     
%     for j = 1:length(tempindex(:,238))
%         if tempindex(j,238) == 20
% %             t20 = [t20; tempindex(j,3)];
%             maxlist = [maxlist; max(tempindex(20:end,3))];% -min(tempindex(:,3))%/tempindex(j,3)
% %             act = [act; tempindex(j,243)];
% %             peaktime = [peaktime;5*(find(tempindex(:,3) == max(tempindex(20:end,3)),1,'First'))-100];
% 
% %             sumvals = sumvals(tempindex(j,239)) + tempmovmean(j);
% %             countvals(tempindex(j,238)) = countvals(tempindex(j,239)) +1;
%         end
%     end
% %     
%     integrated = [integrated; trapz(5*linspace(1,length(tempmovmean),length(tempmovmean)),tempmovmean/6500)];
%     maxlist = [maxlist; max(tempmovmean)];
     %peaktime = [peaktime; 5*find(tempmovmean == max(tempmovmean))];

    %,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    %'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05
%     active = input('active?');
%     cellvals = [cellvals; active];
%     
%     frame = input('fn?');
%     fn = [fn; frame];
end
finallist = horzcat(peaktime,t20,maxlist,act);
% plot(-15:5:70, sumvals./countvals, 'Color', [1 0 0])
% out = horzcat(cellvals,fn);

=======
function[traceclustermat] = graphs(tracks)

maxcell = max(tracks(:,239));
cellvals = [];
fn = [];



hold on
sumvals = zeros(18);
countvals = zeros(18);
integrated = [];
maxlist = [];
peaktime = [];
finallist = [];
act = [];
t20 = [];
hold on
traceclustermat = NaN(max(tracks(:,239)),max(tracks(:,238)));
indexlist = [];

for i = 1:maxcell
    tempindex = tracks(find(tracks(:,239)==i),:);
    %/max(tempindex(:,3))
    tempmovmean = movmean(tempindex(:,3) - min(tempindex(:,3)),3);
    for j = 1:size(tempindex)
        traceclustermat(i,j) = tempmovmean(j);
    end
    
%     if tempindex(1,243) == 1
    plot(tempindex(:,238)*5-100,tempmovmean,'Color', [1 0 0 .15])
%     end
%     if tempindex(1,243) == 0
%         plot(tempindex(:,238)*5-100,tempmovmean,'Color', [0 0 1 .15])
%     end
%     
%     for j = 1:length(tempindex(:,238))
%         if tempindex(j,238) == 20
% %             t20 = [t20; tempindex(j,3)];
%             maxlist = [maxlist; max(tempindex(20:end,3))];% -min(tempindex(:,3))%/tempindex(j,3)
% %             act = [act; tempindex(j,243)];
% %             peaktime = [peaktime;5*(find(tempindex(:,3) == max(tempindex(20:end,3)),1,'First'))-100];
% 
% %             sumvals = sumvals(tempindex(j,239)) + tempmovmean(j);
% %             countvals(tempindex(j,238)) = countvals(tempindex(j,239)) +1;
%         end
%     end
% %     
%     integrated = [integrated; trapz(5*linspace(1,length(tempmovmean),length(tempmovmean)),tempmovmean/6500)];
%     maxlist = [maxlist; max(tempmovmean)];
     %peaktime = [peaktime; 5*find(tempmovmean == max(tempmovmean))];

    %,'MarkerFaceColor','r','MarkerEdgeColor','r',...
    %'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.05
%     active = input('active?');
%     cellvals = [cellvals; active];
%     
%     frame = input('fn?');
%     fn = [fn; frame];
end
finallist = horzcat(peaktime,t20,maxlist,act);
% plot(-15:5:70, sumvals./countvals, 'Color', [1 0 0])
% out = horzcat(cellvals,fn);

>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
