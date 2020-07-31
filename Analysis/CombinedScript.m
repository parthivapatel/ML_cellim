<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
%% Use track program to track centroids
% param=struct(); %structure array
% param.mem = 2; %can leave for 10 frames
% param.dim = 2;
% param.good = 10; %set high when code is working (max is number of frames)
% param.quiet = 1;
% maxdisplacement = 100;

%example: Use ImageAnalysis function as a function of XYposition, then track using track function
%xy01 = ImageAnalysis(1);
%med01 = track(xy01,maxdisplacement,param);
%graphs(med01)

%% Combine dosages and isolate t=0 for each stimulation

% 
% maxlist = []
% lowlist = vertcat(low01,low02,low03,low04,low05,low06,low07,low08,low09,low10,low11,low12,low13,low14,low15);
% lowdata7 = lowlist(find(lowlist(:,238)==7),:);
% lowdata21 = lowlist(find(lowlist(:,238)==21),:);
% maxlist = [graphs(low01);graphs(low02);graphs(low03);graphs(low04);graphs(low05);graphs(low06);graphs(low07);graphs(low08);graphs(low09);graphs(low10);graphs(low11);graphs(low12);graphs(low13);graphs(low14);graphs(low15)];
% lowdataclean = horzcat(lowdata7,maxlist)
% lowdataclean = lowdataclean(find(lowdataclean(:,240)>0),:);
% 
% 
% maxlist = []
% medlist = vertcat(med01,med02,med03,med04,med05,med06,med07,med08,med09,med10,med11,med12,med13,med14,med15);
% meddata7 = medlist(find(medlist(:,238)==7),:);
% meddata21 = medlist(find(medlist(:,238)==21),:);
% maxlist = [graphs(med01);graphs(med02);graphs(med03);graphs(med04);graphs(med05);graphs(med06);graphs(med07);graphs(med08);graphs(med09);graphs(med10);graphs(med11);graphs(med12);graphs(med13);graphs(med14);graphs(med15)];
% meddataclean = horzcat(meddata7,maxlist)
% meddataclean = meddataclean(find(meddataclean(:,240)>0),:);
% 
% 
% maxlist = []
% highlist = vertcat(high01,high02,high03,high04,high05,high06,high07,high08,high09,high10,high11,high12,high13,high15);
% highdata7 = highlist(find(highlist(:,238)==7),:);
% highdata21 = highlist(find(highlist(:,238)==21),:);
% maxlist = [graphs(high01);graphs(high02);graphs(high03);graphs(high04);graphs(high05);graphs(high06);graphs(high07);graphs(high08);graphs(high09);graphs(high10);graphs(high11);graphs(high12);graphs(high13);graphs(high15)];
% highdataclean = horzcat(highdata7,maxlist)
% highdataclean = highdataclean(find(highdataclean(:,240)>0),:);
% 
% maxlist = []
% vhlist = vertcat(vh01,vh02,vh03,vh04,vh05,vh06,vh07,vh08,vh09,vh10,vh11);
% vhdata7 = vhlist(find(vhlist(:,238)==7),:);
% vhdata21 = vhlist(find(vhlist(:,238)==21),:);
% maxlist = [graphs(vh01);graphs(vh02);graphs(vh03);graphs(vh04);graphs(vh05);graphs(vh06);graphs(vh07);graphs(vh08);graphs(vh09);graphs(vh10);graphs(vh11)];
% vhdataclean = horzcat(vhdata7,maxlist)
% vhdataclean = vhdataclean(find(vhdataclean(:,240)>0),:);
% 
% lowdataclean = horzcat(lowdataclean,ones(size(lowdataclean,1),1));
% meddataclean = horzcat(meddataclean,2*ones(size(meddataclean,1),1));
% highdataclean = horzcat(highdataclean,3*ones(size(highdataclean,1),1));
% vhdataclean = horzcat(vhdataclean,4*ones(size(vhdataclean,1),1));

load("MatlabWorkspaceBase.mat")

%Include ikb stain data
ikbstain = vertcat(xy13,xy14,xy15,xy16,xy17,xy18);
ikbstain13 = ikbstain(find(ikbstain(:,238)==13),:);

%Normalize
mmlow = minmax(lowdataclean');
mmmed = minmax(meddataclean');
mmhigh = minmax(highdataclean');
mmvh = minmax(vhdataclean');
mmikb = minmax(ikbstain13');


lowdatanorm = [];
meddatanorm = [];
highdatanorm = [];
vhdatanorm = [];
ikbnorm = [];
for i = 1:241
    if i <237
        lowdatanorm = [lowdatanorm,(lowdataclean(:,i)-mmlow(i,1))/mmlow(i,2)];
        meddatanorm = [meddatanorm,(meddataclean(:,i)-mmmed(i,1))/mmmed(i,2)];
        highdatanorm = [highdatanorm,(highdataclean(:,i)-mmhigh(i,1))/mmhigh(i,2)];
        vhdatanorm = [vhdatanorm,(vhdataclean(:,i)-mmvh(i,1))/mmvh(i,2)];
        ikbnorm = [ikbnorm, (ikbstain13(:,i)-mmikb(i,1))/mmikb(i,2)];
    end
end


alldosedataclean = vertcat(lowdataclean,meddataclean,highdataclean,vhdataclean);
alldosedata = vertcat(lowdatanorm,meddatanorm,highdatanorm,vhdatanorm);
alldosedataclean = [alldosedataclean, alldosedataclean(:,240)-alldosedataclean(:,3)];
alldosedataclean_active = alldosedataclean(find(alldosedataclean(:,242)> 500),:);
activedummy = zeros(size(alldosedataclean,1),1);
activedummy(find(alldosedataclean(:,242)> 200)) = 1;



%% tsne
rng(1)
tsneplotdose = tsne(alldosedata(:,3:15),'Algorithm','barneshut');
figure()
gscatter(tsneplotdose(:,1),tsneplotdose(:,2),alldosedataclean(:,241))
alldosedata_ikb = vertcat(lowdatanorm,meddatanorm,highdatanorm,vhdatanorm,ikbnorm);
alldosedata_norm = vertcat(lowdatanorm,meddatanorm,highdatanorm,vhdatanorm);
alldosedata_norm = horzcat(alldosedata_norm, alldosedataclean(:,241), activedummy);


%rng1
rng(1)
tsneplotdose_ikb_1 = tsne(alldosedata_ikb(:,3:15),'Algorithm','barneshut','perplexity',800,'Exaggeration', 1);

x0=10;
y0=10;

figure()
scatter(tsneplotdose_ikb_1(1:3456,1),tsneplotdose_ikb_1(1:3456,2), 3, log(alldosedataclean(:,240)));
colormap(flipud(hot))
caxis([6 9])
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
hold on

predictmodel = [];
for i = 1:3456
    predictmodel = [predictmodel; trainedModel.predictFcn(horzcat(alldosedataclean(i,3:236),alldosedataclean(i,241)))];
end


%% Embedding with UMAP

figure()
hold on
scatter(embedding(find(predictmodel==1),1),embedding(find(predictmodel==1),2), 4, [1 0 0],'filled');
scatter(embedding(find(predictmodel==0),1),embedding(find(predictmodel==0),2), 4, [.5 .5 .5],'filled');
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('prediction')


% correlation in nucleus texture feature for resistent group
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,alldosedataclean(:,22),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('correlation in nucleus texture feature for resistent group')


% consistency in nucleus (entropy)
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,alldosedataclean(:,23),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('consistency in nucleus (entropy)')


%Otsu dim
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,log(alldosedataclean(:,95)),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('Otsu dim')

%INFORMATION MEASURE OF CORR
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,alldosedataclean(:,29),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('Information measure of correlation')

%leaky
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,log(alldosedataclean(:,3)),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('leaky nfkb')

%leaky std
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,log(alldosedataclean(:,4)),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('leaky nfkb std in nucleus')


%% Embedding breakdown by group
figure()
hold on
plot(embedding(find(alldosedataclean(:,241) == 1),1),embedding(find(alldosedataclean(:,241) == 1),2),'r','LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 2),1),embedding(find(alldosedataclean(:,241) == 2),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 3),1),embedding(find(alldosedataclean(:,241) == 3),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 4),1),embedding(find(alldosedataclean(:,241) == 4),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])
yticks([-50 0 50])


figure()
hold on
plot(embedding(find(alldosedataclean(:,241) == 1),1),embedding(find(alldosedataclean(:,241) == 1),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 2),1),embedding(find(alldosedataclean(:,241) == 2),2),'r','LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 3),1),embedding(find(alldosedataclean(:,241) == 3),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 4),1),embedding(find(alldosedataclean(:,241) == 4),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])
yticks([-50 0 50])



figure()
hold on
plot(embedding(find(alldosedataclean(:,241) == 1),1),embedding(find(alldosedataclean(:,241) == 1),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 2),1),embedding(find(alldosedataclean(:,241) == 2),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 3),1),embedding(find(alldosedataclean(:,241) == 3),2),'r','LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 4),1),embedding(find(alldosedataclean(:,241) == 4),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])
yticks([-50 0 50])

figure()
hold on
plot(embedding(find(alldosedataclean(:,241) == 1),1),embedding(find(alldosedataclean(:,241) == 1),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 2),1),embedding(find(alldosedataclean(:,241) == 2),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 3),1),embedding(find(alldosedataclean(:,241) == 3),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 4),1),embedding(find(alldosedataclean(:,241) == 4),2),'Color','r','LineStyle','none','Marker','o','MarkerSize', 1)
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])
yticks([-50 0 50])

%% Finding high probability features


pvaluestest1 = [];
logfoldchangetest1 = [];
for i = 1:238
    [~,p] = ttest2(alldosedataclean(find(predictmodel(:)==0),i),alldosedataclean(find(predictmodel(:)==1),i));
    logfoldchange = mean(alldosedataclean(find(predictmodel(:)==0),i))/mean(alldosedataclean(find(predictmodel(:)==1),i));
    pvaluestest1 = [pvaluestest1,p];
    logfoldchangetest1 = [logfoldchangetest1,logfoldchange];
end
figure(); hold on
scatter(-log(logfoldchangetest1),log10(-log10(pvaluestest1)),3,'k')
scatter(-log(logfoldchangetest1([99:107])),log10(-log10(pvaluestest1([99:107]))),25,'r')


pvaluessig = find(log10(-log10(pvaluestest1))>1);
[~,indexpvalues] = sort(-log(logfoldchangetest1));

corrmatrix = [];
for i = 1:238
    corrmatrix = [corrmatrix,corr(alldosedataclean(:,i),predictmodel)];
end
[~,indexcorr] = sort(corrmatrix);


corrsorted = corrmatrix(indexcorr);

figure(); hold on
barh(1:20,corrsorted(1:20),'FaceColor',[.2,0.2,0.2],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)
barh(21:41,corrsorted(211:231),'FaceColor',[1,0,0],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)



%% Morphology velocity
medlist_test = medlist(find(medlist(:,238)<10),:);
lowlist_test = lowlist(find(lowlist(:,238)<10),:);

[knnidx,knnweight] = knnsearch(alldosedataclean(:,3:15),medlist_test(:,3:15),'K', 2);
knnemb = [];
for i = 1:size(knnidx,1)
    tempmat = (1./knnweight(i,:))./(sum(1./knnweight(i,:)));
    knnemb = [knnemb; mean(2*tempmat'.*embedding(knnidx(i,:),:))];
end

%figure()
%scatter(embedding(:,1),embedding(:,2),3)
%hold on; scatter(knnemb(:,1),knnemb(:,2),2)

figure()
f1 = scatter(embedding(:,1),embedding(:,2),5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
hold on
countval =0;
for j = 1:max(medlist_test(:,237))
    maxcell = max(medlist_test(find(medlist_test(:,237) == j),239));
    for i = 1:maxcell
        %tempindex = tracks(find(tracks(:,239)==i),:);
        %tempmovmean = movmean(tempindex(:,3) - min(tempindex(:,3)),3);
        plot(knnemb(find(medlist_test(:,239)==i & medlist_test(:,237) == j),1),knnemb(find(medlist_test(:,239)==i & medlist_test(:,237) == j),2))
    end
end

axis([-8 7 -7 7])
axis off

% cell3=medlist_test(find(medlist_test(:,239)==63 & medlist_test(:,237) == 2),:);
% cell4=medlist_test(find(medlist_test(:,239)==16 & medlist_test(:,237) == 4),:);
% cell1=lowlist_test(find(lowlist_test(:,239)==9 & lowlist_test(:,237) == 7),:);
% cell2=medlist_test(find(medlist_test(:,239)==19 & medlist_test(:,237) == 15),:);
% 
% figure()
% clf
% hold on
% plot(-30:5:0,movmean(cell1(1:7,3),3),'r')
% plot(-30:5:0,movmean(cell2(1:7,3),3),'y')
% plot(-30:5:0,movmean(cell3(1:7,3),3),'m')
% plot(-30:5:0,movmean(cell4(1:7,3),3),'g')
% 
% %pixel entropy
% figure()
% clf
% hold on
% plot(-30:5:0,movmean(cell1(1:7,27),3),'r')
% plot(-30:5:0,movmean(cell2(1:7,27),3),'y')
% plot(-30:5:0,movmean(cell3(1:7,27),3),'m')
% plot(-30:5:0,movmean(cell4(1:7,27),3),'g')
% 
% 
% %pixel correlation
% figure()
% clf
% hold on
% plot(-30:5:0,movmean(cell1(1:7,22),3),'r')
% plot(-30:5:0,movmean(cell2(1:7,22),3),'y')
% plot(-30:5:0,movmean(cell3(1:7,22),3),'m')
% plot(-30:5:0,movmean(cell4(1:7,22),3),'g')
% 
% %information measure
% figure()
% clf
% hold on
% plot(-30:5:0,movmean(cell1(1:7,67),3),'r')
% plot(-30:5:0,movmean(cell2(1:7,67),3),'y')
% plot(-30:5:0,movmean(cell3(1:7,67),3),'m')
% plot(-30:5:0,movmean(cell4(1:7,67),3),'g')
% 


numgrids = 25;
transmat = cell(numgrids,numgrids);
R = 9 ;
n  = size(knnemb(:,1)) ;

% - Build grid.
nBinsX = numgrids;
nBinsY = numgrids;
xg     = linspace( -R, R, nBinsX+1 ) ;
yg     = linspace( -R, R, nBinsY+1 ) ;
nCells = nBinsX * nBinsY ;

% - Build set of unique IDs for cells.
xId = sum( bsxfun( @ge, knnemb(:,1), xg(1:end-1) ), 2 ) ;
yId = sum( bsxfun( @ge, knnemb(:,2), yg(1:end-1) ), 2 ) ;
cellId = nBinsY * (xId - 1) + yId ;

[X,Y]  = meshgrid( (xg(1:end-1)+xg(2:end))/2, (yg(1:end-1)+yg(2:end))/2 ) ;

for j = 1:max(medlist_test(:,237))
    maxcell = max(medlist_test(find(medlist_test(:,237) == j),239));
    for i = 1:maxcell
        %tempindex = tracks(find(tracks(:,239)==i),:);
        %tempmovmean = movmean(tempindex(:,3) - min(tempindex(:,3)),3);
        tempx = knnemb(find(medlist_test(:,239)==i & medlist_test(:,237) == j),1);
        tempy = knnemb(find(medlist_test(:,239)==i & medlist_test(:,237) == j),2);
        temptime = medlist_test(find(medlist_test(:,239)==i & medlist_test(:,237) == j),238);
        temptimediff = diff(temptime);
        for k = 1:size(temptimediff)
            if temptimediff(k) == 1
                transmat{yId(find(medlist_test(:,239)==i & medlist_test(:,237) == j & medlist_test(:,238)==temptime(k)))+1,xId(find(medlist_test(:,239)==i & medlist_test(:,237) == j & medlist_test(:,238)==temptime(k)))+1} = [transmat{yId(find(medlist_test(:,239)==i & medlist_test(:,237) == j & medlist_test(:,238)==temptime(k)))+1,xId(find(medlist_test(:,239)==i & medlist_test(:,237) == j & medlist_test(:,238)==temptime(k)))+1}; [tempx(k+1)-tempx(k),tempy(k+1)-tempy(k)]];
            end
        end
    end
end

transmat_x = zeros(numgrids);
transmat_y = zeros(numgrids);
transmatpool = [];

for i = 1:numgrids
    for j = 1:numgrids
        if ~isempty(transmat{i,j})
            if size(transmat{i,j},1) > 7
                transmat_x(i,j) = nanmean(transmat{i,j}(:,1));
                transmat_y(i,j) = nanmean(transmat{i,j}(:,2));
                transmatpool = [transmatpool;sqrt(transmat{i,j}(:,1).^2+transmat{i,j}(:,2).^2)];

            end
        end
    end
end


figure()
hold on
f1 = scatter(embedding(:,1),embedding(:,2),5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
quiver(X-.5,Y-.5,transmat_y,transmat_x,'LineWidth', 1.5,'AutoScaleFactor',1.25,'Color','r','MaxHeadSize', .15)
axis([-8 7 -7 7])
axis off


%% Autocorrelation of prediction

[labels2,score2] = predict(trainedModel1,horzcat(medlist_test(:,3:236),2*ones(9942,1)));

lagshiftmat = NaN(size(medlist_test(:,1),1),10);
for j = 1:max(medlist_test(:,237))
    maxcell = max(medlist_test(find(medlist_test(:,237) == j),239));
    for i = 1:maxcell
        %tempindex = tracks(find(tracks(:,239)==i),:);
        %tempmovmean = movmean(tempindex(:,3) - min(tempindex(:,3)),3);
        tempcelltime = medlist_test(find(medlist_test(:,239)==i & medlist_test(:,237) == j),238);
        tempcellscore = score2(find(medlist_test(:,239)==i & medlist_test(:,237) == j),2);
        tempcelllabel = labels2(find(medlist_test(:,239)==i & medlist_test(:,237) == j));

        lagshiftidx = find(medlist_test(:,239)==i & medlist_test(:,237) == j);

        for k = 2:size(tempcelltime)
            templagshift = tempcelltime - tempcelltime(k);
            for l = 1:size(tempcelltime)
                if templagshift(l) > 0
                    lagshiftmat(lagshiftidx(l),templagshift(l)) = tempcellscore(k);
                end
                                
            end
        end
    end
end


lagshiftmat = horzcat(score2(:,2),lagshiftmat);

autocorrmat = []
for i = 1:8
    [autocorr_r, autocorr_p] = corr(lagshiftmat(:,1),lagshiftmat(:,i), 'rows','pairwise');
    autocorrmat = [autocorrmat; autocorr_r, autocorr_p]
end


[labelsa,scorea] = predict(trainedModel1,horzcat(cell1(1:7,3:236),2*ones(7,1)));
[labelsb,scoreb] = predict(trainedModel1,horzcat(cell2(1:7,3:236),1*ones(7,1)));
[labelsc,scorec] = predict(trainedModel1,horzcat(cell3(1:7,3:236),2*ones(7,1)));
[labelsd,scored] = predict(trainedModel1,horzcat(cell4(1:7,3:236),2*ones(7,1)));

figure(3)
clf
hold on
plot(-30:5:0,movmean(scorea(1:7,2),3),'r')
plot(-30:5:0,movmean(scoreb(1:7,2),3),'y')
plot(-30:5:0,movmean(scorec(1:7,2),3),'m')
plot(-30:5:0,movmean(scored(1:7,2),3),'g')

figure(1)
hold on
cdfplot(transmatpool)
prctile(transmatpool(find(~isnan(transmatpool))),73.73)



%% ikb boxplots

%rowsusedtrainedmodel = [alldosedata(trainedModel_norm1.RowsUsed,:),predictmodel(trainedModel_norm1.RowsUsed), activedummy(trainedModel_norm1.RowsUsed),trainedModel_norm1.Gradient];
%classp4 = classperf(rowsusedtrainedmodel(rowsusedtrainedmodel(239,:)<1,238),rowsusedtrainedmodel(rowsusedtrainedmodel(239,:)<1,237))


coloringlist = csvread('coloringlist.csv');
meanembedding = csvread('meanembedding.csv');


ikbgroup1_14 = ikbstain13(find(ikbstain13(:,237) == 14),:);
ikbgroup1_15 = ikbstain13(find(ikbstain13(:,237) == 15),:);
ikbgroup1_16 = ikbstain13(find(ikbstain13(:,237) == 16),:);
ikbgroup1_17 = ikbstain13(find(ikbstain13(:,237) == 17),:);
ikbgroup1_18 = ikbstain13(find(ikbstain13(:,237) == 18),:);



min14a = [];
for i = 1:size(ikbgroup1_14,1)
    min14a= [min14a,((ikbgroup1_14(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 1),31)).^2 + (ikbgroup1_14(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 1),32)).^2).^(1/2)];
end
[val,min14a_i] = min(min14a);
min14a_i= [min14a_i;val]';

min15a = [];
for i = 1:size(ikbgroup1_15,1)
    min15a= [min15a,((ikbgroup1_15(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 2),31)).^2 + (ikbgroup1_15(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 2),32)).^2).^(1/2)];
end
[val,min15a_i] = min(min15a);
min15a_i= [min15a_i;val]';

min16a = [];
for i = 1:size(ikbgroup1_16,1)
    min16a= [min16a,((ikbgroup1_16(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 3),31)).^2 + (ikbgroup1_16(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 3),32)).^2).^(1/2)];
end
[val,min16a_i] = min(min16a);
min16a_i= [min16a_i;val]';

min17a = [];
for i = 1:size(ikbgroup1_17,1)
    min17a= [min17a,((ikbgroup1_17(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 4),31)).^2 + (ikbgroup1_17(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 4),32)).^2).^(1/2)];
end
[val,min17a_i] = min(min17a);
min17a_i= [min17a_i;val]';

min18a = [];
for i = 1:size(ikbgroup1_18,1)
    min18a= [min18a,((ikbgroup1_18(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 5),31)).^2 + (ikbgroup1_18(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 5),32)).^2).^(1/2)];
end
[val,min18a_i] = min(min18a);
min18a_i= [min18a_i;val]';


linkage14a = ikbstaindata(find(ikbstaindata(:,61) == 1),:);
linkage14b = [];
for i = 1:size(min14a_i)
    linkage14b = [linkage14b; linkage14a(min14a_i(i,1),:)];
end
linkage14c = unique(min14a_i(:,1));

[linkage14d,~] = hist(min14a_i(:,1),1:numel(linkage14a(:,1)));
ismember(min14a_i(:,1),find(linkage14d==1))
linkage14h = linkage14b(ismember(min14a_i(:,1),find(linkage14d==1)),:);


linkage15a = ikbstaindata(find(ikbstaindata(:,61) == 2),:);
linkage15b = [];
for i = 1:size(min15a_i)
    linkage15b = [linkage15b; linkage15a(min15a_i(i,1),:)];
end
linkage15c = unique(min15a_i(:,1));
[linkage15d,~] = hist(min15a_i(:,1),1:numel(linkage15a(:,1)));
ismember(min15a_i(:,1),find(linkage15d==1));
linkage15h = linkage15b(ismember(min15a_i(:,1),find(linkage15d==1)),:);

linkage16a = ikbstaindata(find(ikbstaindata(:,61) == 3),:);
linkage16b = [];
for i = 1:size(min16a_i)
    linkage16b = [linkage16b; linkage16a(min16a_i(i,1),:)];
end
linkage16c = unique(min16a_i(:,1));
[linkage16d,~] = hist(min16a_i(:,1),1:numel(linkage16a(:,1)));
ismember(min16a_i(:,1),find(linkage16d==1))
linkage16h = linkage16b(ismember(min16a_i(:,1),find(linkage16d==1)),:);

linkage17a = ikbstaindata(find(ikbstaindata(:,61) == 4),:);
linkage17b = [];
for i = 1:size(min17a_i)
    linkage17b = [linkage17b; linkage17a(min17a_i(i,1),:)];
end
linkage17c = unique(min17a_i(:,1));
[linkage17d,~] = hist(min17a_i(:,1),1:numel(linkage17a(:,1)));
ismember(min17a_i(:,1),find(linkage17d==1))
linkage17h = linkage17b(ismember(min17a_i(:,1),find(linkage17d==1)),:);

linkage18a = ikbstaindata(find(ikbstaindata(:,61) == 5),:);
linkage18b = [];
for i = 1:size(min18a_i)
    linkage18b = [linkage18b; linkage18a(min18a_i(i,1),:)];
end
linkage18c = unique(min18a_i(:,1));
[linkage18d,~] = hist(min18a_i(:,1),1:numel(linkage18a(:,1)));
ismember(min18a_i(:,1),find(linkage18d==1))
linkage18h = linkage18b(ismember(min18a_i(:,1),find(linkage18d==1)),:);


linkage = [linkage14h;linkage15h;linkage16h;linkage17h;linkage18h];
linkagedummy = [ismember(min14a_i(:,1),find(linkage14d==1));ismember(min15a_i(:,1),find(linkage15d==1));ismember(min16a_i(:,1),find(linkage16d==1));ismember(min17a_i(:,1),find(linkage17d==1));ismember(min18a_i(:,1),find(linkage18d==1))];

ikbstain13_rem = ikbstain13(ikbstain13(:,237)~=13,:);
score4 = score1(ikbstain13(:,237)~=13,2);
score4 = score4(linkagedummy)


score5 = score4>1
score5 = score5*2
score5 = score5 + double(score4<-1)
boxplot(linkage(:,58),score5)




ikbgroup2_14 = ikbgroup2(find(ikbgroup2(:,237) == 14),:);
ikbgroup2_15 = ikbgroup2(find(ikbgroup2(:,237) == 15),:);
ikbgroup2_16 = ikbgroup2(find(ikbgroup2(:,237) == 16),:);
ikbgroup2_17 = ikbgroup2(find(ikbgroup2(:,237) == 17),:);
ikbgroup2_18 = ikbgroup2(find(ikbgroup2(:,237) == 18),:);

max14a = [];
for i = 1:size(ikbgroup2_14,1)
    max14a= [max14a,((ikbgroup2_14(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 1),31)).^2 + (ikbgroup2_14(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 1),32)).^2).^(1/2)];
end
[val,max14a_i] = min(max14a);
max14a_i= [max14a_i;val]';

max15a = [];
for i = 1:size(ikbgroup2_15,1)
    max15a= [min15a,((ikbgroup2_15(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 2),31)).^2 + (ikbgroup2_15(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 2),32)).^2).^(1/2)];
end
[val,max15a_i] = min(max15a);
max15a_i= [max15a_i;val]';

max16a = [];
for i = 1:size(ikbgroup2_16,1)
    max16a= [max16a,((ikbgroup2_16(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 3),31)).^2 + (ikbgroup2_16(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 3),32)).^2).^(1/2)];
end
[val,max16a_i] = min(max16a);
max16a_i= [max16a_i;val]';

max17a = [];
for i = 1:size(ikbgroup2_17,1)
    max17a= [max17a,((ikbgroup2_17(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 4),31)).^2 + (ikbgroup2_17(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 4),32)).^2).^(1/2)];
end
[val,max17a_i] = min(max17a);
max17a_i= [max17a_i;val]';

max18a = [];
for i = 1:size(ikbgroup2_18,1)
    max18a= [max18a,((ikbgroup2_18(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 5),31)).^2 + (ikbgroup2_18(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 5),32)).^2).^(1/2)];
end
[val,max18a_i] = min(max18a);
max18a_i= [max18a_i;val]';





% A14 = [A14,zeros(size(A14,1),1)];
% A15 = [A15,zeros(size(A15,1),1)];
% A16 = [A16,zeros(size(A16,1),1)];
% A17 = [A17,zeros(size(A17,1),1)];
% A18 = [A18,zeros(size(A18,1),1)];

A14(:,61) = zeros(size(A14,1),1);
A15(:,61) = zeros(size(A15,1),1);
A16(:,61) = zeros(size(A16,1),1);
A17(:,61) = zeros(size(A17,1),1);
A18(:,61) = zeros(size(A18,1),1);


for i = 1:size(min14a_i)
    if min14a_i(i,2)<200
        A14(min14a_i(i,1),61) = 1;
    end
end
for i = 1:size(max14a_i)
    if max14a_i(i,2)<200
       A14(max14a_i(i),61) = 2;
    end
end

for i = 1:size(min15a_i)
    if min15a_i(i,2)<200
        A15(min15a_i(i),61) = 1;
    end
end
for i = 1:size(max15a_i)
    if max15a_i(i,2)<200
        A15(max15a_i(i),61) = 2;
    end
end

for i = 1:size(min16a_i)
    if min16a_i(i,2)<200
        A16(min16a_i(i),61) = 1;
    end
end

for i = 1:size(max16a_i)
    if max16a_i(i,2)<200
        A16(max16a_i(i),61) = 2;
    end
end

for i = 1:size(min17a_i)
    if min17a_i(i,2)<200
        A17(min17a_i(i),61) = 1;
    end
end
for i = 1:size(max17a_i)
    if max17a_i(i,2)<200
       A17(max17a_i(i),61) = 2;
    end
end

for i = 1:size(min18a_i)
    if min18a_i(i,2)<200
        A18(min18a_i(i),61) = 1;
    end
end
for i = 1:size(max18a_i)
    if max18a_i(i,2)<200
        A18(max18a_i(i),61) = 2;
    end
end

ikbstaindata_1 = [A14;A15;A16;A17;A18];
ikbstaindata_2 = [A14;A15;A16;A17;A18];


%boxplots for figure 5
%ikbtot,p65nuc,ratio,p65nucfrac,area
%49,52,56,60,2,
testvars = [52,56,49,2];
titlevars = {'p65nuc','ikb:nfkb ratio','total ikb','cell area'};
for i = 1:4
    x1 = ikbstaindata_1(find(ikbstaindata_1(:,61) ==1),testvars(i));
    x2 = ikbstaindata_1(find(ikbstaindata_1(:,61) ==2),testvars(i));
%    x3 = ikbstaindata_2(find(ikbstaindata_2(:,61) ==1),testvars(i));
%    x4 = ikbstaindata_2(ikbstaindata_2(:,61) ==2,testvars(i));

    x = [x1; x2]%; x3; x4];


    g = [zeros(length(x1), 1); ones(length(x2), 1)];%0*ones(length(x3), 1);1*ones(length(x4), 1)];
    figure()    
    boxplot(x, g)
    title(titlevars{i})
end

%% Classifier
trainingdata = alldosedataclean(:,3:236);
function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236', 'column_237'});

predictorNames = {'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_237;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationSVM = fitcsvm(...
    predictors, ...
    response, ...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 3, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'ClassNames', [0; 1]);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2018a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 234 columns because this model was trained using 234 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236', 'column_237'});

predictorNames = {'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_237;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation predictions
[~, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

end


%% Coupling

load("MatlabWorkspaceRevision.mat")


%Predict and compare Different Cell behaviors

traceclustermat = NaN(5000,30);


[max(med01(:,239)),max(med02(:,239)),max(med03(:,239)),max(med04(:,239)),max(med05(:,239)),max(med06(:,239)),max(med07(:,239)),max(med08(:,239)),max(med09(:,239)),max(med10(:,239)),max(med11(:,239)),max(med12(:,239)),max(med13(:,239)),max(med14(:,239)),max(med15(:,239))]

meddataclean1 = meddataclean;
meddataclean1(99:149,239)=meddataclean1(99:149,239)+111;
meddataclean1(150:237,239)=meddataclean1(150:237,239)+111 +59;
meddataclean1(238:304,239)=meddataclean1(238:304,239)+111+59+105;
meddataclean1(305:371,239)=meddataclean1(305:371,239)+111+59+105+73;
meddataclean1(372:462,239)=meddataclean1(372:462,239)+111+59+105+73+89;
meddataclean1(463:554,239)=meddataclean1(463:554,239)+111+59+105+73+89+109;
meddataclean1(555:648,239)=meddataclean1(555:648,239)+111+59+105+73+89+109+100;
meddataclean1(649:708,239)=meddataclean1(649:708,239)+111+59+105+73+89+109+100+120;
meddataclean1(709:765,239)=meddataclean1(709:765,239)+111+59+105+73+89+109+100+120+77;
meddataclean1(766:826,239)=meddataclean1(766:826,239)+111+59+105+73+89+109+100+120+77+64;
meddataclean1(827:884,239)=meddataclean1(827:884,239)+111+59+105+73+89+109+100+120+77+64+70;
meddataclean1(885:958,239)=meddataclean1(885:958,239)+111+59+105+73+89+109+100+120+77+64+70+69;
meddataclean1(959:1048,239)=meddataclean1(959:1048,239)+111+59+105+73+89+109+100+120+77+64+70+69+84;
meddataclean1(1049:1109,239)=meddataclean1(1049:1109,239)+111+59+105+73+89+109+100+120+77+64+70+69+84+105;

meddataclean1 = horzcat(meddataclean1,medclusterid(meddataclean1(:,239),2));

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==1))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==2))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==3))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==4))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==5))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==6))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);


medmltest = meddataclean1(:,3:236);
medmltest1 = horzcat(medmltest,meddataclean1(:,243));
medmltest1_2 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 2),:)
medmltest1_3 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 3),:)
medmltest1_4 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 4),:)
medmltest1_5 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 5),:)
medmltest1_6 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 6),:)

medmltest2 = horzcat(medmltest,meddataclean1(:,244));
medmltest2_3 = medmltest2(find(meddataclean1(:,242) == 2 | meddataclean1(:,242) == 3),:)
medmltest2_4 = medmltest2(find(meddataclean1(:,242) == 2 | meddataclean1(:,242) == 4),:)
medmltest2_5 = medmltest2(find(meddataclean1(:,242) == 2 | meddataclean1(:,242) == 5),:)
medmltest2_6 = medmltest2(find(meddataclean1(:,242) == 2 | meddataclean1(:,242) == 6),:)

medmltest3 = horzcat(medmltest,meddataclean1(:,245));
medmltest3_4 = medmltest3(find(meddataclean1(:,242) == 3 | meddataclean1(:,242) == 4),:)
medmltest3_5 = medmltest3(find(meddataclean1(:,242) == 3 | meddataclean1(:,242) == 5),:)
medmltest3_6 = medmltest3(find(meddataclean1(:,242) == 3 | meddataclean1(:,242) == 6),:)

medmltest4 = horzcat(medmltest,meddataclean1(:,246));
medmltest4_5 = medmltest4(find(meddataclean1(:,242) == 4 | meddataclean1(:,242) == 5),:)
medmltest4_6 = medmltest4(find(meddataclean1(:,242) == 4 | meddataclean1(:,242) == 6),:)

medmltest5 = horzcat(medmltest,meddataclean1(:,247));
medmltest5_6 = medmltest5(find(meddataclean1(:,242) == 5 | meddataclean1(:,242) == 6),:)

medmltest6 = horzcat(medmltest,meddataclean1(:,248));


figure()
hold on
plot(nanmean(traceclustermat(find(medclusterid(:,2)==1),:)),'r')
plot(nanmean(traceclustermat(find(medclusterid(:,2)==2),:)),'k')
plot(nanmean(traceclustermat(find(medclusterid(:,2)==3),:)),'b')
plot(nanmean(traceclustermat(find(medclusterid(:,2)==4),:)),'g')
plot(nanmean(traceclusterinactive),'y')
plot(nanmean(traceclustermat(find(medclusterid(:,2)==6),:)),'m')



%show variables associated with trace differences
boxplot(meddataclean1(:,11),meddataclean1(:,242),'OutlierSize',3,'Symbol', 'ko')

boxplot(meddataclean1(:,94),meddataclean1(:,242),'OutlierSize',3,'Symbol', 'ko')
   
boxplot(meddataclean1(:,12),meddataclean1(:,242),'OutlierSize',3,'Symbol', 'ko')


% Corr plot
corrmatrix3 = [];
for i = 1:235
    corrmatrix3 = [corrmatrix3,corr(medmltest3_4(:,i),medmltest3_4(:,235))];
end
[~,indexcorr3] = sort(corrmatrix3);

corrsorted3 = corrmatrix3(indexcorr3);

figure(); hold on
barh(1:20,corrsorted3(1:20),'FaceColor',[.2,0.2,0.2],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)
barh(21:41,corrsorted3(207:227),'FaceColor',[1,0,0],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)


%Use Classification ensemble to visualize coupling

view(trainedModel.ClassificationEnsemble.Trained{8},'Mode','graph');


%find accuracy through repeated thresholded
figure()
box(alldosedataclean(:,9),activedummy+1)

couplingactive = alldosedataclean(:,9) > 3625.57;
classp = classperf(activedummy,couplingactive)

couplingtrue = activedummy(find(alldosedataclean(:,8)<1.58e06),:);
couplingactive = alldosedataclean(find(alldosedataclean(:,8)<1.58e06),:);
couplingactive = couplingactive(:,9) > 3625.57;
classp1 = classperf(couplingtrue,couplingactive)

couplingtrue = activedummy(find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741),:);
couplingactive = alldosedataclean(find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741),:);
couplingactive = couplingactive(:,9) > 3625.57;
classp2 = classperf(couplingtrue,couplingactive)

couplingtrue = activedummy(find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741 & alldosedataclean(:,4)<114.12),:);
couplingactive = alldosedataclean(find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741 & alldosedataclean(:,4)<114.12),:);
couplingactive = couplingactive(:,9) > 3625.57;

couplingdummy = zeros(3456,1);
couplingdummy1 = find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741 & alldosedataclean(:,4)<114.12 & alldosedataclean(:,9) < 3625.57);
for i = 1:size(couplingdummy1)
    couplingdummy(couplingdummy1(i)) = 1;
end

%Coupled population corr chart

corrmatrix1 = [];
for i = 1:236
    corrmatrix1 = [corrmatrix1,corr(alldosedataclean(:,i),couplingdummy)];
end
[~,indexcorr1] = sort(corrmatrix1);

corrsorted1 = corrmatrix1(indexcorr1);

figure(); hold on
barh(1:20,corrsorted1(1:20),'FaceColor',[.2,0.2,0.2],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)
barh(21:41,corrsorted1(208:228),'FaceColor',[1,0,0],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)

%Visualize coupled population
   
figure()
hold on
plot(embedding(couplingdummy==1,1),embedding(couplingdummy==1,2),'r' ,'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(couplingdummy==0,1),embedding(couplingdummy==0,2), 'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('Coupling')




corrmatrix2 = [];
for i = 1:236
    corrmatrix2 = [corrmatrix2,corr(alldosedataclean(:,i),alldosedataclean(:,3))];
end
[~,indexcorr2] = sort(corrmatrix2);

corrsorted2 = corrmatrix2(indexcorr2);

figure(); hold on
barh(1:20,corrsorted2(1:20),'FaceColor',[.2,0.2,0.2],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)
barh(21:41,corrsorted2(208:228),'FaceColor',[1,0,0],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)


classperf(couplingtrue,couplingactive)


%% test different threshholds for activation
activedummy = zeros(size(alldosedataclean,1),1);
activedummy(find(alldosedataclean(:,242)> 100)) = 1;

figure()
hold on
scatter(embedding(find(activedummy==1),1),embedding(find(activedummy==1),2), 6, [1 0 0],'filled');
scatter(embedding(find(activedummy==0),1),embedding(find(activedummy==0),2), 4, [.5 .5 .5],'filled');
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('prediction')

activedummy = zeros(size(alldosedataclean,1),1);
activedummy(find(alldosedataclean(:,242)> 300)) = 1;

figure()
hold on
scatter(embedding(find(activedummy==1),1),embedding(find(activedummy==1),2), 6, [1 0 0],'filled');
scatter(embedding(find(activedummy==0),1),embedding(find(activedummy==0),2), 4, [.5 .5 .5],'filled');
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('prediction')


activedummy = zeros(size(alldosedataclean,1),1);
activedummy(find(alldosedataclean(:,242)> 500)) = 1;

figure()
hold on
scatter(embedding(find(activedummy==1),1),embedding(find(activedummy==1),2), 6, [1 0 0],'filled');
scatter(embedding(find(activedummy==0),1),embedding(find(activedummy==0),2), 4, [.5 .5 .5],'filled');
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('prediction')


activedummy = zeros(size(alldosedataclean,1),1);
activedummy(find(alldosedataclean(:,242)> 700)) = 1;

figure()
hold on
scatter(embedding(find(activedummy==1),1),embedding(find(activedummy==1),2), 6, [1 0 0],'filled');
scatter(embedding(find(activedummy==0),1),embedding(find(activedummy==0),2), 4, [.5 .5 .5],'filled');
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('prediction')

   


<<<<<<< HEAD
=======
=======
%% Use track program to track centroids
% param=struct(); %structure array
% param.mem = 2; %can leave for 10 frames
% param.dim = 2;
% param.good = 10; %set high when code is working (max is number of frames)
% param.quiet = 1;
% maxdisplacement = 100;

%example: Use ImageAnalysis function as a function of XYposition, then track using track function
%xy01 = ImageAnalysis(1);
%med01 = track(xy01,maxdisplacement,param);
%graphs(med01)

%% Combine dosages and isolate t=0 for each stimulation

% 
% maxlist = []
% lowlist = vertcat(low01,low02,low03,low04,low05,low06,low07,low08,low09,low10,low11,low12,low13,low14,low15);
% lowdata7 = lowlist(find(lowlist(:,238)==7),:);
% lowdata21 = lowlist(find(lowlist(:,238)==21),:);
% maxlist = [graphs(low01);graphs(low02);graphs(low03);graphs(low04);graphs(low05);graphs(low06);graphs(low07);graphs(low08);graphs(low09);graphs(low10);graphs(low11);graphs(low12);graphs(low13);graphs(low14);graphs(low15)];
% lowdataclean = horzcat(lowdata7,maxlist)
% lowdataclean = lowdataclean(find(lowdataclean(:,240)>0),:);
% 
% 
% maxlist = []
% medlist = vertcat(med01,med02,med03,med04,med05,med06,med07,med08,med09,med10,med11,med12,med13,med14,med15);
% meddata7 = medlist(find(medlist(:,238)==7),:);
% meddata21 = medlist(find(medlist(:,238)==21),:);
% maxlist = [graphs(med01);graphs(med02);graphs(med03);graphs(med04);graphs(med05);graphs(med06);graphs(med07);graphs(med08);graphs(med09);graphs(med10);graphs(med11);graphs(med12);graphs(med13);graphs(med14);graphs(med15)];
% meddataclean = horzcat(meddata7,maxlist)
% meddataclean = meddataclean(find(meddataclean(:,240)>0),:);
% 
% 
% maxlist = []
% highlist = vertcat(high01,high02,high03,high04,high05,high06,high07,high08,high09,high10,high11,high12,high13,high15);
% highdata7 = highlist(find(highlist(:,238)==7),:);
% highdata21 = highlist(find(highlist(:,238)==21),:);
% maxlist = [graphs(high01);graphs(high02);graphs(high03);graphs(high04);graphs(high05);graphs(high06);graphs(high07);graphs(high08);graphs(high09);graphs(high10);graphs(high11);graphs(high12);graphs(high13);graphs(high15)];
% highdataclean = horzcat(highdata7,maxlist)
% highdataclean = highdataclean(find(highdataclean(:,240)>0),:);
% 
% maxlist = []
% vhlist = vertcat(vh01,vh02,vh03,vh04,vh05,vh06,vh07,vh08,vh09,vh10,vh11);
% vhdata7 = vhlist(find(vhlist(:,238)==7),:);
% vhdata21 = vhlist(find(vhlist(:,238)==21),:);
% maxlist = [graphs(vh01);graphs(vh02);graphs(vh03);graphs(vh04);graphs(vh05);graphs(vh06);graphs(vh07);graphs(vh08);graphs(vh09);graphs(vh10);graphs(vh11)];
% vhdataclean = horzcat(vhdata7,maxlist)
% vhdataclean = vhdataclean(find(vhdataclean(:,240)>0),:);
% 
% lowdataclean = horzcat(lowdataclean,ones(size(lowdataclean,1),1));
% meddataclean = horzcat(meddataclean,2*ones(size(meddataclean,1),1));
% highdataclean = horzcat(highdataclean,3*ones(size(highdataclean,1),1));
% vhdataclean = horzcat(vhdataclean,4*ones(size(vhdataclean,1),1));

load("MatlabWorkspaceBase.mat")

%Include ikb stain data
ikbstain = vertcat(xy13,xy14,xy15,xy16,xy17,xy18);
ikbstain13 = ikbstain(find(ikbstain(:,238)==13),:);

%Normalize
mmlow = minmax(lowdataclean');
mmmed = minmax(meddataclean');
mmhigh = minmax(highdataclean');
mmvh = minmax(vhdataclean');
mmikb = minmax(ikbstain13');


lowdatanorm = [];
meddatanorm = [];
highdatanorm = [];
vhdatanorm = [];
ikbnorm = [];
for i = 1:241
    if i <237
        lowdatanorm = [lowdatanorm,(lowdataclean(:,i)-mmlow(i,1))/mmlow(i,2)];
        meddatanorm = [meddatanorm,(meddataclean(:,i)-mmmed(i,1))/mmmed(i,2)];
        highdatanorm = [highdatanorm,(highdataclean(:,i)-mmhigh(i,1))/mmhigh(i,2)];
        vhdatanorm = [vhdatanorm,(vhdataclean(:,i)-mmvh(i,1))/mmvh(i,2)];
        ikbnorm = [ikbnorm, (ikbstain13(:,i)-mmikb(i,1))/mmikb(i,2)];
    end
end


alldosedataclean = vertcat(lowdataclean,meddataclean,highdataclean,vhdataclean);
alldosedata = vertcat(lowdatanorm,meddatanorm,highdatanorm,vhdatanorm);
alldosedataclean = [alldosedataclean, alldosedataclean(:,240)-alldosedataclean(:,3)];
alldosedataclean_active = alldosedataclean(find(alldosedataclean(:,242)> 500),:);
activedummy = zeros(size(alldosedataclean,1),1);
activedummy(find(alldosedataclean(:,242)> 500)) = 1;



%% tsne
rng(1)
tsneplotdose = tsne(alldosedata(:,3:15),'Algorithm','barneshut');
figure()
gscatter(tsneplotdose(:,1),tsneplotdose(:,2),alldosedataclean(:,241))
alldosedata_ikb = vertcat(lowdatanorm,meddatanorm,highdatanorm,vhdatanorm,ikbnorm);
alldosedata_norm = vertcat(lowdatanorm,meddatanorm,highdatanorm,vhdatanorm);
alldosedata_norm = horzcat(alldosedata_norm, alldosedataclean(:,241), activedummy);


%rng1
rng(1)
tsneplotdose_ikb_1 = tsne(alldosedata_ikb(:,3:15),'Algorithm','barneshut','perplexity',800,'Exaggeration', 1);

x0=10;
y0=10;

figure()
scatter(tsneplotdose_ikb_1(1:3456,1),tsneplotdose_ikb_1(1:3456,2), 3, log(alldosedataclean(:,240)));
colormap(flipud(hot))
caxis([6 9])
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
hold on

predictmodel = [];
for i = 1:3456
    predictmodel = [predictmodel; trainedModel.predictFcn(horzcat(alldosedataclean(i,3:236),alldosedataclean(i,241)))];
end


%% Embedding with UMAP

figure()
hold on
scatter(embedding(find(predictmodel==1),1),embedding(find(predictmodel==1),2), 4, [1 0 0],'filled');
scatter(embedding(find(predictmodel==0),1),embedding(find(predictmodel==0),2), 4, [.5 .5 .5],'filled');
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('prediction')


% correlation in nucleus texture feature for resistent group
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,alldosedataclean(:,22),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('correlation in nucleus texture feature for resistent group')


% consistency in nucleus (entropy)
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,alldosedataclean(:,23),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('consistency in nucleus (entropy)')


%Otsu dim
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,log(alldosedataclean(:,95)),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('Otsu dim')

%INFORMATION MEASURE OF CORR
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,alldosedataclean(:,29),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('Information measure of correlation')

%leaky
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,log(alldosedataclean(:,3)),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('leaky nfkb')

%leaky std
figure()
scatter(embedding(1:3456,1),embedding(1:3456,2), 3,log(alldosedataclean(:,4)),'filled');
colormap(flipud(hot))
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('leaky nfkb std in nucleus')


%% Embedding breakdown by group
figure()
hold on
plot(embedding(find(alldosedataclean(:,241) == 1),1),embedding(find(alldosedataclean(:,241) == 1),2),'r','LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 2),1),embedding(find(alldosedataclean(:,241) == 2),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 3),1),embedding(find(alldosedataclean(:,241) == 3),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 4),1),embedding(find(alldosedataclean(:,241) == 4),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])
yticks([-50 0 50])


figure()
hold on
plot(embedding(find(alldosedataclean(:,241) == 1),1),embedding(find(alldosedataclean(:,241) == 1),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 2),1),embedding(find(alldosedataclean(:,241) == 2),2),'r','LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 3),1),embedding(find(alldosedataclean(:,241) == 3),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 4),1),embedding(find(alldosedataclean(:,241) == 4),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])
yticks([-50 0 50])



figure()
hold on
plot(embedding(find(alldosedataclean(:,241) == 1),1),embedding(find(alldosedataclean(:,241) == 1),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 2),1),embedding(find(alldosedataclean(:,241) == 2),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 3),1),embedding(find(alldosedataclean(:,241) == 3),2),'r','LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 4),1),embedding(find(alldosedataclean(:,241) == 4),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])
yticks([-50 0 50])

figure()
hold on
plot(embedding(find(alldosedataclean(:,241) == 1),1),embedding(find(alldosedataclean(:,241) == 1),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 2),1),embedding(find(alldosedataclean(:,241) == 2),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 3),1),embedding(find(alldosedataclean(:,241) == 3),2),'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(find(alldosedataclean(:,241) == 4),1),embedding(find(alldosedataclean(:,241) == 4),2),'Color','r','LineStyle','none','Marker','o','MarkerSize', 1)
x0=10;
y0=10;
width=200;
height=200;
set(gcf,'position',[x0,y0,width,height])
yticks([-50 0 50])

%% Finding high probability features


pvaluestest1 = [];
logfoldchangetest1 = [];
for i = 1:238
    [~,p] = ttest2(alldosedataclean(find(predictmodel(:)==0),i),alldosedataclean(find(predictmodel(:)==1),i));
    logfoldchange = mean(alldosedataclean(find(predictmodel(:)==0),i))/mean(alldosedataclean(find(predictmodel(:)==1),i));
    pvaluestest1 = [pvaluestest1,p];
    logfoldchangetest1 = [logfoldchangetest1,logfoldchange];
end
figure(); hold on
scatter(-log(logfoldchangetest1),log10(-log10(pvaluestest1)),3,'k')
scatter(-log(logfoldchangetest1([99:107])),log10(-log10(pvaluestest1([99:107]))),25,'r')


pvaluessig = find(log10(-log10(pvaluestest1))>1);
[~,indexpvalues] = sort(-log(logfoldchangetest1));

corrmatrix = [];
for i = 1:238
    corrmatrix = [corrmatrix,corr(alldosedataclean(:,i),predictmodel)];
end
[~,indexcorr] = sort(corrmatrix);


corrsorted = corrmatrix(indexcorr);

figure(); hold on
barh(1:20,corrsorted(1:20),'FaceColor',[.2,0.2,0.2],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)
barh(21:41,corrsorted(211:231),'FaceColor',[1,0,0],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)



%% Morphology velocity
medlist_test = medlist(find(medlist(:,238)<10),:);
lowlist_test = lowlist(find(lowlist(:,238)<10),:);

[knnidx,knnweight] = knnsearch(alldosedataclean(:,3:15),medlist_test(:,3:15),'K', 2);
knnemb = [];
for i = 1:size(knnidx,1)
    tempmat = (1./knnweight(i,:))./(sum(1./knnweight(i,:)));
    knnemb = [knnemb; mean(2*tempmat'.*embedding(knnidx(i,:),:))];
end

%figure()
%scatter(embedding(:,1),embedding(:,2),3)
%hold on; scatter(knnemb(:,1),knnemb(:,2),2)

figure()
f1 = scatter(embedding(:,1),embedding(:,2),5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
hold on
countval =0;
for j = 1:max(medlist_test(:,237))
    maxcell = max(medlist_test(find(medlist_test(:,237) == j),239));
    for i = 1:maxcell
        %tempindex = tracks(find(tracks(:,239)==i),:);
        %tempmovmean = movmean(tempindex(:,3) - min(tempindex(:,3)),3);
        plot(knnemb(find(medlist_test(:,239)==i & medlist_test(:,237) == j),1),knnemb(find(medlist_test(:,239)==i & medlist_test(:,237) == j),2))
    end
end

axis([-8 7 -7 7])
axis off

% cell3=medlist_test(find(medlist_test(:,239)==63 & medlist_test(:,237) == 2),:);
% cell4=medlist_test(find(medlist_test(:,239)==16 & medlist_test(:,237) == 4),:);
% cell1=lowlist_test(find(lowlist_test(:,239)==9 & lowlist_test(:,237) == 7),:);
% cell2=medlist_test(find(medlist_test(:,239)==19 & medlist_test(:,237) == 15),:);
% 
% figure()
% clf
% hold on
% plot(-30:5:0,movmean(cell1(1:7,3),3),'r')
% plot(-30:5:0,movmean(cell2(1:7,3),3),'y')
% plot(-30:5:0,movmean(cell3(1:7,3),3),'m')
% plot(-30:5:0,movmean(cell4(1:7,3),3),'g')
% 
% %pixel entropy
% figure()
% clf
% hold on
% plot(-30:5:0,movmean(cell1(1:7,27),3),'r')
% plot(-30:5:0,movmean(cell2(1:7,27),3),'y')
% plot(-30:5:0,movmean(cell3(1:7,27),3),'m')
% plot(-30:5:0,movmean(cell4(1:7,27),3),'g')
% 
% 
% %pixel correlation
% figure()
% clf
% hold on
% plot(-30:5:0,movmean(cell1(1:7,22),3),'r')
% plot(-30:5:0,movmean(cell2(1:7,22),3),'y')
% plot(-30:5:0,movmean(cell3(1:7,22),3),'m')
% plot(-30:5:0,movmean(cell4(1:7,22),3),'g')
% 
% %information measure
% figure()
% clf
% hold on
% plot(-30:5:0,movmean(cell1(1:7,67),3),'r')
% plot(-30:5:0,movmean(cell2(1:7,67),3),'y')
% plot(-30:5:0,movmean(cell3(1:7,67),3),'m')
% plot(-30:5:0,movmean(cell4(1:7,67),3),'g')
% 


numgrids = 25;
transmat = cell(numgrids,numgrids);
R = 9 ;
n  = size(knnemb(:,1)) ;

% - Build grid.
nBinsX = numgrids;
nBinsY = numgrids;
xg     = linspace( -R, R, nBinsX+1 ) ;
yg     = linspace( -R, R, nBinsY+1 ) ;
nCells = nBinsX * nBinsY ;

% - Build set of unique IDs for cells.
xId = sum( bsxfun( @ge, knnemb(:,1), xg(1:end-1) ), 2 ) ;
yId = sum( bsxfun( @ge, knnemb(:,2), yg(1:end-1) ), 2 ) ;
cellId = nBinsY * (xId - 1) + yId ;

[X,Y]  = meshgrid( (xg(1:end-1)+xg(2:end))/2, (yg(1:end-1)+yg(2:end))/2 ) ;

for j = 1:max(medlist_test(:,237))
    maxcell = max(medlist_test(find(medlist_test(:,237) == j),239));
    for i = 1:maxcell
        %tempindex = tracks(find(tracks(:,239)==i),:);
        %tempmovmean = movmean(tempindex(:,3) - min(tempindex(:,3)),3);
        tempx = knnemb(find(medlist_test(:,239)==i & medlist_test(:,237) == j),1);
        tempy = knnemb(find(medlist_test(:,239)==i & medlist_test(:,237) == j),2);
        temptime = medlist_test(find(medlist_test(:,239)==i & medlist_test(:,237) == j),238);
        temptimediff = diff(temptime);
        for k = 1:size(temptimediff)
            if temptimediff(k) == 1
                transmat{yId(find(medlist_test(:,239)==i & medlist_test(:,237) == j & medlist_test(:,238)==temptime(k)))+1,xId(find(medlist_test(:,239)==i & medlist_test(:,237) == j & medlist_test(:,238)==temptime(k)))+1} = [transmat{yId(find(medlist_test(:,239)==i & medlist_test(:,237) == j & medlist_test(:,238)==temptime(k)))+1,xId(find(medlist_test(:,239)==i & medlist_test(:,237) == j & medlist_test(:,238)==temptime(k)))+1}; [tempx(k+1)-tempx(k),tempy(k+1)-tempy(k)]];
            end
        end
    end
end

transmat_x = zeros(numgrids);
transmat_y = zeros(numgrids);
transmatpool = [];

for i = 1:numgrids
    for j = 1:numgrids
        if ~isempty(transmat{i,j})
            if size(transmat{i,j},1) > 7
                transmat_x(i,j) = nanmean(transmat{i,j}(:,1));
                transmat_y(i,j) = nanmean(transmat{i,j}(:,2));
                transmatpool = [transmatpool;sqrt(transmat{i,j}(:,1).^2+transmat{i,j}(:,2).^2)];

            end
        end
    end
end


figure()
hold on
f1 = scatter(embedding(:,1),embedding(:,2),5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
quiver(X-.5,Y-.5,transmat_y,transmat_x,'LineWidth', 1.5,'AutoScaleFactor',1.25,'Color','r','MaxHeadSize', .15)
axis([-8 7 -7 7])
axis off


%% Autocorrelation of prediction

[labels2,score2] = predict(trainedModel1,horzcat(medlist_test(:,3:236),2*ones(9942,1)));

lagshiftmat = NaN(size(medlist_test(:,1),1),10);
for j = 1:max(medlist_test(:,237))
    maxcell = max(medlist_test(find(medlist_test(:,237) == j),239));
    for i = 1:maxcell
        %tempindex = tracks(find(tracks(:,239)==i),:);
        %tempmovmean = movmean(tempindex(:,3) - min(tempindex(:,3)),3);
        tempcelltime = medlist_test(find(medlist_test(:,239)==i & medlist_test(:,237) == j),238);
        tempcellscore = score2(find(medlist_test(:,239)==i & medlist_test(:,237) == j),2);
        tempcelllabel = labels2(find(medlist_test(:,239)==i & medlist_test(:,237) == j));

        lagshiftidx = find(medlist_test(:,239)==i & medlist_test(:,237) == j);

        for k = 2:size(tempcelltime)
            templagshift = tempcelltime - tempcelltime(k);
            for l = 1:size(tempcelltime)
                if templagshift(l) > 0
                    lagshiftmat(lagshiftidx(l),templagshift(l)) = tempcellscore(k);
                end
                                
            end
        end
    end
end


lagshiftmat = horzcat(score2(:,2),lagshiftmat);

autocorrmat = []
for i = 1:8
    [autocorr_r, autocorr_p] = corr(lagshiftmat(:,1),lagshiftmat(:,i), 'rows','pairwise');
    autocorrmat = [autocorrmat; autocorr_r, autocorr_p]
end


[labelsa,scorea] = predict(trainedModel1,horzcat(cell1(1:7,3:236),2*ones(7,1)));
[labelsb,scoreb] = predict(trainedModel1,horzcat(cell2(1:7,3:236),1*ones(7,1)));
[labelsc,scorec] = predict(trainedModel1,horzcat(cell3(1:7,3:236),2*ones(7,1)));
[labelsd,scored] = predict(trainedModel1,horzcat(cell4(1:7,3:236),2*ones(7,1)));

figure(3)
clf
hold on
plot(-30:5:0,movmean(scorea(1:7,2),3),'r')
plot(-30:5:0,movmean(scoreb(1:7,2),3),'y')
plot(-30:5:0,movmean(scorec(1:7,2),3),'m')
plot(-30:5:0,movmean(scored(1:7,2),3),'g')

figure(1)
hold on
cdfplot(transmatpool)
prctile(transmatpool(find(~isnan(transmatpool))),73.73)



%% ikb boxplots

%rowsusedtrainedmodel = [alldosedata(trainedModel_norm1.RowsUsed,:),predictmodel(trainedModel_norm1.RowsUsed), activedummy(trainedModel_norm1.RowsUsed),trainedModel_norm1.Gradient];
%classp4 = classperf(rowsusedtrainedmodel(rowsusedtrainedmodel(239,:)<1,238),rowsusedtrainedmodel(rowsusedtrainedmodel(239,:)<1,237))


coloringlist = csvread('coloringlist.csv');
meanembedding = csvread('meanembedding.csv');


ikbgroup1_14 = ikbstain13(find(ikbstain13(:,237) == 14),:);
ikbgroup1_15 = ikbstain13(find(ikbstain13(:,237) == 15),:);
ikbgroup1_16 = ikbstain13(find(ikbstain13(:,237) == 16),:);
ikbgroup1_17 = ikbstain13(find(ikbstain13(:,237) == 17),:);
ikbgroup1_18 = ikbstain13(find(ikbstain13(:,237) == 18),:);



min14a = [];
for i = 1:size(ikbgroup1_14,1)
    min14a= [min14a,((ikbgroup1_14(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 1),31)).^2 + (ikbgroup1_14(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 1),32)).^2).^(1/2)];
end
[val,min14a_i] = min(min14a);
min14a_i= [min14a_i;val]';

min15a = [];
for i = 1:size(ikbgroup1_15,1)
    min15a= [min15a,((ikbgroup1_15(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 2),31)).^2 + (ikbgroup1_15(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 2),32)).^2).^(1/2)];
end
[val,min15a_i] = min(min15a);
min15a_i= [min15a_i;val]';

min16a = [];
for i = 1:size(ikbgroup1_16,1)
    min16a= [min16a,((ikbgroup1_16(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 3),31)).^2 + (ikbgroup1_16(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 3),32)).^2).^(1/2)];
end
[val,min16a_i] = min(min16a);
min16a_i= [min16a_i;val]';

min17a = [];
for i = 1:size(ikbgroup1_17,1)
    min17a= [min17a,((ikbgroup1_17(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 4),31)).^2 + (ikbgroup1_17(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 4),32)).^2).^(1/2)];
end
[val,min17a_i] = min(min17a);
min17a_i= [min17a_i;val]';

min18a = [];
for i = 1:size(ikbgroup1_18,1)
    min18a= [min18a,((ikbgroup1_18(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 5),31)).^2 + (ikbgroup1_18(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 5),32)).^2).^(1/2)];
end
[val,min18a_i] = min(min18a);
min18a_i= [min18a_i;val]';


linkage14a = ikbstaindata(find(ikbstaindata(:,61) == 1),:);
linkage14b = [];
for i = 1:size(min14a_i)
    linkage14b = [linkage14b; linkage14a(min14a_i(i,1),:)];
end
linkage14c = unique(min14a_i(:,1));

[linkage14d,~] = hist(min14a_i(:,1),1:numel(linkage14a(:,1)));
ismember(min14a_i(:,1),find(linkage14d==1))
linkage14h = linkage14b(ismember(min14a_i(:,1),find(linkage14d==1)),:);


linkage15a = ikbstaindata(find(ikbstaindata(:,61) == 2),:);
linkage15b = [];
for i = 1:size(min15a_i)
    linkage15b = [linkage15b; linkage15a(min15a_i(i,1),:)];
end
linkage15c = unique(min15a_i(:,1));
[linkage15d,~] = hist(min15a_i(:,1),1:numel(linkage15a(:,1)));
ismember(min15a_i(:,1),find(linkage15d==1));
linkage15h = linkage15b(ismember(min15a_i(:,1),find(linkage15d==1)),:);

linkage16a = ikbstaindata(find(ikbstaindata(:,61) == 3),:);
linkage16b = [];
for i = 1:size(min16a_i)
    linkage16b = [linkage16b; linkage16a(min16a_i(i,1),:)];
end
linkage16c = unique(min16a_i(:,1));
[linkage16d,~] = hist(min16a_i(:,1),1:numel(linkage16a(:,1)));
ismember(min16a_i(:,1),find(linkage16d==1))
linkage16h = linkage16b(ismember(min16a_i(:,1),find(linkage16d==1)),:);

linkage17a = ikbstaindata(find(ikbstaindata(:,61) == 4),:);
linkage17b = [];
for i = 1:size(min17a_i)
    linkage17b = [linkage17b; linkage17a(min17a_i(i,1),:)];
end
linkage17c = unique(min17a_i(:,1));
[linkage17d,~] = hist(min17a_i(:,1),1:numel(linkage17a(:,1)));
ismember(min17a_i(:,1),find(linkage17d==1))
linkage17h = linkage17b(ismember(min17a_i(:,1),find(linkage17d==1)),:);

linkage18a = ikbstaindata(find(ikbstaindata(:,61) == 5),:);
linkage18b = [];
for i = 1:size(min18a_i)
    linkage18b = [linkage18b; linkage18a(min18a_i(i,1),:)];
end
linkage18c = unique(min18a_i(:,1));
[linkage18d,~] = hist(min18a_i(:,1),1:numel(linkage18a(:,1)));
ismember(min18a_i(:,1),find(linkage18d==1))
linkage18h = linkage18b(ismember(min18a_i(:,1),find(linkage18d==1)),:);


linkage = [linkage14h;linkage15h;linkage16h;linkage17h;linkage18h];
linkagedummy = [ismember(min14a_i(:,1),find(linkage14d==1));ismember(min15a_i(:,1),find(linkage15d==1));ismember(min16a_i(:,1),find(linkage16d==1));ismember(min17a_i(:,1),find(linkage17d==1));ismember(min18a_i(:,1),find(linkage18d==1))];

ikbstain13_rem = ikbstain13(ikbstain13(:,237)~=13,:);
score4 = score1(ikbstain13(:,237)~=13,2);
score4 = score4(linkagedummy)


score5 = score4>1
score5 = score5*2
score5 = score5 + double(score4<-1)
boxplot(linkage(:,58),score5)




ikbgroup2_14 = ikbgroup2(find(ikbgroup2(:,237) == 14),:);
ikbgroup2_15 = ikbgroup2(find(ikbgroup2(:,237) == 15),:);
ikbgroup2_16 = ikbgroup2(find(ikbgroup2(:,237) == 16),:);
ikbgroup2_17 = ikbgroup2(find(ikbgroup2(:,237) == 17),:);
ikbgroup2_18 = ikbgroup2(find(ikbgroup2(:,237) == 18),:);

max14a = [];
for i = 1:size(ikbgroup2_14,1)
    max14a= [max14a,((ikbgroup2_14(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 1),31)).^2 + (ikbgroup2_14(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 1),32)).^2).^(1/2)];
end
[val,max14a_i] = min(max14a);
max14a_i= [max14a_i;val]';

max15a = [];
for i = 1:size(ikbgroup2_15,1)
    max15a= [min15a,((ikbgroup2_15(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 2),31)).^2 + (ikbgroup2_15(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 2),32)).^2).^(1/2)];
end
[val,max15a_i] = min(max15a);
max15a_i= [max15a_i;val]';

max16a = [];
for i = 1:size(ikbgroup2_16,1)
    max16a= [max16a,((ikbgroup2_16(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 3),31)).^2 + (ikbgroup2_16(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 3),32)).^2).^(1/2)];
end
[val,max16a_i] = min(max16a);
max16a_i= [max16a_i;val]';

max17a = [];
for i = 1:size(ikbgroup2_17,1)
    max17a= [max17a,((ikbgroup2_17(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 4),31)).^2 + (ikbgroup2_17(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 4),32)).^2).^(1/2)];
end
[val,max17a_i] = min(max17a);
max17a_i= [max17a_i;val]';

max18a = [];
for i = 1:size(ikbgroup2_18,1)
    max18a= [max18a,((ikbgroup2_18(i,1)-ikbstaindata(find(ikbstaindata(:,61) == 5),31)).^2 + (ikbgroup2_18(i,2)-ikbstaindata(find(ikbstaindata(:,61) == 5),32)).^2).^(1/2)];
end
[val,max18a_i] = min(max18a);
max18a_i= [max18a_i;val]';





% A14 = [A14,zeros(size(A14,1),1)];
% A15 = [A15,zeros(size(A15,1),1)];
% A16 = [A16,zeros(size(A16,1),1)];
% A17 = [A17,zeros(size(A17,1),1)];
% A18 = [A18,zeros(size(A18,1),1)];

A14(:,61) = zeros(size(A14,1),1);
A15(:,61) = zeros(size(A15,1),1);
A16(:,61) = zeros(size(A16,1),1);
A17(:,61) = zeros(size(A17,1),1);
A18(:,61) = zeros(size(A18,1),1);


for i = 1:size(min14a_i)
    if min14a_i(i,2)<200
        A14(min14a_i(i,1),61) = 1;
    end
end
for i = 1:size(max14a_i)
    if max14a_i(i,2)<200
       A14(max14a_i(i),61) = 2;
    end
end

for i = 1:size(min15a_i)
    if min15a_i(i,2)<200
        A15(min15a_i(i),61) = 1;
    end
end
for i = 1:size(max15a_i)
    if max15a_i(i,2)<200
        A15(max15a_i(i),61) = 2;
    end
end

for i = 1:size(min16a_i)
    if min16a_i(i,2)<200
        A16(min16a_i(i),61) = 1;
    end
end

for i = 1:size(max16a_i)
    if max16a_i(i,2)<200
        A16(max16a_i(i),61) = 2;
    end
end

for i = 1:size(min17a_i)
    if min17a_i(i,2)<200
        A17(min17a_i(i),61) = 1;
    end
end
for i = 1:size(max17a_i)
    if max17a_i(i,2)<200
       A17(max17a_i(i),61) = 2;
    end
end

for i = 1:size(min18a_i)
    if min18a_i(i,2)<200
        A18(min18a_i(i),61) = 1;
    end
end
for i = 1:size(max18a_i)
    if max18a_i(i,2)<200
        A18(max18a_i(i),61) = 2;
    end
end

ikbstaindata_1 = [A14;A15;A16;A17;A18];
ikbstaindata_2 = [A14;A15;A16;A17;A18];


%boxplots for figure 5
%ikbtot,p65nuc,ratio,p65nucfrac,area
%49,52,56,60,2,
testvars = [52,56,49,2];
titlevars = {'p65nuc','ikb:nfkb ratio','total ikb','cell area'};
for i = 1:4
    x1 = ikbstaindata_1(find(ikbstaindata_1(:,61) ==1),testvars(i));
    x2 = ikbstaindata_1(find(ikbstaindata_1(:,61) ==2),testvars(i));
%    x3 = ikbstaindata_2(find(ikbstaindata_2(:,61) ==1),testvars(i));
%    x4 = ikbstaindata_2(ikbstaindata_2(:,61) ==2,testvars(i));

    x = [x1; x2]%; x3; x4];


    g = [zeros(length(x1), 1); ones(length(x2), 1)];%0*ones(length(x3), 1);1*ones(length(x4), 1)];
    figure()    
    boxplot(x, g)
    title(titlevars{i})
end

%%
trainingdata = alldosedataclean(:,3:236);
function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236', 'column_237'});

predictorNames = {'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_237;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
classificationSVM = fitcsvm(...
    predictors, ...
    response, ...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 3, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'ClassNames', [0; 1]);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2018a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 234 columns because this model was trained using 234 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236', 'column_237'});

predictorNames = {'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10', 'column_11', 'column_12', 'column_13', 'column_14', 'column_15', 'column_16', 'column_17', 'column_18', 'column_19', 'column_20', 'column_21', 'column_22', 'column_23', 'column_24', 'column_25', 'column_26', 'column_27', 'column_28', 'column_29', 'column_30', 'column_31', 'column_32', 'column_33', 'column_34', 'column_35', 'column_36', 'column_37', 'column_38', 'column_39', 'column_40', 'column_41', 'column_42', 'column_43', 'column_44', 'column_45', 'column_46', 'column_47', 'column_48', 'column_49', 'column_50', 'column_51', 'column_52', 'column_53', 'column_54', 'column_55', 'column_56', 'column_57', 'column_58', 'column_59', 'column_60', 'column_61', 'column_62', 'column_63', 'column_64', 'column_65', 'column_66', 'column_67', 'column_68', 'column_69', 'column_70', 'column_71', 'column_72', 'column_73', 'column_74', 'column_75', 'column_76', 'column_77', 'column_78', 'column_79', 'column_80', 'column_81', 'column_82', 'column_83', 'column_84', 'column_85', 'column_86', 'column_87', 'column_88', 'column_89', 'column_90', 'column_91', 'column_92', 'column_93', 'column_94', 'column_95', 'column_96', 'column_97', 'column_98', 'column_99', 'column_100', 'column_101', 'column_102', 'column_103', 'column_104', 'column_105', 'column_106', 'column_107', 'column_108', 'column_109', 'column_110', 'column_111', 'column_112', 'column_113', 'column_114', 'column_115', 'column_116', 'column_117', 'column_118', 'column_119', 'column_120', 'column_121', 'column_122', 'column_123', 'column_124', 'column_125', 'column_126', 'column_127', 'column_128', 'column_129', 'column_130', 'column_131', 'column_132', 'column_133', 'column_134', 'column_135', 'column_136', 'column_137', 'column_138', 'column_139', 'column_140', 'column_141', 'column_142', 'column_143', 'column_144', 'column_145', 'column_146', 'column_147', 'column_148', 'column_149', 'column_150', 'column_151', 'column_152', 'column_153', 'column_154', 'column_155', 'column_156', 'column_157', 'column_158', 'column_159', 'column_160', 'column_161', 'column_162', 'column_163', 'column_164', 'column_165', 'column_166', 'column_167', 'column_168', 'column_169', 'column_170', 'column_171', 'column_172', 'column_173', 'column_174', 'column_175', 'column_176', 'column_177', 'column_178', 'column_179', 'column_180', 'column_181', 'column_182', 'column_183', 'column_184', 'column_185', 'column_186', 'column_187', 'column_188', 'column_189', 'column_190', 'column_191', 'column_192', 'column_193', 'column_194', 'column_195', 'column_196', 'column_197', 'column_198', 'column_199', 'column_200', 'column_201', 'column_202', 'column_203', 'column_204', 'column_205', 'column_206', 'column_207', 'column_208', 'column_209', 'column_210', 'column_211', 'column_212', 'column_213', 'column_214', 'column_215', 'column_216', 'column_217', 'column_218', 'column_219', 'column_220', 'column_221', 'column_222', 'column_223', 'column_224', 'column_225', 'column_226', 'column_227', 'column_228', 'column_229', 'column_230', 'column_231', 'column_232', 'column_233', 'column_234', 'column_235', 'column_236'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_237;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationSVM, 'KFold', 5);

% Compute validation predictions
[~, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

end


%%

load("MatlabWorkspaceRevision.mat")


%Predict and compare Different Cell behaviors

traceclustermat = NaN(5000,30);


[max(med01(:,239)),max(med02(:,239)),max(med03(:,239)),max(med04(:,239)),max(med05(:,239)),max(med06(:,239)),max(med07(:,239)),max(med08(:,239)),max(med09(:,239)),max(med10(:,239)),max(med11(:,239)),max(med12(:,239)),max(med13(:,239)),max(med14(:,239)),max(med15(:,239))]

meddataclean1 = meddataclean;
meddataclean1(99:149,239)=meddataclean1(99:149,239)+111;
meddataclean1(150:237,239)=meddataclean1(150:237,239)+111 +59;
meddataclean1(238:304,239)=meddataclean1(238:304,239)+111+59+105;
meddataclean1(305:371,239)=meddataclean1(305:371,239)+111+59+105+73;
meddataclean1(372:462,239)=meddataclean1(372:462,239)+111+59+105+73+89;
meddataclean1(463:554,239)=meddataclean1(463:554,239)+111+59+105+73+89+109;
meddataclean1(555:648,239)=meddataclean1(555:648,239)+111+59+105+73+89+109+100;
meddataclean1(649:708,239)=meddataclean1(649:708,239)+111+59+105+73+89+109+100+120;
meddataclean1(709:765,239)=meddataclean1(709:765,239)+111+59+105+73+89+109+100+120+77;
meddataclean1(766:826,239)=meddataclean1(766:826,239)+111+59+105+73+89+109+100+120+77+64;
meddataclean1(827:884,239)=meddataclean1(827:884,239)+111+59+105+73+89+109+100+120+77+64+70;
meddataclean1(885:958,239)=meddataclean1(885:958,239)+111+59+105+73+89+109+100+120+77+64+70+69;
meddataclean1(959:1048,239)=meddataclean1(959:1048,239)+111+59+105+73+89+109+100+120+77+64+70+69+84;
meddataclean1(1049:1109,239)=meddataclean1(1049:1109,239)+111+59+105+73+89+109+100+120+77+64+70+69+84+105;

meddataclean1 = horzcat(meddataclean1,medclusterid(meddataclean1(:,239),2));

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==1))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==2))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==3))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==4))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==5))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);

dummy1=zeros(1109,1);
dummy1(find(meddataclean1(:,242)==6))=1;
meddataclean1 = horzcat(meddataclean1,dummy1);


medmltest = meddataclean1(:,3:236);
medmltest1 = horzcat(medmltest,meddataclean1(:,243));
medmltest1_2 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 2),:)
medmltest1_3 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 3),:)
medmltest1_4 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 4),:)
medmltest1_5 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 5),:)
medmltest1_6 = medmltest1(find(meddataclean1(:,242) == 1 | meddataclean1(:,242) == 6),:)

medmltest2 = horzcat(medmltest,meddataclean1(:,244));
medmltest2_3 = medmltest2(find(meddataclean1(:,242) == 2 | meddataclean1(:,242) == 3),:)
medmltest2_4 = medmltest2(find(meddataclean1(:,242) == 2 | meddataclean1(:,242) == 4),:)
medmltest2_5 = medmltest2(find(meddataclean1(:,242) == 2 | meddataclean1(:,242) == 5),:)
medmltest2_6 = medmltest2(find(meddataclean1(:,242) == 2 | meddataclean1(:,242) == 6),:)

medmltest3 = horzcat(medmltest,meddataclean1(:,245));
medmltest3_4 = medmltest3(find(meddataclean1(:,242) == 3 | meddataclean1(:,242) == 4),:)
medmltest3_5 = medmltest3(find(meddataclean1(:,242) == 3 | meddataclean1(:,242) == 5),:)
medmltest3_6 = medmltest3(find(meddataclean1(:,242) == 3 | meddataclean1(:,242) == 6),:)

medmltest4 = horzcat(medmltest,meddataclean1(:,246));
medmltest4_5 = medmltest4(find(meddataclean1(:,242) == 4 | meddataclean1(:,242) == 5),:)
medmltest4_6 = medmltest4(find(meddataclean1(:,242) == 4 | meddataclean1(:,242) == 6),:)

medmltest5 = horzcat(medmltest,meddataclean1(:,247));
medmltest5_6 = medmltest5(find(meddataclean1(:,242) == 5 | meddataclean1(:,242) == 6),:)

medmltest6 = horzcat(medmltest,meddataclean1(:,248));


figure()
hold on
plot(nanmean(traceclustermat(find(medclusterid(:,2)==1),:)),'r')
plot(nanmean(traceclustermat(find(medclusterid(:,2)==2),:)),'k')
plot(nanmean(traceclustermat(find(medclusterid(:,2)==3),:)),'b')
plot(nanmean(traceclustermat(find(medclusterid(:,2)==4),:)),'g')
plot(nanmean(traceclusterinactive),'y')
plot(nanmean(traceclustermat(find(medclusterid(:,2)==6),:)),'m')



%show variables associated with trace differences
boxplot(meddataclean1(:,11),meddataclean1(:,242),'OutlierSize',3,'Symbol', 'ko')

boxplot(meddataclean1(:,94),meddataclean1(:,242),'OutlierSize',3,'Symbol', 'ko')
   
boxplot(meddataclean1(:,12),meddataclean1(:,242),'OutlierSize',3,'Symbol', 'ko')


% Corr plot
corrmatrix3 = [];
for i = 1:235
    corrmatrix3 = [corrmatrix3,corr(medmltest3_4(:,i),medmltest3_4(:,235))];
end
[~,indexcorr3] = sort(corrmatrix3);

corrsorted3 = corrmatrix3(indexcorr3);

figure(); hold on
barh(1:20,corrsorted3(1:20),'FaceColor',[.2,0.2,0.2],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)
barh(21:41,corrsorted3(207:227),'FaceColor',[1,0,0],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)


%Use Classification ensemble to visualize coupling

view(trainedModel.ClassificationEnsemble.Trained{8},'Mode','graph');


%find accuracy through repeated thresholded
figure()
box(alldosedataclean(:,9),activedummy+1)

couplingactive = alldosedataclean(:,9) > 3625.57;
classp = classperf(activedummy,couplingactive)

couplingtrue = activedummy(find(alldosedataclean(:,8)<1.58e06),:);
couplingactive = alldosedataclean(find(alldosedataclean(:,8)<1.58e06),:);
couplingactive = couplingactive(:,9) > 3625.57;
classp1 = classperf(couplingtrue,couplingactive)

couplingtrue = activedummy(find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741),:);
couplingactive = alldosedataclean(find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741),:);
couplingactive = couplingactive(:,9) > 3625.57;
classp2 = classperf(couplingtrue,couplingactive)

couplingtrue = activedummy(find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741 & alldosedataclean(:,4)<114.12),:);
couplingactive = alldosedataclean(find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741 & alldosedataclean(:,4)<114.12),:);
couplingactive = couplingactive(:,9) > 3625.57;

couplingdummy = zeros(3456,1);
couplingdummy1 = find(alldosedataclean(:,8)<1.58e06 & alldosedataclean(:,3)<741 & alldosedataclean(:,4)<114.12 & alldosedataclean(:,9) < 3625.57);
for i = 1:size(couplingdummy1)
    couplingdummy(couplingdummy1(i)) = 1;
end

%Coupled population corr chart

corrmatrix1 = [];
for i = 1:236
    corrmatrix1 = [corrmatrix1,corr(alldosedataclean(:,i),couplingdummy)];
end
[~,indexcorr1] = sort(corrmatrix1);

corrsorted1 = corrmatrix1(indexcorr1);

figure(); hold on
barh(1:20,corrsorted1(1:20),'FaceColor',[.2,0.2,0.2],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)
barh(21:41,corrsorted1(208:228),'FaceColor',[1,0,0],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)

%Visualize coupled population
   
figure()
hold on
plot(embedding(couplingdummy==1,1),embedding(couplingdummy==1,2),'r' ,'LineStyle','none','Marker','o','MarkerSize', 1)
plot(embedding(couplingdummy==0,1),embedding(couplingdummy==0,2), 'Color',[.7 .7 .7 .1],'LineStyle','none','Marker','o','MarkerSize', 1)
width=400;
height=400;
set(gcf,'position',[x0,y0,width,height])
axis([-8 7 -7 7])
axis off
title('Coupling')




corrmatrix2 = [];
for i = 1:236
    corrmatrix2 = [corrmatrix2,corr(alldosedataclean(:,i),alldosedataclean(:,3))];
end
[~,indexcorr2] = sort(corrmatrix2);

corrsorted2 = corrmatrix2(indexcorr2);

figure(); hold on
barh(1:20,corrsorted2(1:20),'FaceColor',[.2,0.2,0.2],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)
barh(21:41,corrsorted2(208:228),'FaceColor',[1,0,0],...
       'EdgeColor',[1,1,1],'LineWidth',0.00001)


classperf(couplingtrue,couplingactive)






   


>>>>>>> a38de5df5e103f6cadaa2ad8a65cc281ad9984cc
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
