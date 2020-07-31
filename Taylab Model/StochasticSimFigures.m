<<<<<<< HEAD
%Simulation figures:

load('cell1000_tnfsweep.mat')

peakheights = [];
for i = 1:20
    peakheights(i,:) = max(squeeze(XXX(:,i,37:457,8))');
end

figure()
hold on
for i = 1:20
    scatter(TNF_List(i)*ones(1000,1),peakheights(i,:))
end
title('peak heights')
xlabel('log10(tnf)')
ylabel('single cell peak height')


TNF1 = squeeze(XXX(:,1,37,:));
TNF1_norm = zscore(TNF1);
rng(1)
tsnesim = tsne(TNF1_norm,'Algorithm','barneshut','perplexity',20,'Exaggeration', 6);

%5,6,11,12,14,15 (1)
t0_ikb_tot=(TNF1(:,5)+TNF1(:,6)+TNF1(:,11)+TNF1(:,12)+TNF1(:,14)+TNF1(:,15));

% 6,7,8,14,15(3)
t0_nfkb_tot=(TNF1(:,6)+TNF1(:,7)+TNF1(:,8)+TNF1(:,14)+TNF1(:,15));

t0_ikb_nfkb_ratio= (t0_ikb_tot./t0_nfkb_tot);


%Figure4b1 nuc p65
figure()
scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(log(TNF1(:,8))))
colormap(flipud(hot))
caxis([-2 2])
title('nuclear nfkb')

%Figure4b2 ikb nfkb ratio
figure()
scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(t0_ikb_nfkb_ratio))
colormap(flipud(hot))
caxis([-2 2])
title('ikb nfkb ratio')

%Figure4b3 total ikb
figure()
scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(t0_ikb_tot))
colormap(flipud(hot))
caxis([-2 2])
title('total ikb')

% %cyto A20
% figure()
% scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(TNF1(:,9)))
% colormap(flipud(hot))
% caxis([-2 2])
% title('cytoplasmic A20')
% 
% 
% %ikk total
% figure()
% scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(TNF1(:,2)+TNF1(:,3)+TNF1(:,4)))
% colormap(flipud(hot))
% caxis([-2 2])
% title('ikk total')


%Figure 4c
TNFsweep = max(squeeze(XXX(:,:,37:457,8)),[],3);

figure(); hold on
yyaxis left
scatter(log10(TNF1(:,8)),TNFsweep(:,17),3,[.1 .1 .1])
ylabel('peak nfkb')

yyaxis right
scatter(log10(TNF1(:,8)),t0_ikb_nfkb_ratio,3,'r')
ylabel('ikb:nfkb ratio')
xlabel('leaky nfkb')
=======
%Simulation figures:

load('cell1000_tnfsweep.mat')

peakheights = [];
for i = 1:20
    peakheights(i,:) = max(squeeze(XXX(:,i,37:457,8))');
end

figure()
hold on
for i = 1:20
    scatter(TNF_List(i)*ones(1000,1),peakheights(i,:))
end
title('peak heights')
xlabel('log10(tnf)')
ylabel('single cell peak height')


TNF1 = squeeze(XXX(:,1,37,:));
TNF1_norm = zscore(TNF1);
rng(1)
tsnesim = tsne(TNF1_norm,'Algorithm','barneshut','perplexity',20,'Exaggeration', 6);

%5,6,11,12,14,15 (1)
t0_ikb_tot=(TNF1(:,5)+TNF1(:,6)+TNF1(:,11)+TNF1(:,12)+TNF1(:,14)+TNF1(:,15));

% 6,7,8,14,15(3)
t0_nfkb_tot=(TNF1(:,6)+TNF1(:,7)+TNF1(:,8)+TNF1(:,14)+TNF1(:,15));

t0_ikb_nfkb_ratio= (t0_ikb_tot./t0_nfkb_tot);


%Figure4b1 nuc p65
figure()
scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(log(TNF1(:,8))))
colormap(flipud(hot))
caxis([-2 2])
title('nuclear nfkb')

%Figure4b2 ikb nfkb ratio
figure()
scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(t0_ikb_nfkb_ratio))
colormap(flipud(hot))
caxis([-2 2])
title('ikb nfkb ratio')

%Figure4b3 total ikb
figure()
scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(t0_ikb_tot))
colormap(flipud(hot))
caxis([-2 2])
title('total ikb')

% %cyto A20
% figure()
% scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(TNF1(:,9)))
% colormap(flipud(hot))
% caxis([-2 2])
% title('cytoplasmic A20')
% 
% 
% %ikk total
% figure()
% scatter(tsnesim(:,1),tsnesim(:,2),3,zscore(TNF1(:,2)+TNF1(:,3)+TNF1(:,4)))
% colormap(flipud(hot))
% caxis([-2 2])
% title('ikk total')


%Figure 4c
TNFsweep = max(squeeze(XXX(:,:,37:457,8)),[],3);

figure(); hold on
yyaxis left
scatter(log10(TNF1(:,8)),TNFsweep(:,17),3,[.1 .1 .1])
ylabel('peak nfkb')

yyaxis right
scatter(log10(TNF1(:,8)),t0_ikb_nfkb_ratio,3,'r')
ylabel('ikb:nfkb ratio')
xlabel('leaky nfkb')
>>>>>>> e409751e531eb57f3e9bdec15b37d1cdd56cfeee
