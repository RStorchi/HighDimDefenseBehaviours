function[pval] = figure_specificity_all_cl()
%Plot the SI index for locomotion and full set. Returns pval, 
%a 1 x 9 array associated with each comparison shown in th figure.
%The p-values are calculated with sign-tests.  

Kneigh = [1 2 4 8 16 32 64 128];
Nk = numel(Kneigh);
which_ind = 1:9; %all latency & stimuli combinations
Nind = numel(which_ind);
for n = 1:Nind
    [SI(:,n),SIloc(:,n)] = figure_specificity_all_cl_sub1(Kneigh,which_ind(n));    
end
SI = SI(:,[1 4 7 2 5 8 3 6 9]);
SIloc = SIloc(:,[1 4 7 2 5 8 3 6 9]);
Ntrial = size(SI,1);
mSI = mean(SI); sSI = std(SI);  
mSIloc = mean(SIloc); sSIloc = std(SIloc);  
%stats
for n = 1:Nind
    pval(n) = signtest(SI(:,n),SIloc(:,n));
end
%figure
fig = figure; set(fig,'Position',[100 100 600 275]);
h = subplot(1,1,1); hold on; 
cf = 'kkkrrrbbb'; ce = 'bbbkkkrrr';
XTickLab = {'0-1s','1-2s','2-3s','0-1s','1-2s','2-3s','0-1s','1-2s','2-3s'};
XTick = [];
for n = 1:Nind
    nb = (n>3)+(n>6);
    XTick = [XTick n+nb];
    if mSI(n)>mSIloc(n)
        bar(n+nb,mSI(n),'BarWidth',0.6,'LineWidth',2,'EdgeColor',ce(n),'FaceColor',cf(n)); 
        bar(n+nb,mSIloc(n),'BarWidth',0.6,'LineWidth',2,'EdgeColor',ce(n),'FaceColor',0.8*ones(1,3));
    else
        bar(n+nb,mSIloc(n),'BarWidth',0.6,'LineWidth',2,'EdgeColor',ce(n),'FaceColor',0.8*ones(1,3));
        bar(n+nb,mSI(n),'BarWidth',0.6,'LineWidth',2,'EdgeColor',ce(n),'FaceColor',cf(n)); 
    end
end
errorbar(XTick,mSI,sSI/sqrt(Ntrial),'.','LineWidth',2,'Color',0.6*ones(1,3));
errorbar(XTick,mSIloc,sSIloc/sqrt(Ntrial),'.','LineWidth',2,'Color',0.6*ones(1,3));
ylim([0.5 1]);
ylabel('SI','FontSize',14); 
set(h,'XTick',XTick,'XTickLabel',XTickLab);

function[SI,SIloc] = figure_specificity_all_cl_sub1(Kneigh,ind)
Nk = numel(Kneigh); 
for n = 1:Nk
    load(['SIinv_review_K' num2str(Kneigh(n)) '_v1.mat'],'spec','spec_loc');
    spec_all(:,:,n) = squeeze(spec(:,:,ind)); 
    spec_all_loc(:,:,n) = squeeze(spec_loc(:,:,ind));
end
spec = spec_all; spec_loc = spec_all_loc; 
clear spec_all; clear spec_all_loc;
%find best combination of K and Npca
mspec = squeeze(mean(spec,1));
mspec_loc = squeeze(mean(spec_loc,1));
[idx,idy] = find(mspec==max(mspec(:)),1);
[idx_loc,idy_loc] = find(mspec_loc==max(mspec_loc(:)),1);
%resize 
if idy~=idy_loc
    idy_loc = idy;
end
SI = spec(:,idx,idy);
SIloc = spec_loc(:,idx_loc,idy_loc);
%SI = spec(:,idx,idy);
%SIloc = spec_loc(:,idx_loc,idy_loc);



