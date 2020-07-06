function[pval] = figure_decoding_all_scatter_cl(which_alg)
%Displays decoding results 
%INPUT:
%which_alg: type 'knn' for K-Nearest Neighbour, type 'rf' for Random
%Forest
if ~(strcmp(which_alg,'knn')|strcmp(which_alg,'rf'))
    disp('Insert valid decoding algorithm (knn/rf)');
    return
end

Kneigh = [1 2 4 8 16 32 64 128];
Nk = numel(Kneigh);
which_stim = {[0 1],[0 2],[1 2]}; 
Nstim = numel(which_stim);
mSI = []; sSI = [];
mSIloc = []; sSIloc = [];
for n = 1:Nstim
    [mean_out, std_out] = sort_decoding_results_cl(which_stim{n},which_alg);
    mSI = [mSI mean_out(1:2:end)]; mSIloc = [mSIloc mean_out(2:2:end)];
    sSI = [sSI std_out(1:2:end)]; sSIloc = [sSIloc std_out(2:2:end)];
end
Nind = numel(mSI);
%figure
fig = figure; set(fig,'Position',[100 100 350 300]);
h = subplot(1,1,1); hold on; 
set(h,'FontSize',12);
line([50 100],[50 100],'Color',0.66*ones(1,3),'LineWidth',2);
cf = 'kkkrrrbbb'; ce = 'bbbkkkrrr'; cm = 'osdosdosd';
for n = 1:Nind
    plot(mSIloc(n),mSI(n),'Marker',cm(n),'LineWidth',2,'MarkerSize',10,'Color',ce(n),'MarkerFaceColor',cf(n)); 
end
xlim([50 100]); ylim([50 100]);
ylabel('Full Set(%Correct)','FontSize',14); 
xlabel('Locomotion(%Correct)','FontSize',14); 
%statistics
Ntrial = 344;
for n = 1:Nind
    s = round(mSI(n)*Ntrial/100);
    p = mSIloc(n)/100;
    pval(n) = myBinomTest(s,Ntrial,p,'two');
end


