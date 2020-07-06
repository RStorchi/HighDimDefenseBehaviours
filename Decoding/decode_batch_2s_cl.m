function[pval,frac] = decode_batch_2s_cl(calculate, make_fig)
rng(1);
%calculate and display knn decoding results for 0-2s epochs
epoch = 'all'; %0-2s
Nrep = 50; 
%
which_dim{1} = 1:9; which_stim{1} = [0 1];
which_dim{2} = 4; which_stim{2} = [0 1];
%
which_dim{3} = 1:9; which_stim{3} = [0 2];
which_dim{4} = 4; which_stim{4} = [0 2];
%
which_dim{5} = 1:9; which_stim{5} = [1 2];
which_dim{6} = 4; which_stim{6} = [1 2];
%
which_dim{7} = 1:9; which_stim{7} = [0 1 2];
which_dim{8} = [1 4:7]; which_stim{8} = [0 1 2];
which_dim{9} = 4; which_stim{10} = [0 1 2];
%
Ndec = numel(which_dim);
%
if calculate
    correct = []; 
    for n = 1:Ndec
        correct{n} = decode_cl(epoch, which_stim{n}, which_dim{n}, Nrep);
    end
    save('results_stimulus_predict_0_2s_knn_review');
end
%
if make_fig
    %make figures 
    load('results_stimulus_predict_0_2s_knn_review');
    for n = 1:Ndec
        mc_temp = squeeze(mean(100*correct{n},1));
        sc_temp = squeeze(std(100*correct{n},1));
        [ind_best_row, ind_best_col] = find(mc_temp == max(mc_temp(:)),1);
        mc(n) = mc_temp(ind_best_row, ind_best_col);
        sc(n) = sc_temp(ind_best_row, ind_best_col);
    end
    %figure 0-2s pairwise comparisons
    fig = figure; set(fig,'Position',[100 100 250 300]);
    h = subplot(1,1,1); hold on; 
    xval = [1 1 2 2 3 3];
    cf = [0 0 0; 0.8*ones(1,3); 1 0 0; 0.8*ones(1,3); 0 0 1; 0.8*ones(1,3)];
    ce = [0 0 1; 0 0 1; 0 0 0; 0 0 0; 1 0 0; 1 0 0];
    for n = 1:6
        bar(xval(n),mc(n),'BarWidth',0.5,'FaceColor',cf(n,:),'EdgeColor',ce(n,:),'LineWidth',2);
        errorbar(xval(n),mc(n),sc(n),'.','LineWidth',2,'Color',[1 1 1]*0.6);
    end
    xlim([0 4]); ylim([50 100]);
    %figure 0-2s all 3 stimuli comparisons
    fig = figure; set(fig,'Position',[100 100 250 300]);
    h = subplot(1,1,1); hold on; 
    xval = [1 2 1 2];
    cf = [0.85*[1 0 1]; 0.65*[1 0 1]; 0.8*ones(1,3); 0.8*ones(1,3)];
    ce = zeros(6,3);
    n0 = 6;
    for n = 1:2
        bar(xval(n),mc(n+n0),'BarWidth',0.5,'FaceColor',cf(n,:),'EdgeColor',ce(n,:),'LineWidth',2);
        errorbar(xval(n),mc(n+n0),sc(n+n0),'.','LineWidth',2,'Color',[1 1 1]*0.6);
    end
    for n = 3:4
        bar(xval(n),mc(9),'BarWidth',0.5,'FaceColor',cf(n,:),'EdgeColor',ce(n,:),'LineWidth',2);
        errorbar(xval(n),mc(9),sc(9),'.','LineWidth',2,'Color',[1 1 1]*0.6);
    end
    xlim([0 3]); ylim([33.333 100]);
    %%%%%%%%%%%%calculate statistics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %flash vs loom
    Ntrial = 344;
    s = round(mc(1)*Ntrial/100);
    p = mc(2)/100;
    pval.flash_loom = myBinomTest(s,Ntrial,p,'two');
    frac.flash_loom = 100*(mc(1)-mc(2))/(mc(2)-50);
    %flash vs sound
    Ntrial = 344;
    s = round(mc(3)*Ntrial/100);
    p = mc(4)/100;
    pval.flash_sound = myBinomTest(s,Ntrial,p,'two');
    frac.flash_sound = 100*(mc(3)-mc(4))/(mc(4)-50);
    %loom vs sound
    Ntrial = 344;
    s = round(mc(5)*Ntrial/100);
    p = mc(6)/100;
    pval.loom_sound = myBinomTest(s,Ntrial,p,'two');
    frac.loom_sound = 100*(mc(5)-mc(6))/(mc(6)-50);
    %full vs locomotion
    Ntrial = 516;
    s = round(mc(7)*Ntrial/100);
    p = mc(9)/100;
    pval.full_loc = myBinomTest(s,Ntrial,p,'two');
    frac.full_loc = 100*(mc(7)-mc(9))/(mc(9)-100/3);
    %no-shape vs locomotion
    Ntrial = 516;
    s = round(mc(8)*Ntrial/100);
    p = mc(9)/100;
    pval.noshape_loc = myBinomTest(s,Ntrial,p,'two');
    frac.noshape_loc = 100*(mc(8)-mc(9))/(mc(9)-100/3);
    %compare no-shape vs full
    Ntrial = 516;
    s = round(mc(8)*Ntrial/100);
    p = mc(7)/100;
    pval.full_no_shape = myBinomTest(s,Ntrial,p,'two');
    frac.full_no_shape = 100*(mc(7)-mc(8))/(mc(8)-100/3);
end

