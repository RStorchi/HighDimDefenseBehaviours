function[pval,es,rd] = fig_mean_response_cl(is_normalized)
%Shows mean and sem for each measure of postures and movements across the 
%three stimuli. 
%INPUT: 
%is_normalized: type true for normalized data, false for un-normalized ones.
%OUTPUT:
%pval: provides p-values for each pariwise comparison. p-values are
%calculate with shuffle test implemented in test_mean_difference.m.
%es: effect size for distance between mean responses as fraction of 
%distance obtained by shuffling  

%load data
if is_normalized
    load('measures_postures_movements.mat');
else
    load('measures_postures_movements_unnormalized.mat');
    X = X0; clear X0;
end
load('data_3D_all.mat','flag_val');
%
X = X(flag_val(:,1)==1,:,:);
flag_val = flag_val(flag_val(:,1)==1,:);
%
indt = [105:210];
t = (indt-150)/15;
%
stim0 = flag_val(:,5);
stim = stim0;
stim((stim0==1)|(stim0==2)) = 1;
stim((stim0==3)|(stim0==4)) = 2;
which_stim = unique(stim);
%
Nstim = numel(which_stim);
[N,Nt,Ndim] = size(X);
Ntd = Nt/Ndim;
%
fig = figure;
set(fig,'Position',[50 50 700 600]);
cval = [0 0 0; 0 0 1; 1 0 0];
%
for n = 1:Ndim
    h(n) = subplot(3,3,n); hold on;
    for m = 1:Nstim
        ind = find(stim==which_stim(m));
        temp{m} = squeeze(X(ind,indt,n));
        temp1{m} = temp{m}(:,46:75); %0-2s epoch
        temp{m}(125,:) = NaN;
        mtemp(m,:) = nanmean(temp{m});
        stemp(m,:) = nanstd(temp{m})/sqrt(numel(ind));
        errorbar(t,mtemp(m,:),stemp(m,:),'Color',cval(m,:),'LineWidth',2);
    end
    %stats
    [pval(:,:,n),es(:,:,n)] = test_mean_difference_cl(temp1);
    %
    xlim([t(1) t(end)]);
    if is_normalized
        ylim([0.1 0.9]);
    else
        ymin = min(mtemp(:)-stemp(:)); 
        ymax = max(mtemp(:)+stemp(:));
        dy = 2*(ymax-ymin);
        ylim([ymin-0.1*dy ymax+0.1*dy]);
    end
    set(h(n),'FontSize',12); 
    if ((n==7)|(n==8)|(n==9))
        xlabel('Time(s)','FontSize',14);
    end
    if ((n==1)|(n==4)|(n==7))
        if is_normalized
            ylabel('Normalized','FontSize',14); 
        else
            ylabel('Raw','FontSize',14); 
        end
    end
    title(dim_name{n},'FontSize',14);
end
%eff size, pval
fig = figure;
set(fig,'Position',[50 50 300 250]);
cval = [0 0 0; 0 0 1; 1 0 0];
%
for n = 1:Ndim
    h(n) = subplot(3,3,n); hold on;
    count = 1;
    cf = 'krb'; ce = 'bkr';
    for m = 1:Nstim
        for p = m+1:Nstim
            bar(count,es(m,p,n),'BarWidth',0.5,'FaceColor',cf(count),'EdgeColor',ce(count),'LineWidth',2);
            if pval(m,p,n)<0.001
               plot(count,12.5,'*k','MarkerSize',5,'LineWidth',2); 
            end
            count = count+1;
        end
    end
    xlim([0 4]); ylim([0 13]);
    set(h(n),'XTick',[]);
end

%calculate response divergence across stimuli
%flash
temp_es = [squeeze(es(1,2,:)); squeeze(es(1,2,:))];
temp_pval = [squeeze(pval(1,2,:)); squeeze(pval(1,2,:))];
significant = find(temp_pval<0.001);
rd.flash = temp_es(significant);
%loom & sound
temp_es = squeeze(es(2,3,:));
temp_pval = squeeze(pval(2,3,:));
significant = find(temp_pval<0.001);
rd.loom_sound = temp_es(significant);
