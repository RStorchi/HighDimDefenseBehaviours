function[] = calculate_SI_batch_cl(Kneigh,allstim)
if ~allstim
    indS = {[0 1], [0 2], [1 2]}; 
    epoch = {'early', 'intermediate', 'late'}; 
else
    indS = {[0 1 2]};     
    epoch = {'all'}; 
end
Ns = numel(indS); Nx = numel(epoch);
count = 0; 
disp(sprintf('K: %s',num2str(Kneigh)));
for n = 1:Nx
    for m = 1:Ns
        count = count+1;
        [spec(:,:,count),spec_loc(:,:,count)] = calculate_SI_cl(indS{m},Kneigh,epoch{n},0);
        disp(sprintf('Iter: %s',num2str(count)));
    end
end
if ~allstim
    save(['SI_K' num2str(Kneigh) '_cl']);
else
    save(['SI_K' num2str(Kneigh) '_early_inter_cl']);
end