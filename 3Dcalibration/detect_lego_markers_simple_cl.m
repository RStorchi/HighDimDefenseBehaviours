function[] = detect_lego_markers_simple_cl(which_camera,which_height)
%load to image 
x0 = imread(['camera' num2str(which_camera) '_lego_h' num2str(which_height) '.tif']);
%init parameters
Np = 24; Nl = 4; Npl = Np/Nl;
%smooth image
[Nrow0,Ncol0] = size(x0(:,:,1));
x0 = double(x0); 
for n = 1:3
    x0(:,:,n) = filter2(ones(3)/9,x0(:,:,n)+20*ones(Nrow0,Ncol0));
end
%reduce to red channel
x0 = x0(:,:,1);

%zero the contours
fig0 = figure; 
imshow(uint8(x0)); 
title('Click Left Bottom & Right Top Corner','FontSize',12);
[idx_null,idy_null] = ginput(2);
x0_mask = false*zeros(Nrow0,Ncol0);
x0_mask(min(idy_null):max(idy_null),min(idx_null):max(idx_null)) = true;
x0 = x0.*x0_mask;
imshow(uint8(x0)); 

%get all points
Nmask = 15;
x = zeros(Nrow0+2*Nmask, Ncol0+2*Nmask);
x(Nmask+1:Nmask+Nrow0,Nmask+1:Nmask+Ncol0) = x0;

%find markers
for n = 1:Np
    [row_max, col_max] = find(x == max(x(:)),1,'first');
    x(row_max-Nmask+1:row_max+Nmask, col_max-Nmask+1:col_max+Nmask) = 0;
    idx(n) = row_max; idy(n) = col_max; 
end
idx = idx-Nmask; idy = idy-Nmask; 
fig1 = figure; 
imshow(uint8(x0)); hold on;
plot(idy,idx,'xr','MarkerSize',15);

%order markers
idx_final = []; idy_final = [];
for n = 1:Nl
    title(['Two-point click for row ' num2str(n)],'FontSize',12);
    [idx_l,idy_l] = ginput(2);
    %find line equation
    a = diff(idy_l)/diff(idx_l);
    b = idy_l(1)-a*idx_l(1);
    plot([1:Ncol0],a*[1:Ncol0]+b,'b');
    %find points along the line
    idx_est = a*idy+b;
    dist = abs(idx - idx_est)
    [dist_sort, id_sort] = sort(dist,'ascend');
    id_sort = id_sort(1:Npl);
    for m = 1:Npl
        text(idy(id_sort(m)),idx(id_sort(m))+10,num2str(id_sort(m)),'FontSize',14,'Color',[0 1 0]);
    end
    %order points along the line
    [ysort,id_sort2] = sort(idx(id_sort),'descend');
    idx_final = [idx_final idx(id_sort(id_sort2))];
    idy_final = [idy_final idy(id_sort(id_sort2))];
end

%final check
fig2 = figure; 
imshow(uint8(x0)); hold on;
plot(idy_final,idx_final,'xk','MarkerSize',15);
for n = 1:Np
    text(idy_final(n),idx_final(n)+10,num2str(n),'FontSize',14,'FontWeight','bold','Color',[0 1 0]);
end

%save figure with ordered markers
%saveas(fig2,['camera' num2str(which_camera) '_lego_h' num2str(which_height) '_marked_cl.tif'],'tif');

%save position data
%save(['camera' num2str(which_camera) '_lego_h' num2str(which_height) '_data_cl'],'idx_final','idy_final');