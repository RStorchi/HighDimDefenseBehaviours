function[] = visualize_trial_3D_cl(Y, el)
Nframe = size(Y,3);

%transform to standard coordinate
yref = [0 3.5; -1 2; 1 2; 0 1; 0 -3.5];
[~, ~, t] = procrustes(yref, squeeze(Y(:,1:2,1)), 'Scaling', false,'Reflection',false);
for n = 1:Nframe
    y = squeeze(Y(:,1:2,n));
    y = y*t.T + t.c;
    Y(:,1:2,n) = y;
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%generate final images%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
set(fig,'Position',[100 100 1200 500]);
for n = 1:2
    h(n) = subplot(1,2,n);
end
cgray_vec = 0.7*[1:-1/(Nframe-1):0]'; 
cgray = repmat(cgray_vec,1,3);
%
x0 = mean(mean(Y(:,1,:))); dx0 = 10; 
y0 = mean(mean(Y(:,2,:))); dy0 = 10; 
z0 = mean(mean(Y(:,3,:))); dz0 = 10; 
%
for n = 1:Nframe
    set(fig,'Color',[0.8 0.8 0.8]);
    %
    subplot(h(1)); hold on;
    plot3(Y(1,1,n),Y(1,2,n),Y(1,3,n),'.r','MarkerSize',15);
    plot3(Y(2,1,n),Y(2,2,n),Y(2,3,n),'.b','MarkerSize',15);
    plot3(Y(3,1,n),Y(3,2,n),Y(3,3,n),'.b','MarkerSize',15);
    plot3(Y(4,1,n),Y(4,2,n),Y(4,3,n),'.g','MarkerSize',15);
    plot3(Y(5,1,n),Y(5,2,n),Y(5,3,n),'.y','MarkerSize',15);
    line([Y(1,1,n) Y(2,1,n)],[Y(1,2,n) Y(2,2,n)],[Y(1,3,n) Y(2,3,n)],'Color',cgray(n,:),'LineWidth',2);
    line([Y(1,1,n) Y(3,1,n)],[Y(1,2,n) Y(3,2,n)],[Y(1,3,n) Y(3,3,n)],'Color',cgray(n,:),'LineWidth',2);
    line([Y(2,1,n) Y(4,1,n)],[Y(2,2,n) Y(4,2,n)],[Y(2,3,n) Y(4,3,n)],'Color',cgray(n,:),'LineWidth',2);
    line([Y(3,1,n) Y(4,1,n)],[Y(3,2,n) Y(4,2,n)],[Y(3,3,n) Y(4,3,n)],'Color',cgray(n,:),'LineWidth',2);
    line([Y(4,1,n) Y(5,1,n)],[Y(4,2,n) Y(5,2,n)],[Y(4,3,n) Y(5,3,n)],'Color',cgray(n,:),'LineWidth',2);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view([el 90]); 
    xlim([x0-dx0 x0+dx0]); ylim([y0-dy0 y0+dy0]); 
    zlim([z0-dz0 z0+dz0]); 
    %
    subplot(h(2)); hold on;
    plot3(Y(1,1,n),Y(1,2,n),Y(1,3,n),'.r','MarkerSize',15);
    plot3(Y(2,1,n),Y(2,2,n),Y(2,3,n),'.b','MarkerSize',15);
    plot3(Y(3,1,n),Y(3,2,n),Y(3,3,n),'.b','MarkerSize',15);
    plot3(Y(4,1,n),Y(4,2,n),Y(4,3,n),'.g','MarkerSize',15);
    plot3(Y(5,1,n),Y(5,2,n),Y(5,3,n),'.y','MarkerSize',15);
    line([Y(1,1,n) Y(2,1,n)],[Y(1,2,n) Y(2,2,n)],[Y(1,3,n) Y(2,3,n)],'Color',cgray(n,:),'LineWidth',2);
    line([Y(1,1,n) Y(3,1,n)],[Y(1,2,n) Y(3,2,n)],[Y(1,3,n) Y(3,3,n)],'Color',cgray(n,:),'LineWidth',2);
    line([Y(2,1,n) Y(4,1,n)],[Y(2,2,n) Y(4,2,n)],[Y(2,3,n) Y(4,3,n)],'Color',cgray(n,:),'LineWidth',2);
    line([Y(3,1,n) Y(4,1,n)],[Y(3,2,n) Y(4,2,n)],[Y(3,3,n) Y(4,3,n)],'Color',cgray(n,:),'LineWidth',2);
    line([Y(4,1,n) Y(5,1,n)],[Y(4,2,n) Y(5,2,n)],[Y(4,3,n) Y(5,3,n)],'Color',cgray(n,:),'LineWidth',2);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view([el 0]); 
    xlim([x0-dx0 x0+dx0]); ylim([y0-dy0 y0+dy0]); 
    zlim([z0-dz0 z0+dz0]); 
end