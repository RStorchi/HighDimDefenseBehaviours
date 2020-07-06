function[RGB] = generate_equally_saturated_colours(Vval,Sval,N,make_fig)
%e.g.: generate_equally_saturated_colours(1,1,100,true)
H = linspace(0,1,N)';
V = Vval*ones(N,1);
S = Sval*ones(N,1);
HSV = [H S V]; 
RGB = hsv2rgb(HSV);

if make_fig
    fig = figure; 
    set(fig,'Position',[300 300 900 200]);
    hold on;
    for n = 1:N
        plot(n,1,'.','MarkerSize',30,'Color',RGB(n,:));    
    end
    xlim([0 N+1]); ylim([0 2]);
end

