function[] = generate_world_coordinates()
%generate real world coordinates for Lego(R) objects
h0 = [0 3 5 8 10]; 
Nh = length(h0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x,y coordinates camera 1&2
%%%%%%%%%%%%%%%%%%%%%%%%%%
y0 = [1 7 13 20 26 32]; y0 = y0-mean(y0); y0 = y0; Ny0 = length(y0); 
x0 = [1 13 20 32]; x0 = x0-mean(x0); x0 = x0; Nx0 = length(x0);
Np = Nx0*Ny0*Nh;
%load all points
X = [];
for h = 1:Nh
    for n = 1:Nx0
        for m = 1:Ny0
            X = [X [x0(n); y0(m); h0(h); 1]];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x,y coordinates camera 3&4
%%%%%%%%%%%%%%%%%%%%%%%%%%
y0 = [1 7 13 20 26 32]; y0 = y0-mean(y0); y0 = y0; Ny0 = length(y0); 
x0 = [1 13 20 32]; x0 = x0-mean(x0); x0 = x0; Nx0 = length(x0);
Np = Nx0*Ny0*Nh;
%load all points
X1 = [];
for h = 1:Nh
    for n = 1:Nx0
        for m = 1:Ny0
            X1 = [X1 [x0(n); y0(m); h0(h); 1]];
        end
    end
end

%fig
figure; plot3(X(1,:),X(2,:),X(3,:),'.')
%save
%save('world_coordinates','X','X1');