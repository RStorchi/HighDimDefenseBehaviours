function[P] = DLT_simple_cl(x,X)
%This function use direct linear method (svd) to find the coefficients P
%of the camera matrix.
%Use the fact that for any vectot x the cross product  x cprod x = null vector,
%therefore x cprod PX = null vec if P coefficients are correct. 
%The linear system is composed of 3*n equation, each is equation as
%follows:
%given x = [x1 y1 w1]' 
%    [ 0      -w1*X'   y1*X']   [P1]
%A = [ w1*X'   0      -x1*X'] * [P2] = 0  
%    [-y1*X'   x1*X'   0    ]   [P3]
%and the vector p = [P1' P2' P3']' is made so that P1' = first row of P matrix, 
%P2' = 2nd row, P3' third row. The system is (3*12) * (12*1) elements
%I use SVD so that the vector p is constrained to be |p| = 1

%generate full 3N*12 matrix
N = size(x,2);
A = [];
for n = 1:N
    x1 = x(1,n); y1 = x(2,n); w1 = x(3,n);
    A1 = [zeros(1,4)      -w1*X(:,n)'      y1*X(:,n)'    ];
    A2 = [ w1*X(:,n)'     zeros(1,4)      -x1*X(:,n)'    ];
    A3 = [-y1*X(:,n)'      x1*X(:,n)'     zeros(1,4)     ];
    A  = [A; A1; A2; A3];
end



%apply svd
[U,S,V] = svd(A);
%get the V associated with the smallest singular value (the last one)
p = V(:,end); %reduced vector
P = [p(1:4)'; p(5:8)'; p(9:12)'];







