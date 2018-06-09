function x_com = compress_matrix(x,n)
%% This function compress large matrix by compression ratio n.
%  x : oritinal matrix, data ordered in each column
%  n : compression ratio
%  x_com : compressed matrix
%
%	Sihao Sun 17-Apr-2017
%	S.Sun-4@tudelft.nl

%% 
[n_row,n_clm] = size(x);
n_com = floor(n_row/n);
x_com = zeros(n_com,n_clm);
for j = 1:n_clm
    for i = 1:n_com  
        x_com(i,j) = x(i*n,j);
    end
end
end