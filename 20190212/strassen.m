function [ C ] = strassen( A, B )
%STRASSEN Summary of this function goes here
%   Implements the Strassen MMA which achieves an omega of X
%   A and B are square n x n 


% Partition A and B into equal size block matrices
% C = AB
% A = [A1,1,  A1,2], B = [B1,1,  B1,2], C = [C1,1,  C1,2]
%     [A2,1,  A2,2],     [B2,1,  B2,2],     [C2,1,  C2,2]
%
% If matrices are not of type 2^n x 2^n, then fill the missing rows and
% columns with zeros

% Check lengths of A and B
sizeA = size(A);
sizeB = size(B);
if (length(sizeA) ~= 2) || (length(sizeB) ~= 2)
    sprintf('Size of matrices are not 2 dimensional');
    return
elseif dim(A) ~= dim(B)
    sprintf('Matrices A and B are of different size');
    return
end

% Pad with zeros if necessary
% n = log2(a); e.g. 3 = log2(8) = log2(2^3)
% Use identity: log2(a) = log(a)/log(2)
n0 = log2(sizeA(1,1));
if isinteger(n0)
    % We have a valid matrix in the form 2^n
    additional_rows_cols = 0;
    sprintf('Valid matrix size');
else
    total_rows_cols = 2^ceil(n0);
    additional_rows_cols = total_rows_cols-sizeA(1,1);
end
A = [A, zeros(sizeA(1,1),additional_rows_cols); zeros(additional_rows_cols, total_rows_cols)];
B = [B, zeros(sizeB(1,1),additional_rows_cols); zeros(additional_rows_cols, total_rows_cols)];
n1 = log2(dim(A));

% Let's only do 1 division
A11 = A(1:2^(n1-1),1:2^(n1-1));
A12 = A(1:2^(n1-1),2^(n1-1)+1:end);
A21 = A(2^(n1-1)+1:end,1:2^(n1-1));
A22 = A(2^(n1-1)+1:end,2^(n1-1)+1:end);
B11 = B(1:2^(n1-1),1:2^(n1-1));
B12 = B(1:2^(n1-1),2^(n1-1)+1:end);
B21 = B(2^(n1-1)+1:end,1:2^(n1-1));
B22 = B(2^(n1-1)+1:end,2^(n1-1)+1:end);

% if n1 == 1
    M1 = (A11+A22)*(B11+B22);
    M2 = (A21+A22)*(B11);
    M3 = (A11)*(B12-B22);
    M4 = (A22)*(B21-B11);
    M5 = (A11+A12)*(B22);
    M6 = (A21-A11)*(B11+B12);
    M7 = (A12-A22)*(B21+B22);
% else
%     M1 = strassen(A11+A22,B11+B22);
%     M2 = strassen(A21+A22,B11);
%     M3 = strassen(A11,B12-B22);
%     M4 = strassen(A22,B21-B11);
%     M5 = strassen(A11+A12,B22);
%     M6 = strassen(A21-A11,B11+B12);
%     M7 = strassen(A12-A22,B21+B22);
% end


C = [M1+M4-M5+M7, M3+M5;...
    M2+M4, M1-M2+M3+M6];
C = C(1:sizeA(1,1),1:sizeA(1,1));

end
