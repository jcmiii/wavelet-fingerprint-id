function [q] = vertPredictorNhbrs(V, P, D, DA)
%Creates vector Q for a vertical subband V.
%Q holds the pixels used for the vertical linear approximation
%See Buccigrossi paper. 
%IMPORTANT NOTE: His horizonatal band is our vertical band, and vice-versa 
%  V = Vertical subband V_k (upper right)
%  P = Parent V_{k+1}
%  D = Diagonal subband D_k
% DA = Diagonal Aunt subband D_{k+1}

  %Get dimensions of subband
  [m,n]=size(V);
    
  i = 1;  %initialize index counter
  %Build q
  for row = 3 : m - 2     %stay away from edges
    disp('V: ')
    row
    for col = 3 : n - 2
      q(i,:) = [V(row, col-1), V(row-1, col), P(ceil(row/2), ceil(col/2)), D(row, col), V(row, col-2), DA(ceil(row/2), ceil(col/2))];
      i = i+1;
    end;
  end;
  