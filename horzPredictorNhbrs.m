function [q] = horzPredictorNhbrs(H, P, D, DA)
%Creates vector Q for a vertical subband V.
%Q holds the pixels used for the vertical linear approximation
%See Buccigrossi paper. 
%IMPORTANT NOTE: His horizonatal band is our vertical band, and vice-versa 
%  H = Horizontal subband H_k (lower left)
%  P = Parent H_{k+1}
%  D = Diagonal subband D_k
% DA = Diagonal Aunt subband D_{k+1}

  %Get dimensions of subband
  [m,n]=size(H);
    
  i = 1;  %initialize index counter
  %Build q
  for row = 3 : m - 2     %stay away from edges
    disp('H: ')
    row
    for col = 3 : n - 2
      q(i,:) = [H(row-1, col), H(row, col-1), P(ceil(row/2), ceil(col/2)), D(row, col), H(row-2, col), DA(ceil(row/2), ceil(col/2))];
      i = i+1;
    end;
  end;
