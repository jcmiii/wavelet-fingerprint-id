function [q] = diagPredictorNhbrs(D, P, V, H)
%Creates matrix Q for a vertical subband D.
%Q holds the pixels used for the vertical linear approximation
%See Buccigrossi paper. 
%IMPORTANT NOTE: His horizonatal band is our vertical band, and vice-versa 
%  D = Diagonal subband D_k
%  P = Parent D_{k+1}
%  H = Horizontal subband H_k (lower left)
%  V = Vertical subband V_k (upper right)

   %Get dimensions of subband
  [m,n]=size(D);
    
  i = 1;  %initialize index counter
  %Build q
  for row = 3 : m - 2     %stay away from edges
    disp('D: ')
    row
    for col = 3 : n - 2
      q(i,:) = [D(row-1, col), D(row, col-1), P(ceil(row/2), ceil(col/2)), V(row, col), H(row, col), D(row, col-2)];
      i = i+1;
    end;
  end;
