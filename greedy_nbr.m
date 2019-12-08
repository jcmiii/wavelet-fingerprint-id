%Clear all variables from memory
clear

%local variables: edit as needed
numIts = 4;
h=Daub(4);  %coefficients for wavelet transform

%Array holding pic names
picNames={'A1.png','A2.png','A3.png','B1.png','B2.png','B3.png','G1.png','G2.png','G3.png','C1.png','C2.png','C3.png','M1.png','M2.png','M3.png','T1.png','T2.png','T3.png'};
[rows,numPics]=size(picNames);

for picNum = 1 : numPics
    disp(picNames{picNum})
    %Read in an image
    A=ImageRead(picNames{picNum});

    %Perform numIts iterations of wavelet transform
    WA = WT2D(A, h, numIts);
    disp('Completed wavelet transform')

    %Get dimensions of transformed pic (should be same as original)
    [m,n]=size(WA);

    %Get subbands. Store in cell arrays.containers
    for k = 1 : numIts
        V{k} = WA([1:m/2^k], [n/2^k + 1:n/2^(k-1)]);
        D{k} = WA([m/2^k + 1:m/2^(k-1)], [n/2^k + 1:n/2^(k-1)]);
        H{k} = WA([m/2^k + 1:m/2^(k-1)], [1:n/2^k]);
    end
    %ImagePlot(H{3})

    %Get q for V, D and H bands.  Send correct subbands, correct order!
    for k = 1 : numIts - 1
        qV{k} = vertPredictorNhbrs(V{k}, V{k+1}, D{k}, D{k+1})
        qD{k} = diagPredictorNhbrs(D{k}, D{k+1}, H{k}, V{k})
        qH{k} = horzPredictorNhbrs(H{k}, H{k+1}, D{k}, D{k+1})
    end
    disp('q vectors calculated')
    %size(qD{3})

    %Trim edges from subbands and flatten them
    for k = 1 : numIts - 1
        [m,n] = size(V{k});             %same as D{k}, H{k}
        V{k} = V{k}([3:m-2], [3:n-2]);  %remove edges
        V{k} = V{k}';                   %take transpose so entries will line up w/ q
        V{k} = V{k}(:);                 %and flatten 
        D{k} = D{k}([3:m-2], [3:n-2]);
        D{k} = D{k}';                   
        D{k} = D{k}(:);                  
        H{k} = H{k}([3:m-2], [3:n-2]);
        H{k} = H{k}';                   
        H{k} = H{k}(:);                  
    end
    disp('Subbands trimmed and flattened')

    %Calculate weight vectors for each band
    for k = 1 : numIts - 1
         wV{k} = (inv(qV{k}' * qV{k}))*(qV{k}' * V{k});
         wD{k} = (inv(qD{k}' * qD{k}))*(qD{k}' * D{k});
         wH{k} = (inv(qH{k}' * qH{k}))*(qH{k}' * H{k});
    end
    disp('Weight vectors calculated')
    %wV{3}
    %wH{3}

    %Find error between linear approx and actual pixel value
    for k = 1 : numIts - 1
        %Construct lin. approx.
        LV{k} = qV{k}*wV{k};
        LD{k} = qD{k}*wD{k};
        LH{k} = qH{k}*wH{k};
        %calculate errors
        errV{k} = LV{k} - V{k};
        errD{k} = LD{k} - D{k};
        errH{k} = LH{k} - H{k};
    end

    %Build signature vector
    tempSigVec = [];
    for k = 1 : numIts - 1
        tempSigVec = [tempSigVec mean(wV{k})];      %appends to sigVec
        tempSigVec = [tempSigVec var(wV{k})];
        tempSigVec = [tempSigVec skewness(wV{k})];
        tempSigVec = [tempSigVec kurtosis(wV{k})];
        tempSigVec = [tempSigVec mean(errV{k})];
        tempSigVec = [tempSigVec var(errV{k})];
        tempSigVec = [tempSigVec skewness(errV{k})];
        tempSigVec = [tempSigVec kurtosis(errV{k})];
        tempSigVec = [tempSigVec mean(wD{k})];
        tempSigVec = [tempSigVec var(wD{k})];
        tempSigVec = [tempSigVec skewness(wD{k})];
        tempSigVec = [tempSigVec kurtosis(wD{k})];
        tempSigVec = [tempSigVec mean(errD{k})];
        tempSigVec = [tempSigVec var(errD{k})];
        tempSigVec = [tempSigVec skewness(errD{k})];
        tempSigVec = [tempSigVec kurtosis(errD{k})];
        tempSigVec = [tempSigVec mean(wH{k})];
        tempSigVec = [tempSigVec var(wH{k})];
        tempSigVec = [tempSigVec skewness(wH{k})];
        tempSigVec = [tempSigVec kurtosis(wH{k})];
        tempSigVec = [tempSigVec mean(errH{k})];
        tempSigVec = [tempSigVec var(errH{k})];
        tempSigVec = [tempSigVec skewness(errH{k})];
        tempSigVec = [tempSigVec kurtosis(errH{k})];
    end
    sigVec{picNum} = tempSigVec;
end

%Each sigVec lives in R^72. Get stats for 1st, 2nd,... entry for each sigVec
numCoords = 3*(numIts-1)*8    %3 bands, numIts-1 levels, 8 stats each
for i = 1 : numCoords    
    for j = 1 : numPics       
        z(j) = sigVec{j}(i);  %holds ith coord for each sigVec
    end
    coordMean(i) = mean(z);
    coordSTD(i) = std(z);
end

%Standardize each coord of the sigVec's to have mean = 0, std = 1
for i = 1 : numCoords
    for j = 1 : numPics
        sigVec{j}(i) = (sigVec{j}(i) - coordMean(i)) / coordSTD(i);  %holds ith coord for each sigVec
    end
end

%Test for mean = 0, STD = 1
% for i = 1 : numCoords    
%     for j = 1 : numPics       
%         z(j) = sigVec{j}(i);  %holds ith coord for each sigVec
%     end
%     coordMean(i) = mean(z);
%     coordSTD(i) = std(z);
% end
% coordMean
% coordSTD

%Create mutual-distance matrix
for i = 1 : numPics
    for j = 1 : numPics
        mutualDist(i,j)= norm(sigVec{i}-sigVec{j});
    end
end
mutualDist

%Now process out of sample pics
outOfSampPics={'A4.png','B4.png','G4.png','C4.png','M4.png','T4.png'};
[rows,numOutOfSampPics]=size(outOfSampPics);

for picNum = 1 : numOutOfSampPics
    %Read in an image
    A=ImageRead(outOfSampPics{picNum});

    %Perform numIts iterations of wavelet transform
    WA = WT2D(A, h, numIts);
    disp('Completed wavelet transform')

    %Get dimensions of transformed pic (should be same as original)
    [m,n]=size(WA);

    %Get subbands. Store in cell arrays.containers
    for k = 1 : numIts
        V{k} = WA([1:m/2^k], [n/2^k + 1:n/2^(k-1)]);
        D{k} = WA([m/2^k + 1:m/2^(k-1)], [n/2^k + 1:n/2^(k-1)]);
        H{k} = WA([m/2^k + 1:m/2^(k-1)], [1:n/2^k]);
    end
    %ImagePlot(H{3})

    %Get q for V, D and H bands.  Send correct subbands, correct order!
    for k = 1 : numIts - 1
        qV{k} = vertPredictorNhbrs(V{k}, V{k+1}, D{k}, D{k+1});
        qD{k} = diagPredictorNhbrs(D{k}, D{k+1}, H{k}, V{k});
        qH{k} = horzPredictorNhbrs(H{k}, H{k+1}, D{k}, D{k+1});
    end
    disp('q vectors calculated')
    %size(qD{3})

    %Trim edges from subbands and flatten them
    for k = 1 : numIts - 1
        [m,n] = size(V{k});             %same as D{k}, H{k}
        V{k} = V{k}([3:m-2], [3:n-2]);  %remove edges
        V{k} = V{k}';                   %take transpose so entries will line up w/ q
        V{k} = V{k}(:);                 %and flatten 
        D{k} = D{k}([3:m-2], [3:n-2]);
        D{k} = D{k}';                   
        D{k} = D{k}(:);                  
        H{k} = H{k}([3:m-2], [3:n-2]);
        H{k} = H{k}';                   
        H{k} = H{k}(:);                  
    end
    disp('Subbands trimmed and flattened')

    %Calculate weight vectors for each band
    for k = 1 : numIts - 1
         wV{k} = (inv(qV{k}' * qV{k}))*(qV{k}' * V{k});
         wD{k} = (inv(qD{k}' * qD{k}))*(qD{k}' * D{k});
         wH{k} = (inv(qH{k}' * qH{k}))*(qH{k}' * H{k});
    end
    disp('Weight vectors calculated')
    %wV{3}
    %wH{3}

    %Find error between linear approx and actual pixel value
    for k = 1 : numIts - 1
        %Construct lin. approx.
        LV{k} = qV{k}*wV{k};
        LD{k} = qD{k}*wD{k};
        LH{k} = qH{k}*wH{k};
        %calculate errors
        errV{k} = LV{k} - V{k};
        errD{k} = LD{k} - D{k};
        errH{k} = LH{k} - H{k};
    end

    %Build signature vector
    tempSigVec = [];
    for k = 1 : numIts - 1
        tempSigVec = [tempSigVec mean(wV{k})];      %appends to sigVec
        tempSigVec = [tempSigVec var(wV{k})];
        tempSigVec = [tempSigVec skewness(wV{k})];
        tempSigVec = [tempSigVec kurtosis(wV{k})];
        tempSigVec = [tempSigVec mean(errV{k})];
        tempSigVec = [tempSigVec var(errV{k})];
        tempSigVec = [tempSigVec skewness(errV{k})];
        tempSigVec = [tempSigVec kurtosis(errV{k})];
        tempSigVec = [tempSigVec mean(wD{k})];
        tempSigVec = [tempSigVec var(wD{k})];
        tempSigVec = [tempSigVec skewness(wD{k})];
        tempSigVec = [tempSigVec kurtosis(wD{k})];
        tempSigVec = [tempSigVec mean(errD{k})];
        tempSigVec = [tempSigVec var(errD{k})];
        tempSigVec = [tempSigVec skewness(errD{k})];
        tempSigVec = [tempSigVec kurtosis(errD{k})];
        tempSigVec = [tempSigVec mean(wH{k})];
        tempSigVec = [tempSigVec var(wH{k})];
        tempSigVec = [tempSigVec skewness(wH{k})];
        tempSigVec = [tempSigVec kurtosis(wH{k})];
        tempSigVec = [tempSigVec mean(errH{k})];
        tempSigVec = [tempSigVec var(errH{k})];
        tempSigVec = [tempSigVec skewness(errH{k})];
        tempSigVec = [tempSigVec kurtosis(errH{k})];
    end
    outOfSampSigVec{picNum} = tempSigVec;
end

%Standardize each coord of the sigVec's to have mean = 0, std = 1
for i = 1 : numCoords
    for j = 1 : numOutOfSampPics
        outOfSampSigVec{j}(i) = (outOfSampSigVec{j}(i) - coordMean(i)) / coordSTD(i);  %holds ith coord for each sigVec
    end
end

%Find distances to in-sample pics
mutualDist = [];
for i = 1 : numPics
    for j = 1 : numOutOfSampPics
        mutualDist(i,j)= norm(sigVec{i}-outOfSampSigVec{j});
    end
end
mutualDist
