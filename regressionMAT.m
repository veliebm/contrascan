function[OutMat] = regressionMAT(InMat)
    % regression
    % compute linear regression and subtracts regression from input vector
    
    OutMat = zeros(size(InMat));
    
    % make sure InMat has electrodes as rows
        
    dummy = size(InMat);

    X = [ones(1,dummy(2));1:dummy(2)]';
    
    
    for channel = 1 : size(InMat,1)
        
        regkoefs = pinv(X' * X) * (X' * InMat(channel, :)');
        
        
        reggerade = [1:dummy(2)].*regkoefs(2,1) + regkoefs(1,1);
        OutMat(channel, :) = InMat(channel, :) - reggerade;
    end
end
