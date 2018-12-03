function [snr,meanSpike,N,signal,noise] = getSNRPlaidNorm(segmentData)

maxLimit = 10000;

if isempty(segmentData)
    snr = 0;
    N = 0;
    meanSpike = [];
    signal = 0;
    noise  = 0;
    
else
    [~,N] = size(segmentData);
    
    if N > maxLimit
        segmentData(:,maxLimit:end)=[];
        N = maxLimit;
    end
    
    meanSpike = mean(segmentData,2);
    noise2D = segmentData - repmat(meanSpike,1,size(segmentData,2));
    
    signal = max(meanSpike) - min(meanSpike);
    noise = 2*std(noise2D(:));  % Kelley RC et al., 2007 JNS
    
    snr = signal/noise;
end
end