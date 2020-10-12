function J = costFunc(h,expY) 
    n = size(h,1) ; 
    xrand = [1, 1] ; 
    while length(unique(xrand)) ~= length(xrand)
        xrand = round( (expY(:,1) + 0.0005*rand([size(expY,1),1])), 4 ) ; % query points from experimental data with some random variation
    end
    
    yrand = [1, 1] ; 
    while length(unique(yrand)) ~= length(yrand)
        yrand = round( (expY(:,2) + 0.0005*rand([size(expY,1),1])), 4 ) ; % query points from experimental data with some random variation
    end

    xq = interp1(xrand,yrand,h(:,2)) ; 
    yq = interp1(xrand,yrand,h(:,1)) ; 
    
    xq(isnan(xq)) = 0 ; 
    yq(isnan(yq)) = 0 ; 

%     J = ( 1 / (2*n) ) * sum( (h(:,2)-yq).^2 ) ; 
    J = ( 1 / (2*n) ) * ( ([h(:,1) ; h(:,2)] - [xq ; yq])' * ([h(:,1) ; h(:,2)] - [xq ; yq]) ) ; 
%     J = ( 1 / (2*n) ) * ( (h-[xd ; yd])' * (h-[xd ; yd]) ) ; 
end
