function plotSwimmer(nfig,y,yind,dirB)

    len = 0.005 ; 
    Bx1 = y(yind,1) - len*cos(y(yind,3)) ; 
    By1 = y(yind,2) - len*sin(y(yind,3)) ; 
    Ax1 = y(yind,1) + len*cos(y(yind,3)) ; 
    Ay1 = y(yind,2) + len*sin(y(yind,3)) ; 
    
    if dirB == 1 
        color1 = '#a6a8ab' ; 
        color2 = 'k' ; 
    elseif dirB == 0 
        color1 = 'k' ; 
        color2 = '#a6a8ab' ;
    end
    
    figure(nfig)
    plot([Bx1';y(yind,1)'],[By1';y(yind,2)'],'-','LineWidth',5,'Color',color1) ; 
    hold on
    plot([y(yind,1)';Ax1'],[y(yind,2)';Ay1'],'-','LineWidth',5,'Color',color2) ; 
    hold on

end