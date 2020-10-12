function plotExpData(meanposx, meanposy, stdposx, stdposy, nfig, ntraj, wp, color)

    % calculate upper and lower error bounds on mean trajectory data
    errxu = meanposx + stdposx ; 
    errxl = meanposx - stdposx ; 
    erryu = meanposy + stdposy ;
    erryl = meanposy - stdposy ;

    % save vertices and faces for plotting
    v = {} ; 
    f = {} ; 
    
    figure(nfig) 
    
    for i = 1:ntraj
        v1 = [[errxl(wp(i):wp(i+1),1)', fliplr(errxl(wp(i):wp(i+1),1)')]',[erryu(wp(i):wp(i+1),1)',fliplr(erryl(wp(i):wp(i+1),1)')]'] ; % errors in y
        v2 = [[errxu(wp(i):wp(i+1),1)', fliplr(errxu(wp(i):wp(i+1),1)')]',[erryu(wp(i):wp(i+1),1)',fliplr(erryl(wp(i):wp(i+1),1)')]'] ; 
        v3 = [[errxl(wp(i):wp(i+1),1)', fliplr(errxu(wp(i):wp(i+1),1)')]',[erryl(wp(i):wp(i+1),1)',fliplr(erryl(wp(i):wp(i+1),1)')]'] ; % errors in x
        v4 = [[errxl(wp(i):wp(i+1),1)', fliplr(errxu(wp(i):wp(i+1),1)')]',[erryu(wp(i):wp(i+1),1)',fliplr(erryu(wp(i):wp(i+1),1)')]'] ; 
        vlist = [v1 , v2 , v3 , v4]  ; 

        f1 = linspace(1,length(v1),length(v1)) ; 
        f2 = linspace(1,length(v2),length(v2)) ; 
        f3 = linspace(1,length(v3),length(v3)) ; 
        f4 = linspace(1,length(v4),length(v4)) ; 
        flist = [f1 , f2, f3, f4] ; 
        
        patch('Faces',f1,'Vertices',v1,'FaceAlpha',1,'EdgeColor','none','FaceColor',color,'HandleVisibility','off') ; 
        hold on
        patch('Faces',f2,'Vertices',v2,'FaceAlpha',1,'EdgeColor','none','FaceColor',color,'HandleVisibility','off') ; 
        hold on
        patch('Faces',f3,'Vertices',v3,'FaceAlpha',1,'EdgeColor','none','FaceColor',color,'HandleVisibility','off') ; 
        hold on
        patch('Faces',f4,'Vertices',v4,'FaceAlpha',1,'EdgeColor','none','FaceColor',color,'HandleVisibility','off') ; 
        hold on
        
        v{i} = vlist ; 
        f{i} = flist ; 
    end
        
    figure(nfig)
    plot(meanposx,meanposy,'-ok') % plot average trajectory
    
end