function quadrant = checkQuad(x,y)
    if ((x > 0) && (y > 0))
        quadrant = 1 ; 
    elseif ((x < 0) && (y > 0))
        quadrant = 2 ; 
    elseif ((x < 0) && (y < 0))
        quadrant = 3 ; 
    elseif ((x > 0) && (y < 0))
        quadrant = 4 ; 
    else
        quadrant = 1 ; 
    end
end