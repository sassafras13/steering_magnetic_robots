function By = ByHH(x,y,mu0,nturns,I,a,E,K)
    By = ((mu0 * nturns * I) / (2 * pi)) * (x / (y * sqrt((y + a)^2 + x^2))) * (((a^2 + y^2 + x^2) / ((a - y)^2 + x^2)) * E - K) ; 
end