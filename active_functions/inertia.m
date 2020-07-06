function I12 = inertia(r1, r2, m1, m2)
    I1 = (2/5)*m1*(r1^2) ; 
    I2 = (2/5)*m2*(r2^2) + m2*((r1+r2)^2) ; 
    I12 = I1 + I2 ; 
end