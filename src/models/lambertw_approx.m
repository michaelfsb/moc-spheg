function result = lambertw_approx(x)
    % Lambert W function analytic approximation
    
    E = 0.4586887;
    result = (1+E)*log(6/5*x/log(12/5*x/log(1+12/5*x))) ...
        - E*log(2*x/log(1+2*x));
end
