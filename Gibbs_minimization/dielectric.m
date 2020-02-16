% Calculates the dielectric constant for water using the data from Owin 1961

function di = dielectric(T,P)
    t = T-273.15;
    p = P-1;
    a = [-22.5713 -.032066 -.00028568 .0011832 .000027895 -.00000001476 2300.64 -.13476];
    D0 = 4.476150;
    di = exp((-10^6*D0+2*a(1)*p+2*a(2)*p*t+2*a(3)*p*t^2+2*a(4)*p^2+2*a(5)*p^2*t+2*a(6)*p^2*t^2+2*a(7)*t+2*a(8)*t^2)/(-(10^6)));
end