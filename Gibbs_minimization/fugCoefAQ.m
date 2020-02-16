% This function uses the returns the value of the activity coefficient at
% the given temperature and pressure for the set of given aqueous species
% with the specified mole amounts. The Truesdell-Jones equation (described
% in Krissansen-Totton et al. (2016)) is used to calculate aqueous activity
% coefficients. This function has been superceded by the Pitzer equation
% method, which is more accurate for Earth-like systems.

function lnPhiAQ=fugCoefAQ(temperature,pressure,names,n,v)
     
    global om scale_factor

    mass_ocean=om*1.34e9*1000^3*1030; % mass of ocean in kg
    moles_atm=1.7620e+20; % total moles in atmosphere
    n=n.*moles_atm./mass_ocean; %convert normalized moles back to molalities
    
    I=0.5*sum(n.*v.*v./scale_factor);
    A=0.5092; %valid at 25 deg C
    B=0.3283; %valid at 25 deg C
    
    a_array=[];
    b_array=[];
    for i=1:length(v)
        if abs(v(i))>0
            names(i);
            [index1,index2] = searchDataD(names(i));
            a_array=[a_array,index1];
            b_array=[b_array,index2];
         elseif abs(v(i))==0
            a_array=[a_array,0];
            b_array=[b_array,0];
       end
    end
    a_array=a_array';
    b_array=b_array';
    lnPhiAQ=-A*v.^2.*(sqrt(I)./(1+B.*a_array.*sqrt(I)))+b_array.*I;
    
end

