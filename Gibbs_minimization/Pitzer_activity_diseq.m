%%% The function Pitzer_activity_diseq calculates the activity coefficients
%%% of aqueous species and liquid water activity using a simplified version
%%% of the Pitzer equations. 
%%% This function takes as inputs the names of all species, their
%%% abundances, charges, and phases.

function [lnPhiWATER,act_an,act_ct,a_index,c_index]=Pitzer_activity_diseq(names_in,n_in,v_in,l)
    global B0 B1 B2 C0 om scale_factor % Ensure Pitzer equation binary interaction parameters are loaded from database file
    
    n_in=n_in';

    mass_ocean=om*1.34e9*1000^3*1030; % mass of ocean in kg 
    moles_atm=1.7620e+20; % total moles in atmosphere
    n=n_in.*moles_atm./mass_ocean./scale_factor; %convert normalized moles back to molalities
    
    %%% Firstly, reduce name, abundance and charge vectors to aqueous species
    %%% only (select for l=4 phase)
    [row,ve] = find(l==4);
    names=names_in(row);
    names=strrep(names,')',''); % remove brackets from species names
    names=strrep(names,'(',''); % remove brackets from species names
    names = strtrim(names);
    v=v_in(row);
    m=n(row);

    % Find cations
    [row,c_index]=find(v_in>0); % find indices for all postively charged aqueous species
    names_c=names_in(c_index);
    names_c=strrep(names_c,')','');
    names_c=strrep(names_c,'(','');
    names_c = strtrim(names_c); % vector containing cation names
    v_c=v_in(c_index); % vector containing cation charges
    mc=n(c_index); % vector containing cation abundances

    %anions
    [row,a_index]=find(v_in<0); % find indices for all negatively charged aqueous species
    names_a=names_in(a_index);
    names_a=strrep(names_a,')','');
    names_a=strrep(names_a,'(','');
    names_a = strtrim(names_a); % vector containing anion names
    v_a=v_in(a_index); % vector containing cation charges
    ma=n(a_index); % vector containing cation abundances
 
    %%%% The terms below are defined in equation 11 in Krissansen-Totton et
    %%%% al. (2018), originally from Marion and Kargel (2007).
    sum_mi=sum(m);
    I=0.5*sum(m.*v.^2);
    A_phi=0.3915; %This could be modified to include temperature dependence, but is constant here.
    b=1.2; % Constant (see Marion and Kargel (2007)).
    Z=sum(m.*abs(v));
 
    f_gamma=-A_phi*(I.^0.5/(1+b*I.^0.5)+2*log(1+b*I.^0.5)/b); % part of equation 11 in Krissansen-Totton et al. (2018)

     ln_gamma_ct=0*mc; % Create empty arrays for log activity coefficients
     ln_gamma_an=0*ma; % Create empty arrays for log activity coefficients
     F=f_gamma;
    
     % In the portions of code that follow, B0, B1, B2, and C0 will be
     % treated separately. Note that B0, B1, and B2 are abbreviated forms
     % of the terms B(0)_MX, B(1)_MX, and B(2)_MX, referred to in KT16 and
     % KT18. Bracketed numbers should be superscripted. C0 is a special 
     % case (see below).
     
     % Begin B0 section
    ind1=find(ismember(B0(:,1),names_c)); % find all cations in first column in database
    ind2=find(ismember(B0(:,2),names_a)); % find all anions in second column in database
    [C,ia,ib] = intersect(ind1,ind2); % % C contains indices of all cation-anion database rows relevant to binary interactions in this system
    B0_used=[B0(C,1) B0(C,2) B0(C,3) B0(C,3) B0(C,3)]; % Fill B0_used with relevant rows above: cation, anion, parameter, parameter, parameter (last two columns will later be filled with abundances)
    ind1=find(ismember(B0(:,1),names_a)); % find all anions in first column in database
    ind2=find(ismember(B0(:,2),names_c)); % find all cations in second column in database
    [C,ia,ib] = intersect(ind1,ind2); % C contains indices of all anion-cation database rows relevant to binary interactions in this system
    B0_used=[B0_used; B0(C,1) B0(C,2) B0(C,3) B0(C,3) B0(C,3)]; % Fill in B0_used with relevant rows above: anion, cation, parameter, parameter, parameter (last two columns will later be filled with abundances)
    IndexQ = strfind(B0_used, '+');
    B0_sum=0;
    for i=1:length(B0_used(:,1)) % Here we calculate the 2*B0*abundance contribution to the terms in equation 11 in Krissansen-Totton et al. (2018).
        if cell2mat(IndexQ(i,1))>0 %cation-anion contributions
            mi1=find(ismember(names_c,B0_used(i,1)));
            mi2=find(ismember(names_a,B0_used(i,2)));
            B0_used(i,4)={mc(mi1)}; %cation abundance
            B0_used(i,5)={ma(mi2)}; %anion abundance
            ln_gamma_ct(mi1)=ln_gamma_ct(mi1)+2*cell2mat(B0_used(i,5))*eval(cell2mat(B0_used(i,3)));
            ln_gamma_an(mi2)=ln_gamma_an(mi2)+2*cell2mat(B0_used(i,4))*eval(cell2mat(B0_used(i,3)));
        else %anion-cation contribtions
            mi1=find(ismember(names_a,B0_used(i,1)));
            mi2=find(ismember(names_c,B0_used(i,2)));
            B0_used(i,4)={ma(mi1)}; %anion abundance
            B0_used(i,5)={mc(mi2)}; %cation abundance
            ln_gamma_an(mi1)=ln_gamma_an(mi1)+2*cell2mat(B0_used(i,5))*eval(cell2mat(B0_used(i,3)));
            ln_gamma_ct(mi2)=ln_gamma_ct(mi2)+2*cell2mat(B0_used(i,4))*eval(cell2mat(B0_used(i,3)));
        end
        % Next we calculate the B0*abundance_cation*abundance_anion
        % contribution to equation 34 in Krissansen-Totton et al. (2016)
        B0_sum=B0_sum+eval(cell2mat(B0_used(i,3))).*cell2mat(B0_used(i,5)).*cell2mat(B0_used(i,4));
    end

    % Begin B1 section (same general structure as B0)
    ind1=find(ismember(B1(:,1),names_c));
    ind2=find(ismember(B1(:,2),names_a));
    [C,ia,ib] = intersect(ind1,ind2);
    B1_used=[B1(C,1) B1(C,2) B1(C,3) B1(C,3) B1(C,3) B1(C,3)];
    ind1=find(ismember(B1(:,1),names_a));
    ind2=find(ismember(B1(:,2),names_c));
    [C,ia,ib] = intersect(ind1,ind2);
    B1_used=[B1_used; B1(C,1) B1(C,2) B1(C,3) B1(C,3) B1(C,3) B1(C,3)]; 
    IndexC1 = strfind(B1_used, '+2'); % Different parameterizations for 2:2 electrolytes
    IndexC2 = strfind(B1_used, '-2'); % Different parameterizations for 2:2 electrolytes
    IndexQ = strfind(B1_used, '+');
    B1_sum=0;
    for i=1:length(B1_used(:,1))
        if cell2mat(IndexQ(i,1))>0
            mi1=find(ismember(names_c,B1_used(i,1)));
            mi2=find(ismember(names_a,B1_used(i,2)));
            B1_used(i,4)={mc(mi1)}; % fill abundances
            B1_used(i,5)={ma(mi2)}; % fill abundances
        else
            mi1=find(ismember(names_a,B1_used(i,1)));
            mi2=find(ismember(names_c,B1_used(i,2)));
            B1_used(i,4)={ma(mi1)}; % fill abundances
            B1_used(i,5)={mc(mi2)}; % fill abundances
        end

        % define alpha parameters depending on whether 2:2 electrolytes
        % (see Materials and Methods, Krissansen-Totton et al. (2018)).
        if and(cell2mat(IndexC1(i,1))>0,cell2mat(IndexC2(i,2))>0)
            B1_used(i,6)={[1.4]}; % fill alp[ha parameter
        elseif and(cell2mat(IndexC2(i,1))>0,cell2mat(IndexC1(i,2))>0) 
            B1_used(i,6)={[1.4]}; % fill alp[ha parameter
        else
             B1_used(i,6)={[2.0]}; % fill alp[ha parameter
        end
        
        % Compute 2*B1*f()*abundance contribution in equation 11 in Krissansen-Totton
        % et al. (2018).
        if cell2mat(IndexQ(i,1))>0
            ln_gamma_ct(mi1)=ln_gamma_ct(mi1)+2*g_fun(cell2mat(B1_used(i,6))*sqrt(I))*cell2mat(B1_used(i,5))*eval(cell2mat(B1_used(i,3)));    
            ln_gamma_an(mi2)=ln_gamma_an(mi2)+2*g_fun(cell2mat(B1_used(i,6))*sqrt(I))*cell2mat(B1_used(i,4))*eval(cell2mat(B1_used(i,3)));            
        else
            ln_gamma_an(mi1)=ln_gamma_an(mi1)+2*g_fun(cell2mat(B1_used(i,6))*sqrt(I))*cell2mat(B1_used(i,5))*eval(cell2mat(B1_used(i,3)));
            ln_gamma_ct(mi2)=ln_gamma_ct(mi2)+2*g_fun(cell2mat(B1_used(i,6))*sqrt(I))*cell2mat(B1_used(i,4))*eval(cell2mat(B1_used(i,3)));
        end
        
        % Compute abunance_cation*abundance_anion*B1*f'()/I contribution to
        % equation 11 in Krissansen-Totton et al. (2018)
        F=F+cell2mat(B1_used(i,5)).*cell2mat(B1_used(i,4))*eval(cell2mat(B1_used(i,3)))*g_prime(cell2mat(B1_used(i,6))*I.^0.5)./I;

        % Compute contribution to B1*abundance_cation*abundance_anion term
        % in equation 34 in Krissansen-Totton et al. (2016)
        g1=exp(-cell2mat(B1_used(i,6)).*sqrt(I));
        B1_sum=B1_sum+eval(cell2mat(B1_used(i,3))).*g1.*cell2mat(B1_used(i,5)).*cell2mat(B1_used(i,4));
    end
 
    % Begin B2 section (same general structure as B1)
    ind1=find(ismember(B2(:,1),names_c));
    ind2=find(ismember(B2(:,2),names_a));
    [C,ia,ib] = intersect(ind1,ind2);
    B2_used=[B2(C,1) B2(C,2) B2(C,3) B2(C,3) B2(C,3) B2(C,3)]; 
    ind1=find(ismember(B2(:,1),names_a));
    ind2=find(ismember(B2(:,2),names_c));
    [C,ia,ib] = intersect(ind1,ind2);
    B2_used=[B2_used; B2(C,1) B2(C,2) B2(C,3) B2(C,3) B2(C,3) B2(C,3)]; 
    IndexC1 = strfind(B2_used, '+2');
    IndexC2 = strfind(B2_used, '-2');
    IndexQ = strfind(B2_used, '+');
    B2_sum=0;
    for i=1:length(B2_used(:,1))
        if cell2mat(IndexQ(i,1))>0
            mi1=find(ismember(names_c,B2_used(i,1)));
            mi2=find(ismember(names_a,B2_used(i,2)));
            B2_used(i,4)={mc(mi1)};
            B2_used(i,5)={ma(mi2)};
        else
            mi1=find(ismember(names_a,B2_used(i,1)));
            mi2=find(ismember(names_c,B2_used(i,2)));
            B2_used(i,4)={ma(mi1)};
            B2_used(i,5)={mc(mi2)};
        end

        if and(cell2mat(IndexC1(i,1))>0,cell2mat(IndexC2(i,2))>0)
            B2_used(i,6)={12.0};
        elseif and(cell2mat(IndexC2(i,1))>0,cell2mat(IndexC1(i,2))>0) 
             B2_used(i,6)={12.0};
        else
            B2_used(i,6)={0.0};
        end
         
        % equation 11 contributions KT18
        if cell2mat(IndexQ(i,1))>0
            ln_gamma_ct(mi1)=ln_gamma_ct(mi1)+2*g_fun(cell2mat(B2_used(i,6))*sqrt(I))*cell2mat(B2_used(i,5))*eval(cell2mat(B2_used(i,3)));    
            ln_gamma_an(mi2)=ln_gamma_an(mi2)+2*g_fun(cell2mat(B2_used(i,6))*sqrt(I))*cell2mat(B2_used(i,4))*eval(cell2mat(B2_used(i,3)));    
        else
            ln_gamma_an(mi1)=ln_gamma_an(mi1)+2*g_fun(cell2mat(B2_used(i,6))*sqrt(I))*cell2mat(B2_used(i,5))*eval(cell2mat(B2_used(i,3)));
            ln_gamma_ct(mi2)=ln_gamma_ct(mi2)+2*g_fun(cell2mat(B2_used(i,6))*sqrt(I))*cell2mat(B2_used(i,4))*eval(cell2mat(B2_used(i,3)));
        end

        % equation 34 contributions KT16
        g2=exp(-cell2mat(B2_used(i,6)).*sqrt(I));
        B2_sum=B2_sum+eval(cell2mat(B2_used(i,3))).*g2.*cell2mat(B2_used(i,5)).*cell2mat(B2_used(i,4));
        
        % equation 11 contributions KT18
        F=F+cell2mat(B2_used(i,5)).*cell2mat(B2_used(i,4))*eval(cell2mat(B2_used(i,3)))*g_prime(cell2mat(B2_used(i,6))*I.^0.5)./I;

    end


    % Begin C0 section
    ind1=find(ismember(C0(:,1),names_c));
    ind2=find(ismember(C0(:,2),names_a));
    [C,ia,ib] = intersect(ind1,ind2);
    C0_used=[C0(C,1) C0(C,2) C0(C,3) C0(C,3) C0(C,3)];
    ind1=find(ismember(C0(:,1),names_a));
    ind2=find(ismember(C0(:,2),names_c));
    [C,ia,ib] = intersect(ind1,ind2);
    C0_used=[C0_used; C0(C,1) C0(C,2) C0(C,3) C0(C,3) C0(C,3)];
    IndexQ = strfind(C0_used, '+');
    C0_sum=0;
    for i=1:length(C0_used(:,1))
        if cell2mat(IndexQ(i,1))>0
            mi1=find(ismember(names_c,C0_used(i,1)));
            mi2=find(ismember(names_a,C0_used(i,2)));
            C0_used(i,4)={mc(mi1)};
            C0_used(i,5)={ma(mi2)};
            C0_used(i,6)={v_c(mi1)}; % For C0 terms, cation charges are needed
            C0_used(i,7)={v_a(mi2)}; % For C0 terms, anion charges are needed
        else
            mi1=find(ismember(names_a,C0_used(i,1)));
            mi2=find(ismember(names_c,C0_used(i,2)));
            C0_used(i,4)={ma(mi1)};
            C0_used(i,5)={mc(mi2)};
            C0_used(i,6)={v_a(mi1)}; % anion charges
            C0_used(i,7)={v_c(mi2)}; % cation charges
        end
        
        % Contribution to equation 34 in KT16 - see Marion and Kargel
        % (2007) for full definiion of C0
        denom=(2*(abs(cell2mat(C0_used(i,6))*cell2mat(C0_used(i,7))))^0.5);
        C0_sum=C0_sum+cell2mat(C0_used(i,5)).*cell2mat(C0_used(i,4))*Z*eval(cell2mat(C0_used(i,3)))/denom;
       
        % Contribution to equation 11 in KT18
        if cell2mat(IndexQ(i,1))>0
            ln_gamma_ct(mi1)=ln_gamma_ct(mi1)+cell2mat(C0_used(i,5)).*Z*eval(cell2mat(C0_used(i,3)))/denom;    
            ln_gamma_an(mi2)=ln_gamma_an(mi2)+cell2mat(C0_used(i,4)).*Z*eval(cell2mat(C0_used(i,3)))/denom;    
        else
            ln_gamma_an(mi1)=ln_gamma_an(mi1)+cell2mat(C0_used(i,5)).*Z*eval(cell2mat(C0_used(i,3)))/denom;
            ln_gamma_ct(mi2)=ln_gamma_ct(mi2)+cell2mat(C0_used(i,4)).*Z*eval(cell2mat(C0_used(i,3)))/denom;
        end
    end
    
    summation=B0_sum+B1_sum+B2_sum+C0_sum; % double summation in equation 34 in Krissansen-Totton et al. (2016)

    ln_gamma_ct=ln_gamma_ct+v_c.^2*F; % Equation 11 in Krissansen-Totton et al. (2018)      
    ln_gamma_an=ln_gamma_an+v_a.^2*F; % Equation 11 in Krissansen-Totton et al. (2018)

    act_ct=exp(ln_gamma_ct); % return activity coefficients rather than their logarithms
    act_an=exp(ln_gamma_an); % return activity coefficients rather than their logarithms

    osmo=1+(2./sum_mi).*(-A_phi.*I^1.5/(1+b.*I.^0.5)+summation); %osmotic coefficient (see Krissansen-Totton et al. (2016))

    a_w=exp(-osmo*sum_mi/55.50844); % activity of water

    lnPhiWATER=log(a_w); % take natural log before returning water activity
end

function gp = g_prime(x) % This is f' in Marion and Kergel (2007) and Krissansen-Totton et al. (2018).
    if or(x>0,x<0)
        gp=-2*(1 - (1+x+x.^2/2)*exp(-x))/x.^2;
    else
        gp=0.0; %this is limit as x approaches 0
    end
end

function gp = g_fun(x)  % This is f in Marion and Kergel (2007) and Krissansen-Totton et al. (2018).
    if or(x>0,x<0)
        gp=2*(1 - (1+x)*exp(-x))/x.^2;
    else
        gp=1.0; %this is limit as x approaches 0
    end
end

