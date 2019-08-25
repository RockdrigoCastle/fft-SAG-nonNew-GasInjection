function [ r_sw ] = r_chac(t_d,sw,swr,sgr,krwe,nw,krge,ng,mug,muw,fmmob,epdry,fmdry)
    %% Radii of Borehole and Reservoir
    rw = 0.1; %[m] wellbore radius 
    re = 100; %[m] reservoir radius
    
    swe  = @(sw)(sw-swr)/(1-swr-sgr); 
    krw  = @(sw)(krwe*swe(sw).^nw); 
    lambda_w = @(sw)(krw(sw)./muw); 
    krg  = @(sw)(krge*(1-swe(sw)).^ng);
    lambda_g = @(sw)(krg(sw)./mug); 
    FM   = @(sw)(1+fmmob*((0.5+ atan(epdry.*(sw-fmdry))/pi())-(0.5+ atan(epdry.*(swr-fmdry))/pi()))); 
    krgf = @(sw)(krg(sw)./FM(sw)); 
    lambda_f = @(sw)(real(krgf(sw)./mug)); 
    fw   = @(sw)(1./(1+(lambda_f(sw))./(lambda_w(sw)))); 
    dkrw = @(sw)((nw*krwe*swe(sw).^(nw-1))./(1-swr-sgr)); 
    dlambda_w = @(sw)(real(dkrw(sw)./muw)); 
    dkrg = @(sw)(-(krge*ng*(1-swe(sw)).^(ng-1))./(1-swr-sgr)); 
    dlambda_g = @(sw)(dkrg(sw)./mug); 
    dFM  = @(sw)((fmmob*epdry)./(pi*(1+(epdry^2*(sw-fmdry).^2)))); 
    dkrgf = @(sw)((dkrg(sw).*FM(sw)-dFM(sw).*krg(sw))./FM(sw).^2);
    dlambda_f = @(sw)(dkrgf(sw)./mug); 
    dfw = @(sw)((lambda_f(sw).*dlambda_w(sw) - dlambda_f(sw).*lambda_w(sw))./(lambda_f(sw)+lambda_w(sw)).^2); 
    lambda_rt =@(sw)(lambda_w(sw)+lambda_f(sw));
    xD = @(sw)(min(t_d.*dfw(sw),1));
    r  = @(sw)(sqrt(xD(sw).*(re^2-rw^2)+rw^2));
    r_sw=r(sw);
end