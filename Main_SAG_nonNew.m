%% Project: Shear-thickening SAG gas injection
%% Author: Rodrigo Salazar
%% August 2019

close all 
clear all 
clc

digits(64);
%% Simulation Parameters
num_chac=50; %number of characteristics
tdmax=1; %Dimensionless time
tdpoints=1000; %Grid refinement td
ffcurve_resolution=1000;
n_exp=0.34; %Non-Newtonian Exponent
number_rings=10; %
%% fluid properties
muw = 1e-3;
mug = 2e-5;

%% Radii of Borehole and Reservoir
rw = 0.1; %[m] wellbore radius 
re = 100; %[m] reservoir radius

%%Bentheimer Layer I 
%% Corey Parameters
krwe = 0.39; nw  = 2.86; krge = 0.59; ng  = 0.7; 
%% Residual Saturations
swr = 0.25; sgr = 0.2; 
%% Foam Model Parameters 
fmmob = 47700; fmdry = 0.271; epdry =400;



%fmdry discretization
r_vector=logspace(log10(rw),log10(re),number_rings); %logaritmic spacing

xD_vector(1:number_rings)=((r_vector(1:number_rings)).^2-rw^2)/(re^2-rw^2); %translated into xD

fmdry_vec(1:number_rings)=swr+(fmdry-swr)*(r_vector(1:number_rings)/re).^((n_exp-1)/nw); %fmdry discretization

%calculation shock using tangency condition in each ring
sw_shock_vector=[];

for i=1:number_rings
    [sw_shock_vector(i)]=shock_calc(swr,sgr,krwe,nw,krge,ng,mug,muw,fmmob,epdry,fmdry_vec(i));
end



fmdry=fmdry_vec(number_rings);

%initializing variables
n=ffcurve_resolution;
fwvalue=linspace(0.000001,1,ffcurve_resolution);
xsol=[];

%Definition of functions
swe  = @(sw)(sw-swr)/(1-swr-sgr); 
krw  = @(sw)(krwe*swe(sw).^nw); 
lambda_w = @(sw)(krw(sw)./muw); 
krg  = @(sw)(krge*(1-swe(sw)).^ng);
lambda_g = @(sw)(krg(sw)./mug); 
%NZ
FM   = @(sw)(1+fmmob*((0.5+ atan(epdry.*(sw-fmdry))/pi())-(0.5+ atan(epdry.*(swr-fmdry))/pi()))); 
krgf = @(sw)(krg(sw)./FM(sw)); 
lambda_f = @(sw)(real(krgf(sw)./mug)); 
fw = @(sw)(1./(1+(lambda_f(sw))./(lambda_w(sw)))); 
dkrw = @(sw)((nw*krwe*swe(sw).^(nw-1))./(1-swr-sgr)); 
dlambda_w = @(sw)(real(dkrw(sw)./muw)); 
dkrg = @(sw)(-(krge*ng*(1-swe(sw)).^(ng-1))./(1-swr-sgr)); 
dlambda_g = @(sw)(dkrg(sw)./mug); 
dFM  = @(sw)((fmmob*epdry)./(pi*(1+(epdry^2*(sw-fmdry).^2)))); 
dkrgf = @(sw)((dkrg(sw).*FM(sw)-dFM(sw).*krg(sw))./FM(sw).^2);
dlambda_f = @(sw)(dkrgf(sw)./mug); 
dfw = @(sw)((lambda_f(sw).*dlambda_w(sw) - dlambda_f(sw).*lambda_w(sw))./(lambda_f(sw)+lambda_w(sw)).^2); 
lambda_rt =@(sw)(lambda_w(sw)+lambda_f(sw)); 

fw_shock_rings(1:number_rings)=fw(sw_shock_vector(1:number_rings));
sw_shock=sw_shock_vector(number_rings);

%sw and fw values for the characterictics at the outer radius.
m = linspace(swr,sw_shock_vector(number_rings),num_chac+1);
N = length (m);
%fw vector for all fractional flow curves in each ring (including repeated
%characteristics)

fw_vector(1:num_chac+1)=fw(m(1:num_chac+1)); 


%saturations calculated at each ring
sw_ring=zeros(num_chac+1,number_rings);
for j=1:number_rings-1
    for i=1:num_chac+1
       sw_ring(i,j) = sw_from_fw(fw_vector(i),swr,sgr,krwe,nw,krge,ng,mug,muw,fmmob,epdry,fmdry_vec(j));
    end
end


i=1;
tD=linspace(0,tdmax,tdpoints);

figure(4)
plot(m,fw_vector,sw_ring(:,1),fw_vector,sw_ring(:,2),fw_vector)

%arranging saturation in a vector for final calculation
sw_ring_final=zeros(num_chac+1,number_rings);

%last layer assigment
sw_ring_final(:,number_rings)=transpose(m);

%remaining of layers assigment without characteristics faster than the
%shock

for j=1:number_rings-1
    for i=1:num_chac+1
        if fw_vector(i)<=fw_shock_rings(j)
            sw_ring_final(i,j)=sw_ring(i,j);
        else
            sw_ring_final(i,j)=sw_shock_vector(j);
            %sw_ring_final(i,j)=swr;
        end
    end
end

%reseting counting variables
j=1;
i=1;

%Calculating tD from xD 

td_plot_vector=zeros(num_chac+1,number_rings);
for j=1:number_rings-1
    for i=1:num_chac+1
        td_plot_vector(i,j+1)=td_plot(sw_ring_final(i,j),td_plot_vector(i,j),xD_vector(j),xD_vector(j+1),swr,sgr,krwe,nw,krge,ng,mug,muw,fmmob,epdry,fmdry_vec(j));
    end
end

i=1;
j=1;

%Ploting xD vs tD diagram
figure(3)
hold
for j=1:number_rings-1  
    for i=1:N
        plot([td_plot_vector(i,j);td_plot_vector(i,j+1)],[xD_vector(j);xD_vector(j+1)],'Color','Blue')
    end
end
axis([0 1 0 1])
i=1;
j=1;

%Plot mobility leading characteristic
for i=1:number_rings
    lambda_rt_vector(i)=lambda_behind_shock(sw_ring_final(num_chac+1,i),swr,sgr,krwe,nw,krge,ng,mug,muw,fmmob,epdry,fmdry_vec(i));
end
i=1;

figure (5)
plot(td_plot_vector(num_chac+1,1:number_rings),lambda_rt_vector)


%Deleting repeated characteristics


%Calculating non-Newtonian Pressure

for k=1:tdpoints-1
    
    %Calculating mobilities

    td_pressure_plot=tD(k+1);
    lambda_chac_vector=[];
    r_chac_vector=[];
    lambda_chac=0;
    for j=1:number_rings-1
        for i=1:num_chac+1
            if (td_plot_vector(i,j)<td_pressure_plot) && (td_plot_vector(i,j+1)>td_pressure_plot)
                lambda_chac=lambda_behind_shock(sw_ring_final(i,j),swr,sgr,krwe,nw,krge,ng,mug,muw,fmmob,epdry,fmdry_vec(j));
                lambda_chac_vector=[lambda_chac_vector,lambda_chac];
                r_chac_value=r_chac(td_pressure_plot,sw_ring_final(i,j),swr,sgr,krwe,nw,krge,ng,mug,muw,fmmob,epdry,fmdry_vec(j));
                r_chac_vector=[r_chac_vector,r_chac_value];
            end
        end
    end
    i=1;
    j=1;

    %Deleting mobilities and radii of repetead characteristics

    lambda_chac_pressure_vector=[];
    r_chac_pressure_vector=[];

    for i=1:length(lambda_chac_vector)
        if lambda_chac_vector(i)~=lambda_chac_vector(length(lambda_chac_vector))
            lambda_chac_pressure_vector=[lambda_chac_pressure_vector,lambda_chac_vector(i)];
        end
        if r_chac_vector(i)~=r_chac_vector(length(lambda_chac_vector))
            r_chac_pressure_vector=[r_chac_pressure_vector,r_chac_vector(i)];
        end
    end
    lambda_chac_pressure_vector=[lambda_chac_pressure_vector,lambda_chac_vector(length(lambda_chac_vector))];
    r_chac_pressure_vector=[r_chac_pressure_vector,r_chac_vector(length(lambda_chac_vector))];
    i=1;
    j=1;

    % Calculating Dimensionless Pressure Point
    N=length(r_chac_pressure_vector);
    P_nn(i:N-1) = log(r_chac_pressure_vector(i+1:N)./r_chac_pressure_vector(i:N-1)).*(0.5*(1./lambda_chac_pressure_vector(i+1:N)+1./lambda_chac_pressure_vector(i:N-1)));
    P_nn(N)     = log(re./r_chac_pressure_vector(N))*muw;
    PD_nn(k) = sum(P_nn)/(muw*log(re/rw));
    lambda_chac_vector=[];
    r_chac_vector=[];
    lambda_chac_pressure_vector=[];
end
figure(6)
plot(tD(1:tdpoints-1),PD_nn)
max(PD_nn)