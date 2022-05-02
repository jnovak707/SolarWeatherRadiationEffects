

fit_params=[];          %at the end, we want a list of fit parameters for each shielding level so we can interpolate for a desired level
for ii=1:length(depths) %so for each shielding level
        radAccum=[];    %we're going to find how many krads were accumulated...    
        for i=1:length(storm)  %for each storm
            radAccum=horzcat(radAccum,storm(i).stormKrad(1,ii)); %so that we have an array of all the krad for each storm at the current shielding level
        end
        Const = polyfit(log10([storm(:).NOAA_PFU]),log10(radAccum),1);%fit the log-log space curve, so now we relate PFU to amount of radiation expected at that shielding level
        m = Const(1); fit_params(ii,1)=m;%save our fit parameters so we can interpolate with them later
        k = Const(2); fit_params(ii,2)=k;
        %extrapolated_doses = 10.^(m.*log10([1:50000])+(k)); %not neccessary but this is how we would check if the extrapolation was a good fit
end
fit_params=horzcat(fit_params,depths'); %now we have fit parameters that relate PFU to expected krads received at each shielding level based on the krad received from all storms at that shielding level

desired_shielding_m=interp1(fit_params(:,3),fit_params(:,2),350) %and we can use these to interpolate fit parameters for a user-specified shielding level
desired_shielding_k=interp1(fit_params(:,3),fit_params(:,1),350) %so when we stochastically generate a storm's PFU, we can translate that to the expected dose the spacecraft would receive with it's given shielding level (and orbit)

figure()
plot(depths,storm(1).stormKrad)
hold on
plot(depths,storm(2).stormKrad)
set(gca,'YScale','log')
set(gca,'XScale','log')