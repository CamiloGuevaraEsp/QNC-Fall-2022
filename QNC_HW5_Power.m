%Generate distributions for spiking data and pupil diameter 
%spikes = poissrnd(6);

effect_size = 0.1 : 0.1 : 2;

spikes = poissrnd(2,1000);
pupil = normrnd(0,0.5,1000);

%Calculate the correlation coefficient to produce the null distributiom
for i = 1:1000
   coefficient = corrcoef(spikes(i,:),pupil(i,:)); 
   X(i) = coefficient(2); 
end

%Plot the null distribution
nbins = 20;
[counts, edges] = histcounts(X, nbins);
xaxis = edges(1:end-1)+diff(edges);
npdf = counts./trapz(xaxis, counts);
bar(xaxis, npdf);
xlabel('Coefficient')
ylabel('Counts')

%Ussing the standar deviation and the mean calculate the mean value for a
%given effect size 
%effect = (mean_null - mean_treatment)/sd_null
%(effect * sd_null) - mean_null = -mean_treatment
%mean_null - (effect*sd_null) = mean_tratment

for  i = 1:length(effect_size)
    mean_treatment = mean(X) - (effect_size(i) * std(X));
    n(i) = sampsizepwr('t', [mean(X) std(X)],mean_treatment,0.80);
end

%Plot
figure
plot(effect_size,n, '--o', 'LineWidth',2)    
title ('Effect size vs sample size')
xlabel('Effect Size')
ylabel('Sample size')


