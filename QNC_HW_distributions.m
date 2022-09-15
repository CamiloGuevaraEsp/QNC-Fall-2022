%% Simulation of distribution shown in Fig 1(A) from:

%Fischer D, Lombardi DA, MarucciWellman H, Roenneberg T (2017) Chronotypes in
%the US – Influence of age and sex. PLoS ONE 12
%(6): e0178782. https://doi.org/10.1371/journal.pone.0178782

%This paper is the characterization of US chronotipes taking 
%data from the American Time Use Survey  (2003-2014). The chronotype (x
%axis) is based on mid-point of sleep on weekends from a population of
%53689 individuals 

%Use a random seed to get the same results 
rng('default')

%In the paper the mean chronotype is around 8 hours
mu = 8.5;
%To get the range (between ~3 and 14 hours I used the following sigma
sigma = 2;

%Sample size was
N = 53689;

% Get samples
samples =   normrnd(mu, sigma, N, 1);

% plot histogram, pick number of bins
cla reset; hold on;
%Pick 23 bins as in the paper
nbins = 23;
[counts, edges] = histcounts(samples, nbins);
xaxis = edges(1:end-1)+diff(edges);
npdf = counts./trapz(xaxis, counts);
bar(xaxis, npdf);

% Show theoretical pdf in red
plot(xaxis, normpdf(xaxis, mu, sigma), 'r-', 'LineWidth', 2);

% labels, ets
title('Distribution of sleep duration (hours)')
xlabel('Sleep duration');
ylabel('Probability');
legend('Simulated', 'Theoretical')