clear all
%% Exercise 1

n = 10;
p = 0.2; 
k = 0:10; 

fprintf('Exercise 1\n');
for i = 1: length(k)
    p_release(i) = binopdf(k(i),n,p);
    fprintf ('The probability of observing %d quantal releases is %e\n', k(i), p_release(i))
end

%% Exercise 2.1
clear all
n = 14;
p = 0.1;
k = 8;

p_release_eight = binopdf(k,n,p);
fprintf('Exercise 2.1\n');
fprintf ('The probability of observing %d quantal releases with p = %2.1f is %e\n', k, p, p_release_eight)

%% Exercise 2.2 

n = 14;
p = 0.7;
k = 8;

p_release_eight = binopdf(k,n,p);
fprintf('Exercise 2.2\n');
fprintf ('The probability of observing %d quantal releases with p = %2.1f is %2.3f\n', k, p, p_release_eight) 
%% Exercise 2.3
p = 0.1:0.1:1;

for i = 1:length(p)
    binomial_coeff = factorial(n) / (factorial(k) * factorial(n-k));
    eight_likelihood(i) = binomial_coeff * p(i)^k * ((1-p(i)) ^ (n-k));
    
end

%Plot the likelihood function
figure
plot(p, eight_likelihood, '-o'); 
xlabel('Release probability')
ylabel('Likelihood')
title('Exercise 2.3')

fprintf('Exercise 2.3\n');
fprintf ('The most likely release probability to observe release of %d quanta is %3.1f\n', k, p(eight_likelihood == max(eight_likelihood))); 

%% Exercise 3.1

%assuming true release probability of p = 0.1
n=14;
p = 0.1;
k=5;

%binomial_coeff = factorial(n) / (factorial(k) * factorial(n-k));
%p_release_five = binomial_coeff * p^k * ((1-p) ^ (n-k)); 
p_release_five = binopdf(k,n,p);

total_likelihood = p_release_five * eight_likelihood(1); 
total_log_likelihood = log(p_release_five) + log(eight_likelihood(1));

fprintf('Exercise 3.1\n')
fprintf ('For a true release probability of %2.1f\n', p)
fprintf ('Total Likelihood = %e\n', total_likelihood)
fprintf('Total log likelihood = %3.1f\n', total_log_likelihood)

%% Exercise 3.2
%Now compute the likelihood function for five releases

p = 0.1:0.1:1;

for i = 1:length(p)
    five_likelihood(i) = binopdf(k,n,p(i));
end

%Plot the likelihood function
figure
plot(p, five_likelihood,'-o');
xlabel('Release probability')
ylabel('Likelihood')
title('Exercise 3.2')

fprintf('Exercise 3.2\n')
fprintf ('The most probable release probability to observe release of %d quanta is %3.1f\n', k, p(five_likelihood == max(five_likelihood))); 

% Calculate the total likelihood for each probability value in both
% experiments

total_likelihood = eight_likelihood .* five_likelihood;
total_log_likelihood = log(five_likelihood) + log (eight_likelihood); 


fprintf ('Maximum Likelihood considering both experiments = %3.3f\n', max(total_likelihood))
fprintf('Maximun log likelihood considering both experiments = %3.3f\n', max(total_log_likelihood))
fprintf('release probability with highest likelihood = %3.3f\n', p(total_likelihood == max(total_likelihood)))

%% Exercise 3.3
% Increase resolution of probabilties
clear all 
n = 14;
k = [5 8];
delta = [0.1 0.01 0.001 0.0001 0.00001 0.000001];
fprintf('Exercise 3.3\n') 
for i = 1:length(delta)
   p = 0.1: delta(i):1;
   for j = 1:length(k)
    binomial_coeff(j) = factorial(n) / (factorial(k(j)) * factorial(n-k(j)));
   end
    for ii = 1:length(p)
    five_likelihood(ii) = binomial_coeff(1) * p(ii)^k(1) * ((1-p(ii)) ^ (n-k(1)));
    eight_likelihood(ii) = binomial_coeff(2) * p(ii)^k(2) * ((1-p(ii)) ^ (n-k(2)));
    end
    total_likelihood = eight_likelihood .* five_likelihood;
    total_log_likelihood = log(five_likelihood) + log (eight_likelihood); 
    probabilities(i) = p(total_likelihood == max(total_likelihood));
    fprintf('release probability with highest likelihood at %2.6f probability steps = %3.6f\n', delta(i), probabilities(i))
end
figure
plot(log(delta),probabilities,'-o')
xlabel('log of probability resolution')
ylabel('Maximum Likelihood')
title('Exercise 3.3')

%% Exercise 3.4
clear all
%Increase sample size, to draw  possible outcomes, based on the result from 3.2 i will assume that the
%true release probability is 0.46

p = 0.1;
n = 14;       % number of "trials" per "experiment"
NumExperiments = 1:10:100;
p = 0.1 : 0.001 : 1;

fprintf('Exercise 3.4\n')
for i = 1:length(NumExperiments) 
    outcomes = binornd(n,0.46,NumExperiments(i),1);
    for j = 1: length(p)
        test(:,j) = binopdf(outcomes,n,p(j));
    end
    test_likelihood = prod(test);
    probabilities(i) = p(test_likelihood == max(test_likelihood));
    fprintf('release probability with highest likelihood  with %d samples = %3.6f\n', NumExperiments(i), probabilities(i))
    clear test 
end
figure
plot(NumExperiments,probabilities,'-o')
xlabel('Number of samples')
ylabel('Maximum Likelihood')
title('Exercise 3.4')

%% Exercise 4

clear all
%Write the experimental data as a vector
data = [2*ones(1,3) 4*ones(1,10) 5*ones(1,19) 6*ones(1,26) 7*ones(1,16) 8 *ones(1,16) 9*ones(1,5) 10*ones(1,5)];

%Stimate p-hat 
phat = mle(data,'Distribution','binomial','NTrials',14)
fprintf('Exercise 4\n')
fprintf('p-hat = %2.3f\n',phat)

%Plot histogram
edges = -0.5:14.5;
counts = histcounts(data, edges);

% Show a bar plot of the simulated bionimal distribution
clf;
xs = edges(1:end-1)+diff(edges)/2;
bar(xs, counts);

ylabel('Count');

% Normalize it to make it a pdf. Here counts (the x-axis of the histogram)
%  is a DISCRETE variable, so we just have to add up the values
bar(xs, counts./sum(counts));

%Simulate the theoretical distribution using the stimated phat and compare
%with the actual data distribution
Y = binopdf(xs,14,phat);
hold on;
plot(xs,Y,'ro-', 'LineWidth', 2, 'MarkerSize', 10);
legend('Simulated', 'Theoretical')
xlabel('Number of quanta')
ylabel('Probability')
title('Exercise 4')


%% Exercise 5
 
p = 0.1:0.1:1;
k = 7;
n= 14;

for i = 1:length(p)
    seven_likelihood(i) = binopdf(k,n,p(i));
end


fprintf('Exercise 5\n')
fprintf ('The most probable release probability to observe release of %d quanta is %3.1f\n', k, p(seven_likelihood == max(seven_likelihood))); 

%Under the null hypothesis (i.e. true release probability = 0.3) the
%probability of 7 releases is:
p =0.3;

p_release_seven = binopdf(k,n,p);

fprintf ('The probability of observing %d quanta under the Null Hypothesis is: %2.3f\n', k, p_release_seven) 

%As the stimate of the true release probability in the temperature
%condition (0.5) is different than the true release probability under the
%Null Hypothesis (0.3) it can be concluded that the temperature have an
%effect on synaptic release
