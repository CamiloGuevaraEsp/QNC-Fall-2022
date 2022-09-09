%% QNC Homework 1 Frequentist Versus Bayesian Approaches

%% Exercise #1

%From a frequentist approach, the p-value is the probability of obtaining
%the data if the hypothesis is true. And from such approach it is no
%possible to make statistical inferences from single data, therefore a 
%p-value or a significancy statement can not be assigned to the result of a single test.  
%% Exercise #2 

%As the false positive rate is
p_false_positive = 0.05;
%And as the test does not have a false negative rate then the probability of being infected 
% after receiving a positive tests is:  :
p_infected_positive = 1- p_false_positive; 

%However, we should know the incidence of the condition (as in exercise 2.1) to make accurate
%predicitions of probability, so, saying that the probability of being
%infected after a positive test is 0.95 is misleading. 

%% Exercise #2.1 Bayes approach 

%Let's call p the priors: 
p = 0:0.1:1; 
%Define the false positive and the true positive rates:
p_positive_infected = 1;
p_false_positive = 0.05;

for i = 1:length(p)
    p_infected = p(i); 
    p_healthy = 1-p_infected; 
    %p(data) is calculated using the law of total probability
    p_data = (p_positive_infected * p_infected) + (p_false_positive * p_healthy);
    %Use Bayes' theorem
    p_infected_data(i) = (p_positive_infected * p_infected)/p_data;    
end


