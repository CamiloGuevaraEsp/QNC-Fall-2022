clear all 
close all
%Simulate the samples
sample_size = 50;
mu = 1;
sigma = 2;
num_simulations = 1000; 
alpha = 0.05;
% Apply bonferroni correction
bonferroni_p_value = alpha/num_simulations;
%Run the simulations
idx = [];
while length(idx) == 0;
    for i = 1:num_simulations
    X1 = normrnd(mu, sigma, sample_size, 1);
    X2 = normrnd(mu, sigma, sample_size, 1);
    [~,p] = ttest2(X1,X2);
    p_values(i) = p;
    end
    %Apply Benjamini-Hochber correction
    %Sort p-values in ascending order and assume a FDR of 0.05
    sorted_p = sort(p_values); 
    FDR = 0.05;
    %Calculate critical value for each p-value
    for i = 1: length(p_values)
        sorted_p(2,i) = (i/num_simulations) * FDR; 
    end
    idx = find(sorted_p(sorted_p(1,:) < sorted_p(2,:))); 
end
fprintf('Bonferroni corrected p-value = %2.5f\n',bonferroni_p_value);
fprintf('Percentage of significant results with Bonferroni correction : %2.3f\n', (sum(sorted_p(1,:) < bonferroni_p_value)/num_simulations) *100);
fprintf('Benjamini-Hochber corrected p-value = %2.3f \n',sorted_p(1,idx(end)))
fprintf('Percentage of significant results with Benjamini-Hochber correction : %2.8f\n', ((length(idx))/num_simulations) *100);
%%
clear all
%Simulate the samples
sample_size = 50;
mu1 = 1;
mu2 = 2;
sigma = 2;
num_simulations = 1000; 
alpha = 0.05;
% Apply bonferroni correction
bonferroni_p_value = alpha/num_simulations;
%Run the simulations
idx = [];
while length(idx) == 0;
    for i = 1:num_simulations
    X1 = normrnd(mu1, sigma, sample_size, 1);
    X2 = normrnd(mu2, sigma, sample_size, 1);
    [~,p] = ttest2(X1,X2);
    p_values(i) = p;
    end
    %Apply Benjamini-Hochber correction
    %Sort p-values in ascending order and assume a FDR of 0.05
    sorted_p = sort(p_values); 
    FDR = 0.05;
    %Calculate critical value for each p-value
    for i = 1: length(p_values)
        sorted_p(2,i) = (i/num_simulations) * FDR; 
    end
    idx = find(sorted_p(sorted_p(1,:) < sorted_p(2,:))); 
end  
fprintf('Different means Bonferroni corrected p-value = %2.5f\n',bonferroni_p_value);
fprintf('Percentage of significant results with Bonferroni correction : %2.3f\n', (sum(sorted_p(1,:) < bonferroni_p_value)/num_simulations) *100);
fprintf('Different means Benjamini-Hochber corrected p-value = %2.3f \n',sorted_p(1,idx(end)));
fprintf('Percentage of significant results with Benjamini-Hochber correction : %2.8f\n', ((length(idx))/num_simulations) *100);
%% Now run 10 simulations incresing the difference in means
clear corrected_p
%Define some constants
sigma = 2;
num_simulations = 1000; 
alpha = 0.05;
delta = 0:0.1:2.6;
sample_size = 50;
bonferroni_p_value = alpha/num_simulations;
%Simulate the samples
for i = 1:length(delta)
    mu1 = 1;
    mu2 = mu1 + delta(i);
    clear sorted_p idx p_values
    idx = [];
    while length(idx) == 0;
        for k = 1:num_simulations
        X1 = normrnd(mu1, sigma, sample_size, 1);
        X2 = normrnd(mu2, sigma, sample_size, 1);
        [~,p] = ttest2(X1,X2);
        p_values(k) = p;
        end
        %Apply Benjamini-Hochber correction
        %Sort p-values in ascending order and assume a FDR of 0.05
        sorted_p = sort(p_values); 
        FDR = 0.05;
        %Calculate critical value for each p-value
        for j = 1: length(p_values)
            sorted_p(2,j) = (j/num_simulations) * FDR; 
        end
        idx = find(sorted_p(sorted_p(1,:) < sorted_p(2,:)));
    end
    corrected_p(i) = sorted_p(1,idx(end));
    percentage_significant(i) = (length(idx)/num_simulations) *100;
    percentage_bonferroni(i) = (sum(sorted_p(1,:) < bonferroni_p_value)/num_simulations) *100;
end
figure
plot(delta,percentage_significant,'-o','LineWidth',1.5);
hold on
plot(delta,percentage_bonferroni,'-o','LineWidth',1.5);
xlabel('delta means');
ylabel('% Significant values');
legend('Benjamini-Hochber', 'Bonferroni');
%As the difference in mean increase, the percentage of significant values
%increases, this increase is higher in the Benjamini-Hochber correction
%compared with the Bonferroni correction 



