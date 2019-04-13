clc
clear all
close all

%Classifying the data into training and testing data
Data_Set= readtable('house_prices_data_training_data.csv');
training_Data = table2array(Data_Set(1:500,4:21));
[row_train,col_train]=size(training_Data );
testing_data=table2array(Data_Set(501:17999,4:21));
[row_test,col_test]=size(testing_data);

%Normalizing both training and testing data 
Normalized_training_Data=var_Normalise(training_Data);
Normalized_testing_data=var_Normalise(testing_data);

[m, n] = size(Normalized_training_Data);
[m_testing, n] = size(Normalized_testing_data);



%Calculating the variance for both training and testing data
for i=1:18
    for j=1:18
    Varfeature_training_data(i)=var(Normalized_training_Data(:,j));
    Varfeature_testing_data(i)=var(Normalized_testing_data(:,j));
    end
end


%Calculating the mean for both training and testing data
mu_training_data = 1/m * sum(Normalized_training_Data);
mu_testing_data = 1/m_testing * sum(Normalized_testing_data);


%calculating the CDF for each example (feature) in the testing and training data set
%using its mean and variance calculated above

for i=1:18

pd_training_data = makedist('Normal','mu',mu_training_data(i),'sigma', Varfeature_training_data(i));
y_training_data(:,i) = cdf(pd_training_data,Normalized_training_Data(:,i));

pd_testing_data = makedist('Normal','mu',mu_testing_data (i),'sigma',Varfeature_testing_data(i));
y_testing_data(:,i) = cdf(pd_testing_data,Normalized_testing_data(:,i));
end


% calculating the  training data epsilons (Probailties)
for i=1:row_train
    
    probabilties_training_data(i,:)=sum(y_training_data(i,:));
   
    
end

% calculating the  testing data epsilons (Probailties)
for i=1:row_test
  
    probabilties_testing_data(i,:)=sum(y_testing_data(i,:));
    
end


%threshold of the epsilon from training data 
min_prob= min( probabilties_training_data);
max_prob=max( probabilties_training_data);


% counting the number of anomalies from the testing data by comparing it
% probabilties to the threshold calculated above

   anomaly_counter=0;
for i=1:row_test
    if probabilties_testing_data(i,:)< min_prob ||  probabilties_testing_data(i,:) > max_prob
        anomaly_counter=anomaly_counter+1;
    end
end
    
        



    