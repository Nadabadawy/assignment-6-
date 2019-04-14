close all
clear
clc

dataset = readtable('house_prices_data_training_data.csv');
Data = table2array(dataset(1:17999, 4:21));
[m n]=size(Data);
data = var_Normalise(Data);

mean_data = mean(data);
std_data = std(data);
pdf_data=zeros(1,18);
for i=1:18
    for j=1:17999
  pdf_data(i,j) = normpdf(data(j,i),mean_data(i),std_data(i));
    end
end 
a=prod(pdf_data);

for i=1:length(a)
if a(i)<0.0000000000001
 anomly(i)=0;
end
if a(i)>0.0000000000001
    anomly(i)=1;
end
end


number_of_ones = sum(anomly);
