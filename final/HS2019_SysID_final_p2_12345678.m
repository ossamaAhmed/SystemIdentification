function [p2_a_ex1, p2_b_ex1,p2_a_ex2, p2_b_ex2, p2_mse_ex1, p2_mse_ex2] = HS2019_SysID_final_p2_12345678()
%% Solution for Problem 1

%% Generate data

% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

% generate data
[p2_u1,p2_u2,p2_u_cv,p2_y1,p2_y2,p2_y_cv] = HS2019_SysID_final_p2_GenerateData(LegiNumber);

%% General instructions for solution

% Change the filename of this function, both in the function definition
% above and in the filename in the folder

% Use the variables p2_y and p2_u to solve the problem. 

% Modify your code in the next sections, and return the variables
% requested.

% If you skip one part of the problem, return the empty vectors as already
% provided in the code

%% Part 1
% Compute parameter estimates using instruments generated by least squares
% regression of the data

p2_a_ex1 = [NaN,NaN,NaN]';
p2_b_ex1 = [NaN,NaN]';
%% Part 2
% Compute parameter estimates using delayed inputs as instruments

p2_a_ex2 = [NaN,NaN,NaN]'; 
p2_b_ex2 = [NaN,NaN]'; 
%% Part 3
% Compute mean squared errors using the cross-validation data set

p2_mse_ex1 = NaN; 
p2_mse_ex2 = NaN;
end

