function [sqrtP_predict] = CalculateErrorCovarianceSqrt(P_predict)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
sqrtP_predict = chol(P_predict,'lower');
end

