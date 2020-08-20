function [xq,approx,beta,modelfun] = xregFitDistortion(sbx, x, y, fitType)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
if nargin<4
    fitType = 'cos2';
end

xq = 1:size(sbx,2);

if strcmp(fitType,'poly3')
    p = polyfit(x,y,3);
    approx = polyval(p,xq);
elseif strcmp(fitType,'cos2')
    modelfun = @(b,x)(b(1)-b(1)*(cos(1/b(2) * (x - b(3))))).^2.*sign(x-b(3));
    beta0 = [10, 600, 300];
    beta = nlinfit(x,y,modelfun,beta0);
    approx=modelfun(beta,xq);
elseif strcmp(fitType,'cos')
    modelfun = @(b,x)(b(1)-b(1)*cos(1/b(2) * (x - b(3)))).*sign(x-b(3));
    beta0 = [10, 600, 300];
    beta = nlinfit(x,y,modelfun,beta0);
    approx=modelfun(beta,xq);
elseif strcmp(fitType,'coscos')
    modelfun = @(b,x)(b(1)-b(1)*cos(1/b(2) * (x - b(3))).*cos(1/b(4) * (x - b(3))).*sign(x-b(3)));
    beta0 = [10, 600, 300,600];
    beta = nlinfit(x,y,modelfun,beta0);
    approx=modelfun(beta,xq);
end

end

