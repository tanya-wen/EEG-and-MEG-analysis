y = mean(accuracy_vals)-50;
y = round(y+abs(min(y)));
% y = mean(accuracy_vals);
% y = round(y);

for i = 1:numel(y)
A{i}=ones(1,y(i))*i;
end

y=[A{:}];

[ests,fVal,exitFlag,solverOutput] = exgauss_fit(y);
nBin        = ceil(sqrt(numel(y)));
[N,binCtr]  = hist(y,nBin);
hist(y,nBin); hold on
binWidth    = range(binCtr)/nBin;
f           = exgauss_pdf(y,ests);
fNorm       = numel(y)*f*binWidth;
plot(y, fNorm, 'Color','r');