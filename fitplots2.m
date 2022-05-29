function [positions nrows ncols ]=fitplots2(N,orientation,spacing)
% djm 10/11/10
if length(N)==2
    nrows=N(1);
    ncols=N(2);
else
    if ~exist('orientation','var') || isempty(orientation) || strcmp(orientation,'portrait')
        nrows=ceil(sqrt(N));
        ncols=floor(sqrt(N));
        if nrows*ncols<N,
            if nrows*(ncols+1)<N, nrows=nrows+1;
            else ncols=ncols+1;
            end
        end
    else
        nrows=floor(sqrt(N));
        ncols=ceil(sqrt(N));
        if nrows*ncols<N,
            if (nrows+1)*ncols<N, ncols=ncols+1;
            else nrows=nrows+1;
            end
        end
    end
end

if ~exist('spacing','var'), spacing='fixed';end
if isnumeric(spacing), fatness=spacing;
else
    switch spacing
        case 'fixed', fatness=6;
        case 'sparse', fatness=4;
            
    end
end

% gap is 1/6 width of each plot
xgap=1/(ncols*(fatness+1)+1);
ygap=1/(nrows*(fatness+1)+1);
width=xgap*fatness;
height=ygap*fatness;

% work out all possible positions
positions={}; %zeros(N,4);
p=0;
for y=1:nrows
    for x=1:ncols
        thispos=[...
            xgap+(x-1)*(xgap+width), ...
            1-y*(ygap+height), ...
            width, height];
            positions{y,x}=thispos;
    end
end

return