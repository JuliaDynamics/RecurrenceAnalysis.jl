function benchmark_rqa()
% Do it with max and euclidean norms
benchmark('euclidean');
benchmark('maxnorm');
end

% Measure the times (in ms) of evaluating an expression n times
function [t, res] = measuretime(f, n, varargin)
    t = zeros(1,n);
    for i=1:n
        t0 = tic;
        res = f(varargin{:});
        t(i) = 1000*toc(t0);
    end
end

% Function that will be measured
function y = fun_rqa(x,metric)
    delay = 6;
    embed = 3;
    radius = 1.2;
    y = crqa(x,x,embed,delay,radius,[],[],2,2,1,metric,'nonormalize','nogui');
end

% Analyse 12 series from 250 to 3000 points 
% (With variable metric)
function benchmark(metric)
    m = dlmread('rossler.txt');
    for r=1:12
        x=m(1:250*r,2*r-1);
        [tt, res] = measuretime(@fun_rqa,5,x,metric);
        t = median(tt);
        % Write table of results
        f = fopen(sprintf('benchmark_rqa_matlab_%s.txt',metric),'a');
        fprintf(f, '%d\t%f\t', r,t);
        % "RR","DET","L","Lmax","ENT","LAM","TT"
        for k=1:7
            fprintf(f, '%f\t', res(1,k));
        end
        fprintf(f, '\n');
        fclose(f);
    end
end
