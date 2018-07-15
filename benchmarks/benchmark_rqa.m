% Measure the times (in ms) of evaluating an expression n times
function [t, res] = measuretime(f, n)
    t = zeros(1,n);
    for i=1:n
        t0 = tic;
        res = f();
        t(n) = 1000*toc;
    end
end

% Function that will be measured
function y = fun_rqa(x,metric)
    delay = 6;
    embed = 3;
    radius = 1.2;
    y = crqa(x,x,embed,delay,radius,[],[],2,2,metric,'nonormalize','nogui');
end

% Analyse 12 series from 250 to 3000 points 
% (With variable metric)
function benchmark(metric)
    m = dlmread('rossler.txt');
    for r=1:12
        x=m(1:250*r,2*r-1);
        [tt, res] = measuretime(fun_rqa(x,metric),5);
        t = median(tt);
        % Write table of results
        f = fopen(sprintf('benchmark_rqa_julia_%s.txt',metric),'a');
        fprintf(f, '%d\t%f\t', r,t);
        % "RR","DET","L","Lmax","ENT","LAM","TT"
        for k=1:7
            fprintf(f, '%f\t', res(1,k));
        end
        fprintf(f, '\n');
        end
    end
end
% Do it with max and euclidean norms
benchmark('euclidean');
benchmark('maxnorm');
