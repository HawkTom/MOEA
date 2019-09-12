% The implemention of SPEA2
clear;
rng(100);
f = @(x)(MOP(x));

% initilization
D = 10;
Xmax = ones(1, nPop) * 1;
Xmin = ones(1, nPop) * 0;
lu = [Xmin;Xmax];

MaxFE = 10000;
nPop = 100;

% population initilization
Pop = repmat((Xmax - Xmin), nPop, 1) .* rand(nPop, D) + repmat(Xmin, nPop, 1);
Fit = f(Pop);
fe = nPop;












while fe < maxFE
    
    % Fitness Assignment
    
    % Environment Selection 
    
    % 
    
    break;
    
end
