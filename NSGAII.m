% The implemention for NSGA-II 
clear;
rng(100)
f = @(x)(MOP(x));

% initilization
D = 10;
Xmin = ones(1, D)  * 0;
Xmax = ones(1, D)  * 1;
lu = [Xmin; Xmax];

MaxFE = 25000;
nPop = 100;

% hyperparameter for GA
pc=0.7;                                 % crossover percentage
nc=2*round(pc*nPop/2);      % number of offsprings (also Parnets)
pm=0.3;                                % mutation percentage
nm=nPop - nc;           % number of mutants

gamma=0.4;                          % extra range factor for crossover
mu=0.1;  

% population
Pop = repmat(Xmin, nPop, 1) + repmat((Xmax - Xmin), nPop, 1) .* rand(nPop, D);
Fit = f(Pop);
fe = nPop;

while fe < MaxFE
    
    ssr = randperm(nPop);
    parent = Pop(ssr, :);
    parentfit = Fit(ssr, :);
    
    % select parents for crossover
    spc = ssr(1: nc);
    % select parents for mutation
    spm = ssr(nc+1: end);
    
    % Crossover
    parentc = Pop(spc, :);
    offspringc = zeros(size(parentc));
    for k=1:nc/2
        
        p1=parentc(2*k-1, :);
        p2=parentc(2*k, :);
        
        [offspringc(2*k-1, :), offspringc(2*k, :)]=Crossover(p1, p2, gamma, Xmax, Xmin);
        
    end
    
    
    % Mutation
    parentm = Pop(spm, :);
    offspringm = zeros(size(parentm));
    for k=1:nm 
        
        p=parentm(k, :);
        
        offspringm(k, :)=Mutate(p,mu, Xmax, Xmin);
        
    end

    % new offspring
    offspring = [offspringc; offspringm];
    offspringfit = f(offspring);
    fe = fe + size(offspring, 1);
    
    % Merge the parent and offspring
    Pop = [parent; offspring];
    Fit = [parentfit; offspringfit];
    
    
    % non-domniated-sort
     [F, F1] = NonDominatedSorting(Pop, Fit);

    % crowding distance assignment
    CCDI = CalcCrowdingDistance(Fit, F);

    % sort population
    indexDC = SortPopulation(F1, CCDI);
    
    % truncation 
    Pop = Pop(indexDC(1:nPop), :);
    Fit = Fit(indexDC(1:nPop), :);
    
    % Plot results;
    plot(Fit(:, 1), Fit(:, 2), 'r*', 'MarkerSize', 8);
    xlabel('1^{st} Objective');
    ylabel('2^{nd} Objective');
    title('Iterated Solutions');
    axis([0 1 0 5])
    pause(0.00001);
end


% non-dominated sorting
function [F, F1] = NonDominatedSorting(pop, fit)
    [n, D] = size(pop);
    S = cell(n, 1);
    F = {[]};
    nDom = zeros(n, 1);
    rrank = zeros(n, 1);
    for p=1:n
        S{p} = [];
        for q = 1:n
            if all(fit(p, :) <= fit(q, :))
                S{p} = [S{p} q];
            elseif all(fit(q, :) <= fit(p, :))
                nDom(p) = nDom(p) + 1;
            end
        end
        if nDom(p) == 0
            rrank(p) = 1;
            F{1} = [F{1} p];
        end
    end
    
    i =1;
    while true
        Q = [];
        for p=F{i}
            for q = S{p}
                nDom(q) = nDom(q) - 1;
                if nDom(q) == 0
                    rrank(q) = i + 1;
                    Q = [Q q];
                end
            end
        end
        
        if isempty(Q)
            break;
        end
        i = i + 1;
        F{i} = Q;

    end
    
    F1 = zeros(n, 1);
    nfront = numel(F);
    for j=1:nfront
        for k = F{j}
            F1(k) = j;
        end
    end
    
end

function CCDI = CalcCrowdingDistance(fit, F)
    
    nfront = numel(F);
    ojn = size(fit, 2);
    CCDI = zeros(size(fit, 1), 1);
    
    for j=1:nfront
        nsol = numel(F{j});
        CDI = zeros(1, nsol);
        tmp = fit(F{j}, :);
        for oj = 1:ojn
            [value, index] = sort(tmp(:, oj));
            fm_max = value(end);
            fm_min = value(1);
            CDI(index(1)) = Inf;
            CDI(index(end)) = Inf;
            for i=2:nsol-1
                dd = ( fit(F{j}(index(i+1)), oj) - fit(F{j}(index(i-1)), oj)  ) / (fm_max - fm_min);
                CDI(index(i)) = CDI(index(i)) + dd;
            end
        end
        
        for k=1:nsol
            CCDI(F{j}(k)) = CDI(k);
        end
        
    end

end

function index = SortPopulation(F, CDI)
    [~,  index] = sortrows([F CDI], [1, -2]);
end

% crossover operator
function [y1, y2]=Crossover(x1,x2,gamma,VarMax,VarMin)

    alpha=unifrnd(-gamma,1+gamma, size(x1));
    
    y1=alpha.*x1+(1-alpha).*x2;
    y2=alpha.*x2+(1-alpha).*x1;
    
    y1=max(y1,VarMin);
    y1=min(y1,VarMax);
    
    y2=max(y2,VarMin);
    y2=min(y2,VarMax);

end

% mutation operator
function y=Mutate(x,mu,VarMax,VarMin)
       
    % uniform mutation 
    r = rand(size(x)) >= mu;
    y = unifrnd(VarMin, VarMax, size(x));
    y(r) = x(r);

end


