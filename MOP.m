function y = MOP(x)

n = size(x, 1);
y = zeros(n, 2);
for i=1:n
    y(i, :) = ZDT1(x(i, :));
end

end


function y = ZDT1(x)
    y1 =  x(1);
    g = 1 + 9 * sum(x(2:end)) / ( numel(x)-1);
    y2 = g * (1 - (x(1)/g)^2);
    y =[y1 y2];
end