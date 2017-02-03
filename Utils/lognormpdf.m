function p = lognormpdf(x,mu,sigma)
% Dealing with log probability, use this when the normal pdf comes up

p = -min(log(sqrt(2*pi)*sigma)+(x - mu).^2/(2*sigma.^2),realmax);

end
