function [y] = dixmaana(x)
	y1 = 0;	y2 = 0; y3 = 0; y4 = 0;	y = 0;
	n = length(x);
	m = n/3;
	a = 1; b = 0.0625; c = 0.0625; d = 0.0625;
	k1 = 2; k4 = 2;

	for i = 1:n
		y1 = y1 + a*x(i)^2*(i/n)^k1;
	end

	for i = 1:n-1
		y2 = y2 + b*x(i)^2*(x(i+1)+x(i+1)^2)^2;
	end

	for i = 1:2*m
		y3 = y3 + c*x(i)^2*x(i+m)^4;
	end

	for i = 1:m
		y4 = y4 + d*x(i)*x(i+2*m)*(i/n)^k4;
	end
	
	y = 1 + y1 + y2 + y3 + y4;
end
	
	
