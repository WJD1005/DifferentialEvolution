function [ret,p] = ranksumtest(X, Y)

	[p, h, stats] = ranksum(X, Y, 'method', 'approximate');
	if h == 0
		ret = '=';
	elseif stats.zval < 0
		ret = '-';
	elseif stats.zval > 0
		ret = '+';
	else
		ret = '?';
	end

end

