    function y = geomean(x)
        if any(x<0); y = -inf; end
        y = exp( sum(log(max(x,realmin)))/length(x) );
    end