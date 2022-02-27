clear ;

function [f, g, ind]=cost(x, ind)
    xref = [1; 2; 3];
    f = 0.5 * norm(x - xref)^2;
    g = x - xref;
endfunction

// Simplest call
x0 = [1; -1; 1];
[fopt, xopt] = optim(cost, x0)
