clear;
function dx = kimi(x)
    dx = zeros(2,1);
    I=1;
    m=1;
    g=9.8;
    l=1;
    c=0.8;
    
    
    dx(1)=x(2);
    dx(2)=-(m*g*l*sin(x(1))+c*x(2))/I;
endfunction
