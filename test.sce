clear
function ans=b(x)
    ans=x^2;
endfunction



function ans=a(f,x)
    ans=f(2)*2;
endfunction

disp(a(b,2));
