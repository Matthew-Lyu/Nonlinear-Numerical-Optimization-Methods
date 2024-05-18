function [xk, val, k] = trust_region_method(fun, gfun, Hess, x0)
    n =length(x0); x=x0; dta=1;
    eta1=0.1; eta2=0.75; dtabar=2.0;
    tau1=0.5; tau2=2.0; epsilon=1e-6;
    k=0; Bk = Hess(x);
    while (k<100000)
        gk = gfun(x);
        if(norm(gk)<epsilon)
            break;
        end
        [d,val,lam,ik]=trust_method_quadratic_model(gk,Bk,dta);
        deltaq = -qk(x,d);
        deltaf = fun(x)-fun(x+d);
        rk = deltaf/deltaq;
        if(rk<=eta1)
            dta=tau1*dta;
        else 
            if (rk >= eta2 && norm(d)==dta)
                dta=min(tau2*dta,dtabar);
            else
                dta = dta;
            end
        end
        if(rk > eta1)
            x=x+d;
            Bk=Hess(x);
        end
        k=k+1;
    end
    xk=x;
    val=fun(xk);

    % 信赖域子问题目标函数
    function qd= qk(x,d)
        gkk = gfun(x);  Bkk=Hess(x);
        qd=gkk'*d + 0.5*d'*Bkk*d;
    end
end

  