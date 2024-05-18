%%% 利用光滑牛顿法求解信赖域子问题 %%%
function  [d, val, lam, k] = trust_region_quadratic_model(gk, Bk, dta)
    % 光滑牛顿法
    % 输入：
    %   gk: 当前点的梯度
    %   Bk: 当前点的 Hessian 矩阵
    %   delta_k: 信赖域半径
    % 输出：
    %   d: 子问题的解
    %   val: 目标函数在 d 处的值
    %   lam: 最优拉格朗日乘子
    %   k: 迭代次数

    n = length(gk); 
    gamma=0.05;
    epsilon=1.0e-6;  
    rho=0.6;  
    sigma=0.2;
    mu0=0.05;  
    lam0=0.05;
    d0=ones(n,1);  
    u0=[mu0, zeros(1,n+1)]';
    z0=[mu0, lam0,d0']';
    k=0; %k为迭代次数
    z = z0; mu=mu0; lam=lam0; d=d0;
    while (k <= 150)
        dh=dah(mu,lam,d,gk,Bk,dta);
        if(norm(dh)<epsilon)
            break;
        end
        A=JacobiH(mu,lam,d,Bk,dta);
        b=beta(mu,lam,d,gk,Bk,dta,gamma)*u0-dh;
        B=inv(A);   dz=B*b;
        dmu=dz(1); dlam=dz(2); dd=dz(3:n+2);
        m=0;  mk=0;
        while (m<20)
            dhnew=dah(mu+rho^m*dmu,lam+rho^m*dlam,d+rho^m*dd,gk,Bk,dta);
            if(norm(dhnew)<=(1-sigma*(1-gamma*mu0)*rho^m)*dh)
                mk=m;
                break;
            end
            m=m+1;
        end
        alpha=rho^mk;
        mu=mu+alpha*dmu;
        lam=lam+alpha*dlam;
        d=d+alpha*dd;
        k=k+1;
    end
    val=gk'*d+0.5*d'*Bk*d;
    
    function p=phi(mu,a,b)
        p=a+b-sqrt((a-b)^2+4*mu);
    end

    function dh=dah(mu,lam,d,gk,Bk,dta)
        n1 = length(d);
        dh(1)=mu;  dh(2)=phi(mu,lam, dta^2-norm(d)^2);
        mh=(Bk+lam*eye(n1))*d+gk;
        for (i=1:n1)
            dh(2+i)=mh(i);
        end
        dh=dh(:);
    end

    function bet=beta(mu,lam,d,gk,Bk,dta,gamma)
        dhh=dah(mu,lam,d,gk,Bk,dta);
        bet=gamma*norm(dhh)*min(1,norm(dhh));
    end

    % 定义雅可比矩阵
    function A=JacobiH(mu,lam,d,Bk,dta)
        n2=length(d);
        A=zeros(n2+2,n2+2);
        pmu=-4*mu/sqrt((lam+norm(d)^2-dta^2)^2+4*mu^2);
        thetak=(lam+norm(d)^2-dta^2)/sqrt((lam+norm(d)^2-dta^2)^2+4*mu^2);
        A = [1,             0,            zeros(1,n2);
            pmu,          1-thetak,     -2*(1+thetak)*d';
            zeros(n2,1),    d,        Bk+lam*eye(n2)];
    end
end