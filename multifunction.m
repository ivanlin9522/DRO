function [phat, CIup, CIdown] = multifunction(KaFlag, Discrete, OneSide, SampleSize)

%KaFlag==0-->q_n; KaFlag==1-->X_1,0.95; KaFlag==2-->X_k-1,0.95
k = 5; %k-discretized exponential
n = SampleSize;
alpha = 0.05;
%num is used to count the case when lowerlimit<=Z0(x)<=upperlimit; 
%or Z0(x)<=upperlimit
num = 0;
replication = 1000;
%next we calculate Z0(x);
Largesize =1000000;
KsaiContinuous_true = random('exp',20,[Largesize,1]);
KsaiDiscrete_true = (floor(KsaiContinuous_true./(50/k)+1)*(50/k));
KsaiDiscrete_true(KsaiDiscrete_true>=50)=50;
if Discrete == 1
    Ksai_true = KsaiDiscrete_true; 
elseif Discrete ==0
    Ksai_true = KsaiContinuous_true;
else
    fprint('uh-uh, Discrete can only be 0 or 1');
end


for i=1:replication
    nnum=1;
    %sample n i.i.d. data Ksai from the k-discretized exp(1/20)
    KsaiContinuous = random('exp',20,[n,1]);
    KsaiDiscrete = (floor(KsaiContinuous./(50/k)+1)*(50/k));
    KsaiDiscrete(KsaiDiscrete>=50)=50;
    if Discrete == 1
        Ksai = KsaiDiscrete; 
    elseif Discrete ==0
        Ksai = KsaiContinuous;
    else
        fprint('uh-uh, Discrete can only be 0 or 1');
    end

    %now we calculate the q_n
    if KaFlag==0
        while true
            result = str2double(second_dev_cov_int(10,5,4,3,40,Ksai,n));
            if real(result)>0
                break;
            end
            Ksai = random('exp',20,[n,1]);
        end
        fun = @(u)(1 - chi2cdf(u,1) + result*exp(-u/2)/pi - alpha);    
        q_n = fzero(fun,[0,20]);
        eta = q_n/(2*n);
    
    elseif KaFlag==1
        eta = icdf('Chisquare',1-alpha,1)/(2*n);

    elseif KaFlag==2
        eta = icdf('Chisquare',1-alpha,k-1)/(2*n);
    
    else
        fprint('uh-uh, Kaflag can only be 0 or 1 or 2');
    end
    
    for j=1:20
        %next we calculate upperlimit
        x = 2.5*j;
        TrueMean = mean(h_function(x,Ksai_true));
        
        fun1=@(variable)sum(phiStar((h_function(x,Ksai) + variable(2))./variable(1))*(variable(1)/n))+variable(1)*eta-variable(2);
        variable0=[10000,1];
        A=[-1,0];
        b=0;
        [temp, upperlimit] = fmincon(fun1,variable0, A, b);

        if OneSide==1
            if upperlimit >= TrueMean
                nnum = nnum * 1;
            else
                nnum = 0;
                break;
            end
        %next we calculate lowerlimit
        elseif OneSide==0
            fun2=@(variable)sum(phiStar((h_function(x,Ksai)*(-1)+variable(2))./variable(1))*(variable(1)/n))+variable(1)*eta-variable(2);
            variable0=[10000,1];
            A=[-1,0];
            b=0;
            [temp, Minuslowerlimit] = fmincon(fun2,variable0, A, b);
            lowerlimit=-Minuslowerlimit;
    
            if (upperlimit>TrueMean)&&(lowerlimit<TrueMean)
                nnum = nnum * 1;
            else
                nnum = 0;
                break;
            end
        else
            fprint('uh-uh, OneSide can only be 0 or 1');
        end
    end
    if nnum == 1
        num = num + 1;
    end
end
phat = num / replication;
CIup = phat + icdf('Normal',1-alpha/2,0,1)*sqrt(phat*(1-phat)/replication);
CIdown = phat - icdf('Normal',1-alpha/2,0,1)*sqrt(phat*(1-phat)/replication);
