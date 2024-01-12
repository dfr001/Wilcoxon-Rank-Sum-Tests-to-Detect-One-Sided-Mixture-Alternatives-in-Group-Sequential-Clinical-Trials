import numpy as np
from scipy.stats import norm
from scipy.stats import logistic
from scipy.stats import laplace
from scipy.stats import t
from scipy.stats import multivariate_normal
from scipy.stats import rankdata
from scipy.optimize import minimize_scalar
import scipy.integrate as integrate
import math

def grp_seq_mix(num_stages=1,alpha=0.05,beta=0.2,rho=None,theta=1,K=1,mu=0,sigma=1,ratio=1,f="normal",method="sar",norm_approx=True,B=10000):
    k = num_stages
    m_start = 1
    std=True
    if method=="sr": method="rerank"
    pi1 = np.append([alpha*1/k**rho],np.diff(alpha*(np.arange(1,1+k)/k)**rho))
    pi2 = np.append([beta*1/k**rho],np.diff(beta*(np.arange(1,1+k)/k)**rho))
    d = f
    arg1=mu; arg2=sigma
    if f=="normal":
        def rpsi(n, arg1=0, arg2=1): return np.random.normal(loc=arg1, scale=arg2, size=n)
        def ppsi(x, arg1=0, arg2=1): return norm.cdf(x, loc=arg1, scale=arg2)
        def dpsi(x, arg1=0, arg2=1): return norm.pdf(x, loc=arg1, scale=arg2)
        def qpsi(q, arg1=0, arg2=1): return norm.ppf(q, loc=arg1, scale=arg2)
    elif f=="logistic":
        def rpsi(x, arg1=0, arg2=1): return np.random.logistic(loc=arg1, scale=arg2*math.sqrt(3)/math.pi, size=n)
        def ppsi(x, arg1=0, arg2=1): return logistic.cdf(x, loc=arg1, scale=arg2*math.sqrt(3)/math.pi)
        def dpsi(x, arg1=0, arg2=1): return logistic.pdf(x, loc=arg1, scale=arg2*math.sqrt(3)/math.pi)
        def qpsi(q, arg1=0, arg2=1): return logistic.ppf(q, loc=arg1, scale=arg2*math.sqrt(3)/math.pi)
    elif f=="laplace":
        def rpsi(n, arg1=0, arg2=1): return np.random.laplace(loc=arg1, scale=arg2/math.sqrt(2), size=n)
        def ppsi(x, arg1=0, arg2=1): return laplace.cdf(x, loc=arg1, scale=arg2/math.sqrt(2))
        def dpsi(x, arg1=0, arg2=1): return laplace.pdf(x, loc=arg1, scale=arg2/math.sqrt(2))
        def qpsi(q, arg1=0, arg2=1): return laplace.ppf(q, loc=arg1, scale=arg2/math.sqrt(2))
    elif f=="t":
        def rpsi(n, arg1=0, arg2=1): return np.random.standard_t(df=3, size=n)*arg2/math.sqrt(3) + arg1
        def ppsi(x, arg1=0, arg2=1): return t.cdf(x, df=3, loc=arg1, scale=arg2/math.sqrt(3))
        def dpsi(x, arg1=0, arg2=1): return t.pdf(x, df=3, loc=arg1, scale=arg2/math.sqrt(3))
        def qpsi(q, arg1=0, arg2=1): return t.ppf(q, df=3, loc=arg1, scale=arg2/math.sqrt(3))
    pxy = integrate.quad(lambda y: ppsi(y)*((1-theta)*dpsi(y)+theta*dpsi(y-K)),float("-inf"),float("inf"))[0]
    pxyy = integrate.quad(lambda y: (1-(1-theta)*ppsi(y)-theta*ppsi(y-K))**2*dpsi(y),float("-inf"),float("inf"))[0]
    pxxy = integrate.quad(lambda y: (ppsi(y))**2*((1-theta)*dpsi(y)+theta*dpsi(y-K)),float("-inf"),float("inf"))[0]
    # jeske/yao formula for first estimate
    def h1(y, theta, K):
        return ppsi(y)*((1-theta)*dpsi(y) + theta*dpsi(y-K))
    gamval = integrate.quad(h1,float("-inf"),float("inf"),args=(theta,K))[0]
    def h2(y, theta, K):
        return ((1-(1-theta)*ppsi(y)-theta*ppsi(y-K))**2)*dpsi(y)
    xi1 = integrate.quad(h2,float("-inf"),float("inf"),args=(theta,K))[0] - gamval**2
    def h3(y, theta, K):
        return (ppsi(y))**2*((1-theta)*dpsi(y) + theta*dpsi(y-K))
    xi2 = integrate.quad(h3,float("-inf"),float("inf"),args=(theta,K))[0] - gamval**2
    # lambda = 1/2
    m_approx = ((norm.ppf(1-alpha)/math.sqrt(6) + norm.ppf(1-beta)*math.sqrt(xi1+xi2))/(gamval-1/2))**2
    m_approx2 = (((norm.ppf(1-alpha)*math.sqrt((ratio+1)/(12*ratio))) + norm.ppf(1-beta)*math.sqrt(xi1+xi2/ratio))/(gamval-1/2))**2
    
    delta = K*sigma
    if norm_approx:
        if num_stages==1:
            m = 0
            p = 1
            while p > pi2:
                m += 1
                n = math.ceil(ratio*m)
                mu = n*(m+n+1)/2
                v = m*n*(m+n+1)/12
                r = norm.ppf(1-pi1)
                p = norm.cdf(r,loc=(n*(m*pxy + (n+1)/2) - mu)/math.sqrt(v),
                               scale=math.sqrt((m*n*pxy*(1-pxy) + m*n*(n-1)*(pxyy - pxy**2) + m*(m-1)*n*(pxxy-pxy**2)))/math.sqrt(v))
            d = dict();
            d['m'] = m
            d['n'] = n
            d['r'] = np.around(r,4)
            d['Type I'] = np.around(1-norm.cdf(r),4)
            d['Type II'] = np.around(p,4)
            return d
        else:
            print('2+ stages')
            delta = K*sigma
            m = 0
            p = 1
            while p >= pi2[-1]:
                m = m+1
                n = math.ceil(ratio*m)
                r = np.zeros(k)
                a = np.zeros(k-1)

                # establish distributions
                if method=="sar" and std:
                    k_range = np.arange(1,1+k)
                    mu = n*(m+n+1)/2
                    v = m*n*(m+n+1)/12
                    null_mu = np.zeros(k)
                    null_var = v/k_range
                    null_cov_num = k_range*k_range[:,None]
                    null_cov_den = np.zeros((k,k))
                    for i1 in k_range:
                        null_cov_den[:,i1-1] = i1
                        null_cov_den[i1-1,:] = i1
                    null_cov = null_cov_num**0.5/null_cov_den
                    alt_mu = (m*n*pxy + n*(n+1)/2 - mu)/(v**(1/2)/(k_range**(1/2)))
                    alt_v = (m*n*pxy*(1-pxy) + m*n*(n-1)*(pxyy - pxy**2) + m*(m-1)*n*(pxxy-pxy**2))/v
                    alt_cov = null_cov*alt_v
                elif method=="rerank" and std:
                    k_range = np.arange(1,1+k)
                    mu = k_range*n*(k_range*(m+n)+1)/2
                    v = k_range**2*m*n*(k_range*(m+n)+1)/12
                    null_mu = np.zeros(k)
                    null_var = np.ones(k)
                    null_cov = np.zeros((k,k))
                    for i1 in k_range:
                        for i2 in k_range:
                            s = min(i1,i2)
                            sp = max(i1,i2)
                            null_cov[i1-1,i2-1] = ((s**2*m*n)*(1/2 - (sp*(m+n)-1)*0.5**2 + (sp*n-1)*1/3 + (sp*m-1)*1/3))/math.sqrt(v[i1-1]*v[i2-1])
                    alt_mu = (k_range*n*(k_range*m*pxy + (k_range*n+1)/2) - mu)/v**(1/2)
                    alt_var = ((k_range*m)*(k_range*n)*(pxy - pxy**2*((k_range*m)+(k_range*n)-1) + ((k_range*n)-1)*pxyy + ((k_range*m)-1)*pxxy))/v
                    alt_cov = np.zeros((k,k))
                    for i1 in k_range:
                        for i2 in k_range:
                            s = min(i1,i2)
                            sp = max(i1,i2)
                            alt_cov[i1-1,i2-1] = ((s**2*m*n)*(pxy - (sp*(m+n)-1)*pxy**2 + (sp*n-1)*pxyy + (sp*m-1)*pxxy))/math.sqrt(v[i1-1]*v[i2-1])
                
                # critical values
                r[0] = norm.ppf(1-pi1[0], loc=null_mu[0], scale=math.sqrt(null_cov[0,0]))
                a[0] = norm.ppf(pi2[0], loc=alt_mu[0], scale=math.sqrt(alt_cov[0,0]))
                for i1 in np.arange(1,k):
                    # finding upper critical values r
                    # lambda x: multivariate_normal.cdf(np.append([r[:i1-1],float("inf")]),mean=np.zeros(k)[:i1],cov=temp[:i1,:i1],lower_limit=np.append([a[:i1-1],x]))
                    r[i1] = minimize_scalar(lambda x: (multivariate_normal.cdf(np.append(r[:i1],float("inf")), mean=null_mu[:i1+1], cov=null_cov[:i1+1,:i1+1], lower_limit=np.append(a[:i1],x)) - pi1[i1])**2, method="brent")["x"]
                    # finding lower critical values a, need k-1 of these since there is only one critical value for final stage
                    if i1!=np.arange(1,k)[-1]:
                        a[i1] = minimize_scalar(lambda x: (multivariate_normal.cdf(np.append(r[:i1],x), mean=alt_mu[:i1+1], cov=alt_cov[:i1+1,:i1+1], lower_limit=np.append(a[:i1],float("-inf"))) - pi2[i1])**2, method="brent")["x"]
                p = multivariate_normal.cdf(r, mean=alt_mu, cov=alt_cov, lower_limit=np.append(a,float("-inf")))

            # get probabilities associated with critical values
            p_r = 1-norm.cdf(r[0],loc=null_mu[0], scale=math.sqrt(null_cov[0,0]))
            p_a = norm.cdf(a[0], loc=alt_mu[0], scale=math.sqrt(alt_cov[0,0]))
            for i1 in np.arange(1,k):
                if i1!=np.arange(1,k)[-1]:
                    p_r = np.append(p_r, multivariate_normal.cdf(np.append(r[:i1],float("inf")), mean=null_mu[:i1+1], cov=null_cov[:i1+1,:i1+1], lower_limit=np.append(a[:i1],r[i1])))
                    p_a = np.append(p_a, multivariate_normal.cdf(np.append(r[:i1],a[i1]), mean=alt_mu[:i1+1], cov=alt_cov[:i1+1,:i1+1], lower_limit=np.append(a[:i1],float("-inf"))))
                else:
                    p_r = np.append(p_r, multivariate_normal.cdf(np.append(r[:i1],float("inf")), mean=null_mu[:i1+1], cov=null_cov[:i1+1,:i1+1], lower_limit=np.append(a[:i1],r[i1])))
                    p_a = np.append(p_a, p)
            p_r = np.append(p_r, sum(p_r))
            p_a = np.append(p_a, sum(p_a))
            d = dict();
            d['m'] = m
            d['n'] = n
            d['r'] = np.around(r,4)
            d['a'] = np.around(a,4)
            d['Type I'] = np.around(p_r,4)
            d['Type II'] = np.around(p_a,4)
            return d
    else: # simulate distribution instead of normal approximation
        if k==1: # just one stage
            delta = K*sigma
            if m_approx2>1: m_start = math.floor(m_approx2-10)
            def inc_func(incr=1, m=m_start):
                p = 1
                while(p>=pi2):
                    if m<15: m+=1
                    else: m+=incr
                    n = math.ceil(ratio*m)
                    w = np.zeros(B)
                    w_alt = np.zeros(B)
                    for i in np.arange(B):
                        x_null = rpsi(n=k*m, arg1=arg1, arg2=arg2)
                        y_null = rpsi(n=k*n, arg1=arg1, arg2=arg2)
                        x_alt = rpsi(n=k*m, arg1=arg1, arg2=arg2)
                        y_alt = rpsi(n=k*n, arg1= arg1 + delta*np.random.default_rng().choice(2,size=k*n,p=[1-theta,theta]), arg2=arg2)
                        w[i] = sum(rankdata([x_null, y_null])[m:])
                        w_alt[i] = sum(rankdata([x_alt, y_alt])[m:])

                    p_r = 0
                    index_max = m*n + n*(n+1)/2
                    u = sorted(np.unique(w), reverse=True)
                    index = index_max+1
                    while p_r <= pi1:
                        p_rej = p_r
                        index -= 1
                        p_r = sum(w>=index)/B
                        if index==1: break
                    if index==index_max and p_r>pi1:
                        r = float("inf")
                    elif p_r<=pi1:
                        r = index
                    else:
                        r = index+1
                    if std==True:
                        mu = n*((m+n)+1)/2; sig = math.sqrt(m*n*((m+n)+1)/12)
                        r = (r-mu)/sig
                        w = (w-mu)/sig; w_alt = (w_alt - mu)/sig
                    p = sum(w_alt<r)/B
                p_rej = sum(w>=r)/B
                p_acc = sum(w_alt<r)/B
                return [m,n,r,p_rej,p_acc]

            incr_size = 10
            res = inc_func(incr=incr_size)
            print(res)
            print(res[0])
            if res[0]>=15:
                res = inc_func(incr = 1, m = res[0]-incr_size-1)
                d = dict();
                d['m'] = res[0]
                d['n'] = res[1]
                d['r'] = res[2]
                d['Type I'] = res[3]
                d['Type II'] = res[4]
                return d
            else:
                d = dict();
                d['m'] = res[0]
                d['n'] = res[1]
                d['r'] = res[2]
                d['Type I'] = res[3]
                d['Type II'] = res[4]
                return d
        else: # 2+ stages
            k_range = np.arange(1,1+k)
            delta = K*sigma
            p_start = 1
            # use normal approx for one stage to get close initially
            if m_approx2/k-10>m_start: m_start = math.floor(m_approx2/k - 10)
            def incr_arm_size(incr = 1, p = 1, m = m_start):
                while p>pi2[-1]:
                    p_previous = p
                    if m < 15: m += 1
                    else: m += incr
                    print(m)
                    n = math.ceil(ratio*m)
                    r = np.zeros(k)
                    a = np.zeros(k)
                    w = np.zeros((B,k))
                    w_alt = np.zeros((B,k))
                    for i in np.arange(B):
                        # generate random observations
                        x_null = rpsi(n=(k,m), arg1=arg1, arg2=arg2)
                        y_null = rpsi(n=(k,n), arg1=arg1, arg2=arg2)
                        x_alt = rpsi(n=(k,m), arg1=arg1, arg2=arg2)
                        y_alt = np.reshape(rpsi(n=k*n, arg1= arg1 + delta*np.random.default_rng().choice(2,size=k*n,p=[1-theta,theta]), arg2=arg2), (k,n))
                        # rank data
                        if method=="sar":
                            w[i,] = np.cumsum(np.apply_along_axis(lambda x: sum(rankdata(x)[m:]), axis=1, arr= np.column_stack((x_null,y_null))))/k_range
                            w_alt[i,] = np.cumsum(np.apply_along_axis(lambda x: sum(rankdata(x)[m:]), axis=1, arr= np.column_stack((x_alt,y_alt))))/k_range
                        elif method=="rerank":
                            for i in k_range:
                                w[i,] = sum(rankdata(np.append(x_null[:i,],y_null[:i,]))[i*m:])
                                w_alt[i,] = sum(rankdata(np.append(x_alt[:i,],y_alt[:i,]))[i*m:])
                    # ts = test statistic
                    if method=="sar" and std:
                        # sar standardized
                        mu = n*(m*0.5 + (n+1)/2)
                        sig = math.sqrt(m*n*(1/2 - 0.5**2*((m+n)-1) + (n-1)*1/3 + (m-1)*1/3))
                        ts = (w-mu)/(sig/(k_range)**0.5)
                        ts_alt = (w_alt-mu)/(sig/(k_range)**0.5)
                    elif method=="sar" and not std:
                        # sar not standardized
                        ts = w
                        ts_alt = w_alt
                    elif method=="rerank" and std:
                        # rerank standardized
                        mu = k_range*n*(k_range*m*0.5 + (k_range*n+1)/2)
                        sig = (k_range**2*m*n*(1/2 - 0.5**2*(k_range*(m+n)-1) + (k_range*n-1)*1/3 + (k_range*m-1)*1/3))**0.5
                        ts = (w-mu)/sig
                        ts_alt = (w_alt-mu)/sig
                    elif method=="rerank" and not std:
                        # rerank not standardized
                        ts = w; ts_alt = w_alt

                    p_r = np.ones(k)
                    p_a = np.ones(k)
                    # get stage 1 upper critical value
                    p = 0
                    index = -1
                    u = sorted(np.unique(ts[:,0]), reverse=True)
                    while p <= pi1[0]:
                        p_rej = p
                        index += 1
                        p = sum(ts[:,0]>=u[index])/B
                        if index==(len(u)-1): print("break");break

                    if (index==0) & (sum(ts[:,0]>=u[index])/B>pi1[0]): # stopped on the first value
                        r[0] = float("inf")
                    elif sum(ts[:,0]>=u[index])/B<=pi1[0]:
                        r[0] = u[index]
                    else:
                        r[0] = u[index-1]
                    p_r[0] = sum(ts[:,0]>=r[0])/B

                    # get stage 1 lower critical value
                    p_acc = 0
                    u_alt = sorted(np.unique(ts_alt[:,0]))
                    index = -1
                    while p_acc <= pi2[0]:
                        index += 1
                        p_acc = sum(ts_alt[:,0]<=u_alt[index])/B
                    
                    if (index==0) & (sum(ts_alt[:,0]<=u_alt[index])/B>pi2[0]):
                        a[0] = float("-inf")
                    elif sum(ts_alt[:,0]<=u_alt[index])/B <= pi2[0]:
                        a[0] = u_alt[index]
                    else:
                        a[0] = u_alt[index-1]
                    p_a[0] = sum(ts_alt[:,0]<=a[0])/B
                    
                    # subset observations based on critical values of first stage
                    s = ts[(ts[:,0]<r[0]) & (ts[:,0]>a[0]),]
                    s_alt = ts_alt[(ts_alt[:,0]<r[0]) & (ts_alt[:,0]>a[0]),]
                    for i in np.arange(1,k):
                        p_rej = 0
                        u = sorted(np.unique(s[:,i]), reverse = True)
                        # print(u)
                        index = -1
                        while p_rej <= pi1[i]:
                            index += 1
                            p_rej = sum(s[:,i]>=u[index])/B
                            # if index==len(u): break
                        if (index==0) and (p_rej>pi1[i]):
                            r[i] = float("inf")
                        elif p_rej<=pi1[i]:
                            r[i] = u[index]
                        else:
                            r[i] = u[index-1]
                                
                        p_r[i] = sum(s[:,i]>=r[i])/B

                        p_acc = 0
                        u_alt = sorted(np.unique(s_alt[:,i]))
                        index = -1
                        while p_acc <= pi2[i]:
                            index += 1
                            p_acc = sum(s_alt[:,i] <= u_alt[index])/B
                            # if index==len(u_alt): break
                        if (index==0) and (p_acc > pi2[i]):
                            a[i] = float("-inf")
                        elif p_acc<=pi2[i]:
                            a[i] = u_alt[index]
                        else:
                            a[i] = u_alt[index-1]
                        p_a[i] = sum(s_alt[:,i]<=a[i])/B
                        if i!=np.arange(1,k)[-1]:
                            s = s[(s[:,i]<r[i]) & (s[:,i]>a[i]),]
                            s_alt = s_alt[(s_alt[:,i]<r[i]) & (s_alt[:,i]>a[i]),]
                    p = sum(s_alt[:,-1]<r[-1])/B
                    p_a[-1] = p
                return [m, n, np.around(r,4), np.around(a,4), np.around(np.append(p_r, sum(p_r)),4), np.around(np.append(p_a, sum(p_a)),4)]
            incr_size = 10
            res1 = incr_arm_size(incr = incr_size, p = 1, m = m_start)
            if res1[0]>=15:
                res = incr_arm_size(incr = 1, m = res1[0]-incr_size-1)
                d = dict();
                d['m'] = res[0]
                d['n'] = res[1]
                d['r'] = res[2]
                d['a'] = res[3]
                d['Type I'] = res[4]
                d['Type II'] = res[5]
                return d
            else:
                d = dict();
                d['m'] = res1[0]
                d['n'] = res1[1]
                d['r'] = res1[2]
                d['a'] = res1[3]
                d['Type I'] = res1[4]
                d['Type II'] = res1[5]
                return d
