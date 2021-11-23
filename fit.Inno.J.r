load.pars = function() {

    h.table = '
1	-0.002	0.325	0.174	0.107	0.076	0.046	0.032	0.022	1.583	1.610	1.648	1.633	1.543	1.412	0.995	0.05
2	-0.001	0.373	0.172	0.092	0.051	0.032	0.017	0.012	1.460	1.652	1.702	1.728	1.814	1.657	1.661	0.05
3	0.003	0.418	0.156	0.077	0.028	0.014	0.003	0.007	1.423	1.823	1.903	1.843	1.310	0.951	1.921	0.04
4	0.004	0.430	0.134	0.083	0.016	0.012	0.012	0.005	1.355	2.011	2.114	1.261	2.292	1.135	1.849	0.04
5	-0.018	0.424	0.040	0.002	0.028	0.019	0.004	0.015	1.494	-3.112	0.199	0.404	2.548	0.643	1.396	0.04
6	0.000	0.412	0.076	0.028	0.021	0.011	0.011	0.005	1.616	2.225	0.787	0.811	0.952	-0.018	1.117	0.05
7	-0.002	0.486	0.093	0.063	0.063	0.045	0.045	0.024	1.421	2.070	1.581	1.535	1.591	1.910	2.081	0.05
8	-0.002	0.477	0.106	0.067	0.060	0.046	0.017	0.015	1.416	2.026	1.623	1.550	1.776	2.188	2.045	0.04
9	-0.003	0.424	0.114	0.067	0.039	0.030	0.012	0.012	1.218	1.731	2.331	2.652	-3.129	-2.609	-2.159	0.03
10	-0.004	0.413	0.133	0.067	0.047	0.033	0.016	0.015	1.215	1.802	2.115	2.409	2.390	2.912	-2.382	0.04
'
    con = textConnection(h.table)
    I15.pars <<- read.table(con)
    close(con)
    p.bins = '
1	1-3	12	3	3
2	3-5	40	35	44
3	5-7	26	29	29
4	7-9.5	18	23	29
5	9.5-10.5	11	13	13
6	10.5-13.5	19	21	25
7	13.5-15.5	24	22	23
8	15.5-20	17	19	18
9	20-30	14	17	19
10	30-100	16	19	17
'
    con = textConnection(p.bins)
    I15.bins = read.table(con)
    close(con)
    str.bins = as.character(I15.bins[,2])
    str.split = strsplit(str.bins,'-')
    I15.p.bins <<- matrix(NA,nrow=nrow(I15.pars),ncol=3)
    for (i in 1:nrow(I15.pars)) {
        I15.p.bins[i,1] <<- i
        I15.p.bins[i,2] <<- as.numeric(str.split[[i]][1])
        I15.p.bins[i,3] <<- as.numeric(str.split[[i]][2])
    }
}

get.amps = function(period) {
    idx = I15.p.bins[,2] < period & I15.p.bins[,3] > period
    if (sum(idx) != 1) {
        stop(paste0(' >> Error: Did not get the bin for period = ',period))
    }
    idx = I15.p.bins[idx,1]
    amps = I15.pars[idx,3:9]
    amps = as.numeric(as.character(amps))
    return(amps)
}

get.a0 = function(period) {
    idx = I15.p.bins[,2] < period & I15.p.bins[,3] > period
    idx = I15.p.bins[idx,1]
    a0 = I15.pars[idx,2]
    return(a0)
}

get.phis = function(period) {
    idx = I15.p.bins[,2] < period & I15.p.bins[,3] > period
    idx = I15.p.bins[idx,1]
    phis = I15.pars[idx,10:16]
    phis = as.numeric(as.character(phis))
    return(phis)
}

calt0 = function(x,PHI,a0) {
    t0 = a0
    for (i in 1:7)
        t0 = t0 + amps[i]*cos(2*pi*i*(x+PHI) + phis[i])
    return(t0)
}

calt = function(x,PHI,M,L,a0) {
    t = M + L*calt0(x,PHI,a0)
    return(t)
}

fit.Inno15 = function(phase, mag, err,
    ini.P = -99, ini.L = -99, ini.M = -99) {
    amps <<- get.amps(period)
    a0 <<- get.a0(period)
    phis <<- get.phis(period)
    nobs = length(phase)
    Dmat = matrix(1, nrow=nobs, ncol=2)
    Y = matrix(mag,nrow=nobs,ncol=1)
    SIGMAS = err / mean(err)
    S = matrix(0, nrow=nobs, ncol=nobs)
    for (idiag in 1:nobs) {
        S[idiag,idiag] = SIGMAS[idiag]^2
    }
    iS = solve(S)
    if (ini.P == -99) {
        nPHIs = 1000
        PHIs = seq(from = 0, to = 1, length=nPHIs)
        chisqrs = rep(NA,nPHIs)
        Ms = rep(NA,nPHIs)
        Ls = rep(NA,nPHIs)
        for (iPHI in 1:nPHIs) {
            PHI = PHIs[iPHI]
            X = phase
            Dmat[,2] = calt0(X,PHI,a0)
            tD = t(Dmat)
            beta = solve(tD %*% iS %*% Dmat) %*% (tD %*% iS %*% Y)
            chisqrs[iPHI] = sum((mag - calt(X,PHI,beta[1],beta[2],a0))^2/err^2)
            Ms[iPHI] = beta[1]
            Ls[iPHI] = beta[2]
        }
        ## plot(PHIs,Ls)
        ## print(beta)
        ## abline(v=PHIs[which.min(chisqrs)])
        ## As expected, there might be some best-fit with L<0, which should be taken care of:
        ## L' = -L, M' = M, and PHI' = PHI + 0.5 and re-test
        mid = cbind(Ms, Ls, PHIs, chisqrs)
        subp = mid[mid[,2]>0,]
        subn = mid[mid[,2]<0,]

        if (min(mid[,4]) == min(subp[,4])) {
            idx = which.min(subp[,4])
            M_guess = subp[idx,1]
            L_guess = subp[idx,2]
            PHI_guess = subp[idx,3]
        } else {
            idx = which.min(subp[,4])
            M_p = subp[idx,1]
            L_p = subp[idx,2]
            PHI_p = subp[idx,3]
            chisqr_p = subp[idx,4]
            idx = which.min(subn[,4])
            M_test = subn[idx,1]
            dM_test = sd(mid[,1])
            L_test = -1 * subn[idx,2]
            dL_test = 0.1
            PHI_test = subn[idx,3] + 0.5
            dPHI_test = 0.01
            if (PHI_test > 1) PHI_test = PHI_test - 1
            M_range = seq(from = M_test-dM_test, to = M_test+dM_test, length = 100)
            L_range = seq(from = max(0,L_test-dL_test), to = L_test+dL_test, length = 100)
            PHI_range = seq(from = PHI_test-dPHI_test, to = PHI_test+dPHI_test, length = 20)
            tmp2 = matrix(NA,ncol=4,nrow=length(M_range)*length(L_range)*length(PHI_range))
            itmp2 = 1
            for (mtest in M_range) {
                for (ltest in L_range) {
                    for (phitest in PHI_range) {
                        chisqr = sum((mag - calt(X,phitest,mtest,ltest,a0))^2/err^2)
                        tmp2[itmp2,1] = mtest
                        tmp2[itmp2,2] = ltest
                        tmp2[itmp2,3] = phitest
                        tmp2[itmp2,4] = chisqr
                        itmp2 = itmp2 + 1
                    }
                }
            }
            idx = which.min(tmp2[,4])
            chisqr_n = tmp2[idx,4]
            if (chisqr_n < chisqr_p) {
                M_guess = tmp2[idx,1]
                L_guess = tmp2[idx,2]
                PHI_guess = tmp2[idx,3]
            } else {
                idx = which.min(subp[,4])
                M_guess = subp[idx,1]
                L_guess = subp[idx,2]
                PHI_guess = subp[idx,3]
            }
        }
        y = mag
        x = phase
        weights = 1/SIGMAS^2
        ## plot(x, y)
        ## xc = seq(0,1,1e-2)
        ## yc = calt(xc, PHI_guess, M_guess, L_guess,a0)
        ## lines(xc, yc)
        
        fit = nls(
            y ~ calt(x, par_phi, par_m, par_l, a0),
            start = list(
                par_phi=PHI_guess+1e-6,
                par_m = M_guess,
                par_l = L_guess),
            weights=weights,
            control = nls.control(minFactor=1/102400))
        nlspars = summary(fit)$parameters[,1:2]
        nls_M = nlspars['par_m',1]
        nls_L = nlspars['par_l',1]
        nls_PHI = nlspars['par_phi',1]
        nls_eM = nlspars['par_m',2]
        nls_eL = nlspars['par_l',2]
        nls_ePHI = nlspars['par_phi',2]
        ret = c(nls_M, nls_L, nls_PHI, nls_eM, nls_eL, nls_ePHI)
        return(ret)
    }
}
