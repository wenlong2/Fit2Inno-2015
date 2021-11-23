load.pars = function() {

    h.table = '
1	-0.013	0.386	0.135	0.057	0.014	0.009	0.004	0.003	1.338	2.020	2.381	2.423	-3.111	2.000	1.658	0.05
2	-0.012	0.388	0.135	0.057	0.014	0.008	0.003	0.002	1.332	2.017	2.406	2.435	-3.096	1.943	1.614	0.05
3	-0.003	0.432	0.120	0.066	0.020	0.007	0.005	0.006	1.291	2.104	2.702	2.437	2.389	-2.527	2.398	0.04
4	-0.003	0.424	0.111	0.066	0.009	0.004	0.002	0.003	1.253	2.241	2.900	2.695	-3.061	1.498	3.054	0.04
5	0.003	0.462	0.069	0.016	0.027	0.023	0.017	0.012	1.456	3.058	0.011	2.419	-1.323	1.461	3.091	0.05
6	-0.002	0.480	0.072	0.019	0.017	0.006	0.006	0.004	1.396	2.878	2.012	2.754	2.726	1.223	1.922	0.03
7	0.000	0.485	0.056	0.031	0.013	0.010	0.006	0.004	1.415	2.667	2.528	2.478	-2.969	-2.867	-1.826	0.03
8	-0.001	0.476	0.064	0.039	0.013	0.007	0.010	0.005	1.369	2.498	2.920	3.009	-2.705	-1.954	-0.790	0.03
9	-0.003	0.451	0.094	0.046	0.023	0.013	0.008	0.001	1.330	2.328	-3.079	-2.651	-1.692	-0.741	1.184	0.03
10	-0.006	0.451	0.092	0.035	0.011	0.006	0.003	0.001	1.314	2.505	-2.720	-1.883	-0.768	1.090	-2.898	0.03
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
