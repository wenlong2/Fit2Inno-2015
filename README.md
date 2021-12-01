# Fit2Inno2015
Fit Inno+2015 NIR Cepheid templates to light curve measurements

Example:

```R
source('fit.Inno.H.r')
load.pars()
f.mlc = mdir + id + '_ukt.dat'
dat = read.table(f.mlc)
ylab = bquote(italic(.(band))~' [mag]')
t = dat[,1] - t0
y = dat[,2]
e = dat[,3]
Q = q0 ## Cepheid period or other folding rules
x = (t %% Q) / Q
col = rep(1, nrow(dat))
plot(x, y, xlab='Phase',ylab=ylab, ylim=rev(range(y))+c(0.15,-0.1)/2, main=id, xlim=c(0,1), pch=19, col=col)
arrows(x, y-e, x, y+e, length=0, col=col)
pars = -1
try(pars <- fit.Inno15(x, y, e))
if (pars[1] == -1) stop(paste(id, 'failed to fit Inno+ 15'))
M = pars[1]
L = pars[2]
PHI = pars[3]
xc = seq(0, 1, 0.01)
yc = calt(xc, PHI, M, L, a0)
res = calt(x, PHI, M, L, a0) - y
lines(xc, yc, col=2, lwd=1)
```
