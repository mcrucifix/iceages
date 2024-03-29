require(iceages)
data(models)
require(multitaper)
times = seq(-500,0,0.1)
deltat = c(0.001, 0.01, 0.1)

# uses 12 components of precession and obliquity
Astro <- read_astro(3,0)

# require(palinsol)
# B78O <- lapply(times*1e4, ber78)
# eccs <- as.numeric(sapply(B78O, function (lambda) lambda['ecc']))
# 
# B90O <- lapply(times*1e4, ber90)
# eccs90 <- as.numeric(sapply(B90O, function (lambda) lambda['ecc']))
# 
# 
# LARO <- lapply(times*1e4, la04)
# eccsla <- as.numeric(sapply(LARO, function (lambda) lambda['ecc']))

model = models$duffing_vdp_d
#model$spar['gamma'] = 0.20
#model$spar['omega'] = 1.0

sol = list()

kappas=seq(-1,2,0.1)

deltat = 0.1
sols <- lapply(kappas,
         function(kappa){
            spar <- model$spar
            spar['kappa'] <- kappa
            print(spar)
            propagate_d (model, spar, init=c(0.1,0.1), 
                             times=times, deltat=deltat, Astro=Astro)
})

i=0
# plot(times[-1],diff(sol[[1]][,2]), type='l', xlim=c(-20,-10), ylim=c(-4,4))
# lines(times[-1],diff(sol[[1]][,1]), type='l', xlim=c(-20,-10), ylim=c(-4,4))

plot(times/100,sols[[1]][,1], type='l', xlim=c(-1,0))
lapply(sols, function(sol) lines(times/100,sol[,1], type='l', col='red'))
#lines(times/100,sol[[3]][,1], type='l', col='blue')

spec <- lapply(sols, function(sol)
        Mod(fft(sol[,1]))[-1][0:1000]) 


freqs <- 1./(length(times)*deltat)*seq(1000)
periods <- 1./freqs * 10     # 1 time unit = 10ka


AmpMax <- sapply(spec, max) 
PerMax <- periods[sapply(spec, which.max)]




