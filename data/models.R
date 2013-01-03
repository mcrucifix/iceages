# Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject

# the following conditions:

# The above copyright notice and this permission notice shall be
# incluudedin all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR

# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# ------------------------------------------------------------------
# R Code developed for R version 2.15.2 (2012-10-26) -- "Trick or Treat"
# ------------------------------------------------------------------ 


# --------------------------------------
# OVERHEAD

models <- list( 

vdp_d = list (func='vdp_f', name = 'vanderpol' , 
               spar = c(ALPHA=30.0, BETA=0.7, GAMMAPRE = 0.7,
                       GAMMAOBL = 0.7, omega = 3.2, ASYM= 0.0),
               initgrid = list(x=seq(-1.2,1.2,0.2), y=seq(-1.2,1.2,0.2))) ,

vdp_s = list (func='vdp_s', name = 'vanderpol' , 
               spar = c(ALPHA=30.0, BETA=0.7, GAMMAPRE = 0.7,
                       GAMMAOBL = 0.7, omega = 3.2, SIGMA= 0.2),
               initgrid = list(x=seq(-1.2,1.2,0.2), y=seq(-1.2,1.2,0.2))) ,


s90_d = list (func='s90_f', name = 'Saltzman - Maasch 1990' , 
               spar = c(gamma = 1.0, omega = 1.0, w=0.5) ,
               initgrid = list(x=seq(-1.0,1.0,0.5), y=seq(-1.0,1.0,0.5), z=seq(-1.0,1.0,0.5))),

s91_d = list (func='s91_f', name = 'Saltzman - Maasch 1991' , 
               spar = c(gamma = 1.0, omega = 1.0) ,
               initgrid = list(x=seq(-1.0,1.0,0.5), y=seq(-1.0,1.0,0.5), z=seq(-1.0,1.0,0.5))),

pp04_d = list (func='pp04_f', name = 'Saltzman - Maasch 1991' , 
               spar = c(gamma = 1.0, omega = 1.0, alpha=0.15, d=0.27, aa=0.4, bb=0.7) ,
               initgrid = list(x=seq(-1.0,1.0,0.5), y=seq(-1.0,1.0,0.5), z=seq(-1.0,1.0,0.5))),


t06_d = list (func='t06_f', name = 'Tziperman et al. 2006' , 
               spar = c(gamma = 1.0, omega = 1.0, p0=0.26, s=0.23, sm=0.03) ,
               initgrid = list(x=seq(0,60,2), y=c(0,1))),


i80_d = list (func='i80_f', name = 'Imbrie - Imbrie 1980' , 
               spar = c(gamma = 1.0, omega = 1.0) ,
               initgrid = seq(-1,1,0.1)) ,


i11_d = list (func='i11_f', name = 'Imbrie et al. 2011' , 
               spar = c(gamma = 1.0, omega = 1.0) ,
               initgrid = list(x=seq(-1.2,1.2,0.2), y=seq(-1.2,1.2,0.2))) ,

i11_altered_d = list (func='i11_altered_f', name = 'Imbrie et al. _altered_ 2011' , 
               spar = c(gamma = 1.0, omega = 1.0) ,
               initgrid = list(x=seq(-1.2,1.2,0.2), y=seq(-1.2,1.2,0.2))) ,

pp12_d = list (func='pp12_f', name = 'Parrenin and Paillard, 2012' , 
               spar = c(gamma = 1.0, omega = 1.0) ,
               initgrid = list(x=seq(0,120,3), y=c(0,1))), 

cr12_d = list (func='cr12_f', name = 'Crucifix et al. 2012' , 
               spar = c(alpha = 30., beta0 = 0.1894, beta1 = 0, 
                        beta2 = 0., delta = 0.2706, 
                        gammapre = 0.1894, gammaobl = 0.1624,
                        omega = 1.000), 
                        initgrid = list(x=seq(-1.5,1.5,0.5), y=seq(-1.5, 1.5, 0.5))),

cr12_s = list (func='cr12_s', name = 'Crucifix et al. 2012' , 
               spar = c(alpha = 30., beta0 = 0.1894, beta1 = 0, 
                        beta2 = 0., delta = 0.2706, 
                        gammapre = 0.1894, gammaobl = 0.1624,
                        omega = 1.000, sigmax=0., sigmay=0.), 
                        initgrid = list(x=seq(-1.5,1.5,0.5), y=seq(-1.5, 1.5, 0.5)))


)
