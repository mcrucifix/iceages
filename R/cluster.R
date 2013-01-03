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


### simple clustering using the distance function supplied by R 

cluster = function(X, cminratio=0.01, tol=0.01)
{
if (is.vector(X)) X = matrix(X)
N = nrow(X)
cmin      = floor(cminratio * N)
d = dist(X)

cluster = rep(NA,N)

nc = 0
ii = 1

for (ij in 1:N)
{

  if (is.na(cluster[ij] )) 
   { nc = nc+1
cluster[ij] <- nc
   }

  P = which(d[ii:(ii+N-ij-1)] < tol)
  if (ij < N) cluster[P+ij] <- cluster[ij]

  ii = ii + N - ij
}
## CC will contain the number of clusters for each time step
CC = sapply(1:nc, function(i) length(which(cluster==i)))

keep = which(CC > cmin)

CCkeep = CC[keep]

nck=length(keep)
ClusterCenters = t(sapply(1:nck, function(i) {
     apply(X[which(cluster == keep[i]),,drop=FALSE],2,mean) }))

cluster_reordered = cluster
for (i in 1:nck) cluster_reordered[which(cluster_reordered == keep[i])] = i

if (ncol(X) == 1) ClusterCenters = t(ClusterCenters)

rm('d') # make sure to release memory used by 'd'

list(nc=nc,nck=nck,belong=cluster_reordered,kernel=ClusterCenters)

}


