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

read_astro <-
function(nap=1, nao=1)
 { 
 pfile <- file.path(path.package("iceages"),"orig","precession.dat")
 ofile <- file.path(path.package("iceages"),"orig","obliquity.dat")
 Astro <-  .Call("c_readinsol", nap, nao, pfile, ofile)
 }
