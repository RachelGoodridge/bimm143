Class 6
================
Rachel Goodridge
January 24, 2019

### Section 1: Reading files again :-)

Here we are going to try using the **read.table()** function and friends to read some example flat files.

First lets try adding a simple plot

``` r
plot(1:10)
```

![](class06_files/figure-markdown_github/unnamed-chunk-1-1.png)

Back to file reading...

``` r
read.table("https://bioboot.github.io/bimm143_W19/class-material/test1.txt",header=TRUE,sep=",")
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

For this common CSV format we can use **read.csv()**

``` r
file1 <- "https://bioboot.github.io/bimm143_W19/class-material/test1.txt"
data1 <- read.csv(file1)
data1
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
file2 <- "https://bioboot.github.io/bimm143_W19/class-material/test2.txt"
data2 <- read.table(file2,header=TRUE,sep="$")
data2
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
file3 <- "https://bioboot.github.io/bimm143_W19/class-material/test3.txt"
data3 <- read.table(file3)
data3
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

### Section 2: R Functions

My first silly function

``` r
add <- function(x, y=1) {
 # Sum the input x and y
 x + y
}
```

Lets try using this function

``` r
add(1)
```

    ## [1] 2

``` r
add(1,100)
```

    ## [1] 101

``` r
add( c(1, 2, 3) )
```

    ## [1] 2 3 4

``` r
add(c(1, 2, 3),4)
```

    ## [1] 5 6 7

``` r
#add(1, 2, 2)
```

``` r
#add(x=1, y=“b”)
```

When you find yourself doing the same thing three or more times, make a function instead.

Instead of this:

``` r
#df$a <- (df$a - min(df$a)) / (max(df$a) - min(df$a))
#df$b <- (df$b - min(df$a)) / (max(df$b) - min(df$b))
#df$c <- (df$c - min(df$c)) / (max(df$c) - min(df$c))
#df$d <- (df$d - min(df$d)) / (max(df$a) - min(df$d))
```

Make a function:

``` r
rescale <- function(x) {
  rng <- range(x)
  (x-rng[1])/(rng[2]-rng[1])
}
```

Test on a small example where we know what the answer should be!

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
#How would you get this function to run?
#rescale(c(1,2,NA,3,10))
```

``` r
#What should your function do here?
#rescale(c(1,10,"string"))
```

Lets make a new version of the function to get rid of the NA error.

``` r
rescale2 <- function(x,na.rm=TRUE) {
  rng <- range(x,na.rm=na.rm)
  (x-rng[1])/(rng[2]-rng[1])
}
```

``` r
rescale2(c(1,2,NA,3,10))
```

    ## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000

Version 3 of the same rescale function with plot

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
 rng <-range(x, na.rm=na.rm)
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
 return(answer)
}
```

``` r
rescale3(c(1:6,NA,10),plot=TRUE)
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"

![](class06_files/figure-markdown_github/unnamed-chunk-21-1.png)

    ## [1] "I can see it in ..."

    ## [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556        NA
    ## [8] 1.0000000

### Working with the bio3d package

To install this package I used the command **install.packages("bio3d")**

Lets get using it.

``` r
library(bio3d)
```

``` r
# Read a PDB file from the database
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1
```

    ## 
    ##  Call:  read.pdb(file = "4AKE")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 3459,  XYZs#: 10377  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 3312  (residues/Calpha atoms#: 428)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 147  (residues: 147)
    ##      Non-protein/nucleic resid values: [ HOH (147) ]
    ## 
    ##    Protein sequence:
    ##       MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT
    ##       DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI
    ##       VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG
    ##       YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILGMRIILLGAPGA...<cut>...KILG
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

Try out the provided code...

``` r
s1 <- read.pdb("4AKE") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

    ## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): C:
    ## \Users\Rachel\AppData\Local\Temp\RtmpSyObJh/4AKE.pdb exists. Skipping
    ## download

``` r
s2 <- read.pdb("1AKE") # kinase no drug
```

    ##   Note: Accessing on-line PDB file
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
s3 <- read.pdb("1E4Y") # kinase with drug
```

    ##   Note: Accessing on-line PDB file

``` r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
plotb3(s1.b, sse=s1.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
plotb3(s2.b, sse=s2.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-24-2.png)

``` r
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
```

![](class06_files/figure-markdown_github/unnamed-chunk-24-3.png)
