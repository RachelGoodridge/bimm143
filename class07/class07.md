Class 07
================
Rachel Goodridge
January 29, 2019

### Functions revisited

``` r
source("http://tinyurl.com/rescale-R")
```

Lets try the **rescale()** function out

``` r
rescale(c(1,5,10))
```

    ## [1] 0.0000000 0.4444444 1.0000000

Lets try the **rescale2()** with the **stop()** function catch for non-numeric input

``` r
rescale2 <- function(x, na.rm=TRUE, plot=FALSE, ...) {
 if( !is.numeric(x) ) {
 stop("Input x should be numeric", call.=FALSE)
 }
 rng <-range(x, na.rm=na.rm)

 answer <- (x - rng[1]) / (rng[2] - rng[1])
 if(plot) {
 plot(answer, ...)
 }
 return(answer)
}
```

``` r
# rescale2(c(1:5,"string"))
```

Lets write a function to find "NA"

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

#positions of "NA"
which(is.na(x))
```

    ## [1] 3 5

``` r
which(is.na(y))
```

    ## [1] 1 3

``` r
#defining "NA" or not
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
#how many are "NA"
sum(is.na(x))
```

    ## [1] 2

``` r
sum(is.na(y))
```

    ## [1] 2

``` r
#how many are "NA" in the same place
sum(is.na(x)&is.na(y))
```

    ## [1] 1

Put it all together in a function now

``` r
both_na <- function(x, y) {
 sum( is.na(x) & is.na(y) )
}

both_na(x,y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)

both_na(x, y1)
```

    ## [1] 2

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

Lets write a new function to get rid of the unhelpful error and confusion

code: both\_na2 &lt;- function(x, y) { if(length(x) != length(y)) { stop("Input x and y should be vectors of the same length")} sum( is.na(x) & is.na(y) )

both\_na2(x, y1) both\_na2(x, y2)

Lets make a more comprehensive function now

``` r
both_na3 <- function(x, y) {
 if(length(x) != length(y)) {
 stop("Input x and y should be vectors of the same length")
 }

 na.in.both <- ( is.na(x) & is.na(y) )
 na.number <- sum(na.in.both)
 na.which <- which(na.in.both)
 message("Found ", na.number, " NA's at position(s):",
 paste(na.which, collapse=", ") )

 return( list(number=na.number, which=na.which) )
}

both_na3(x, y1)
```

    ## Found 2 NA's at position(s):2, 3

    ## $number
    ## [1] 2
    ## 
    ## $which
    ## [1] 2 3
