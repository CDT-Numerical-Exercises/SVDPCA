# Chapter 5 – SVD and PCA 

The implementation of these problems had some interesting quirks due to how GSL handles SVD. Between the two problems, the data had to be represented in different ways, either row or column vectors, in order to satisfy GSL's SVD implementation. 

GSL's SVD implementation requires the matrix have at least as many rows as columns. For problem 1, we had 6D data and 2000 entries, so row vectors work here, producing a 2000x6 matrix. For problem 2, each image contains 19840 pixels and we have 10 images, so we need to use column vectors to ensure we have a 19840x10 matrix, instead of a 10x19840 matrix. 

Since the two representations are just transposes of each other, the SVD equations can be rearranged, and we find the result is that it swaps U and V; U contains the eigenvectors when using row vectors, V contains the eigenvectors when using column vectors. We just need to make sure we return the right one dependent on which format we are forced to use.

## Problem 1

```
== Best Eigenvectors ==
Component 1: [ -0.96841 -0.247092 0.0330351 0.0052721 -0.00255354 0.000994016 ]
Component 2: [ -0.247373 0.968836 -0.00365454 -0.00847422 -0.00250763 -0.00850569 ]
== Discarded Eigenvectors ==
Component 3: [ -0.0310953 -0.0118235 -0.99912 -0.0158303 -0.0183586 0.00806929 ]
Component 4: [ -0.00401035 -0.00745114 0.00416897 -0.858453 0.511955 -0.0295318 ]
Component 5: [ -0.00179713 0.00564518 -0.0243308 0.512535 0.856381 -0.0573834 ]
Component 6: [ -0.00111451 0.00870428 0.00673972 0.00411858 0.0645293 0.997846 ]
```

### Scree Plot

![scree plot](../assets/problem1_scree.png)

The scree plot allows us to justify reducing the dimensionality of the data down to two dimensions. The scree plot shows that the first two components account for ~85% of the total variance; since these two components explain almost all of the variance in the data, it is appropriate to reduce the dataset down to these two dimensions.

### Projected Data

![plot of data projected to 2D](../assets/problem1_projection.png)

## Problem 2

Images omitted for data protection reasons.
