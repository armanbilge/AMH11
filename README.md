# AMH11

A Java implementation of the matrix exponential method described by Al-Mohy and
Higham (2011), which computes the action of a matrix exponential on a vector
without explicitly forming the matrix exponential.

> Al-Mohy, A. H., & Higham, N. J. (2011). Computing the action of the matrix
> exponential, with an application to exponential integrators. *SIAM Journal on
> Scientific Computing, 33*(2), 488–511.
> doi:[10.1137/100788860](http://dx.doi.org/10.1137/100788860)

This implementation is derived from their original
[MATLAB codes](http://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector)
and uses the [MTJ library](https://github.com/fommil/matrix-toolkits-java).
