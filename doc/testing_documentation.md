# Testing documentation

## Unit tests
Testing was done with the gtest framework. 

## Performance report

The performance report can be run with `make performance` command.

### Quicksort

The performance of different quicksort implementations was
compared by generating arrays of random integers. The integers
for the array of size n was sampled from 0-n.

The median of medians quicksort (`mmQuickSort`) was about 4 times slower than
the median of three quicksort (`quickSort`). This is why elsewhere in the
program only the median of three version is used.

The tests were repeated 100 times.

Summary:
```
n = 100
std::qsort:              0.00011329ms
std::sort:               5.132e-05ms
quickSort:               4.773e-05ms
mmQuickSort:             0.00010164ms

n = 100000
std::qsort:              0.167069ms
std::sort:               0.0669119ms
quickSort:               0.0816432ms
mmQuickSort:             0.301048ms

n = 10000000
std::qsort:              17.2194ms
std::sort:               7.82591ms
quickSort:               10.438ms
mmQuickSort:             44.7186ms
```

While the `quickSort` was faster than the C qsort,
it was still slightly slower than the C++ standard library sort.
(which is, unlike median of three `quickSort`, guaranteed to have
O(n log n) time complexity) 

### Kd tree
#### Structure
The structure of the generated kd tree was compared with the
results of I. Wald, et al., (2006). 

The testing was only done with the "bunny.obj" model. (As a note,
the armadillo.obj model that I found had different number of triangles
compared with the armadillo object that the authors of the paper had).

The results were:
```
                     I. Wald, et al.         My implementation
number of leaves     338k                    404k
non-empty leaves     159k                    163k
cost                 926                     967.135
```

I was not able to determine if these different results were
due to a bug in my code or differences in the implementation details.

However, the cost function (estimated number of work per one query) was
quite close to the one by I. Wald et al. I also noticed that
by making some perturbations to my code, such as limiting the
splits to only medians, the cost rose to over 4000. Unfortunately
I didn't have time to properly test this.

I also tested the code by changing the EPS constant from 1e-10 to 1-15 in geometry.h
(Unfortunatley this requires modifying the constant and
is currently not included in `make performance`) and got the following
results:

```
(EPS = 1e-15)
number of leaves     462k
non-empty leaves     188k
cost                 967.546
```

With EPS = 1e-8, the results were:
```
(EPS = 1e-8)
number of leaves    376k
non-empty leaves    151k
cost                967.136
```

Also I wasn't completely sure if the equation 3 in the paper was correct
(or if I read it correctly).

I interpreted the footnote "In theory, adding a constant to simulate the setup cost for traversing the
leaf (i.e., a leaf cost of 1KT +NKI ) should be more accurate"
that the leafs should only contribute to the intersection cost and
that their traveral cost should be 0.

I also interpreted "That the cost of intersecting N triangles is roughly NKI , i.e.,
directly linear in the number of triangles" that
the cost of a leaf L should be: `SA(L)/SA(scene) * k_i * N(L)` where
`N(L)` is the number of triangles in the leaf.

Thus the cost function that I calculated was:
```
sum_(x \in inner_nodes)(k_t * SA(x)/SA(scene)) + sum_(y \in leaves)(k_i * SA(y)/SA(scene) * N(y))
```

#### Random queries from a model
Random intersection queries were done with models/armadillo/armadillo.obj.
The model contains about 200 000 triangles.

The starting point of the ray was sampled randomly within the
axis-aligned bounding box of the object and the direction of the
ray was sampled randomly.

The test was repeated 10 times.
Results:
```
Kd tree time:       0.32054ms
Naive time:         402.698ms
```

#### Rendering
For rendering, the performance was tested by rendering
a models/bunny/bunny\_hires.obj The model contains about 69 000 triangles.

Excluding the time used to build the kd tree,
the rendering time was about 0.05 s with the kd tree
and 21 s without it. This is about 400 x improvement compared with the
naive implementation.


# References

I. Wald and V. Havran, "On building fast kd-Trees for Ray Tracing, and on doing that in O(N log N)," 2006 IEEE Symposium on Interactive Ray Tracing, Salt Lake City, UT, 2006, pp. 61-69, doi: 10.1109/RT.2006.280216.
