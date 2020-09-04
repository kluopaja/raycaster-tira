# Design document

A ray caster with an acceleration scheme.

Ray casting or tracing (recursive) requires solving the following problem:
Given a set of triangles T and a ray r, find the first intersection
of ray r with any of the triangles in T.

Naively this can be done by simply going though every triangle t in T,
finding the intersection section point between t and r and finally 
picking the nearest of these intersection points.

In practice, this can be sped up with various data structures.

## Algorithm design and data structures

The primary data structure implemented here will be the k-d tree.
This was chosen based on the literature [1]. To provide some comparison,
a non-accelerated version as well as an octree will be implemented.

k-d trees are k dimensional binary trees. Every node of a k-d tree corresponds
to some k-dimensional
rectangle. Every non-leaf node also divides the corresponding rectangle
into two rectangles which recursively correspond to the children nodes. [2]

The implementation of building of the k-d tree can vary depending on the
exact way each rectangle is divided into subtrees and when the recursion is
stopped. [1] Different ways to do this will be compared.

The surface area heuristic will be implemented using an O(N log N) algorithm 
desribed by Wald and Havran [1].

Some tree traversal strategy will be implemented to answer the queries.
Many algorithms have been proposed to improve the practical efficiency of
this step.[3] It would be interesting to implement couple of these.

It should be noted that k-d trees do not improve the worst case time complexity
of finding the closest intersection for a ray.

## Inputs and outputs

First the program will receive a list of triangles and the number of 
queries. The program will also receive details about used
heuristics/data structures.

Every query will specify a ray. After receiving a query, the program
will immediately start to process it and output the result as soon as
the program has finished processing it.

## Language:
English,
C++

## Opinto-ohjelma:
Molekyylibiotieteiden kandiohjelma

## References
[1] Wald, Ingo, and Vlastimil Havran. "On building fast kd-trees for ray tracing, and on doing that in O (N log N)." 2006 IEEE Symposium on Interactive Ray Tracing. IEEE, 2006.

[2] Bentley, Jon Louis. "Multidimensional binary search trees used for associative searching." Communications of the ACM 18.9 (1975): 509-517.

[3] Hapala, Michal, and Vlastimil Havran. "Kd‐tree traversal algorithms for ray tracing." Computer Graphics Forum. Vol. 30. No. 1. Oxford, UK: Blackwell Publishing Ltd, 2011.
