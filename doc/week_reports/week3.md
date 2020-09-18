### Week report 3

I worked with the O(N log^2 N) building algorithm and tree traversal.
The kd tree implementation seems functional now. I also tested
kd tree queries against a brute force implementation. 

I implemented quicksort to replace the standard library std::sort.
(This was required for the course)
I also implemented the quickselect with median of medians to be used
as a pivot selection strategy for the quicksort and to ensure
O(N log N) worst case time complexity. Unfortunately this turned
out to be even slower than expected (~7 x slower than std::sort for large
random arrays) so
currenlty "median of three" pivot selection strategy is used
(still ~2.5 x slower than std::sort for large random arrays).

After this, I tried to optimize the kd tree traversal by replacing pointers
to Triangle objects by copying Triangle objects to the kd tree. To better
understand if this was helpful I would need to compare the kd tree
performance on real triangle data.

I also properly learnt the O(N log N) algorithm for the surface area heuristic kd tree
construction but didn't have time to implement it yet. 
Also it should be noted that the "O(N log N)" algorithm is
not always O(N log N).

Time: 20 h ?

## Future

I will think more about the assumptions required to make the
O(N log N) algorithm really O(N log N).

I will implement some sort of interface to test the algorithms with 
real 3d models. This might take some time but if I have time I will finally start
implementing the O(N log N) tree building algorithm.

At some point (not next week) I will also try to optimize the quicksort implementation.
