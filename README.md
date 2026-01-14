1. Cluster has been detected by breadth-fast algorithm
2. We extract the largest cluster Smax and its corresponding coordinates.

To identify if Smax percolates or not,
3. We take the simulation box(original_box) that contains the largest cluster only.
4. Replicate the box in all three directions, so now the box size is now 3*L
5. We now have a new_box that is 27*origianl_box.
6. We do a cluster computation in the new_box again with the same algorithm.
7. a) If Smax spanns in one direction we end up detecting a cluster in new_box with newSmax= 3*Smax.
   b) Similarly if it is in two directions newSmax=9*Smax
   c) If all three directions, it will be newSmax=27*Smax

** I think about extracting the direction of percolation. The way I did above, it is difficult to identify the direction. However, if we do the replication one direction at a time and check wheather it spans, we can identify that. In that case the code needs to modify a LOT as the boxlength will be 3L,L,and L if we replicate in x-direction and so on.

** If it is important to identify the direction, I can give a try!


 Breadth-fast algorith (It is available everuwhere, still I am trying to explain...)

 1. Start from a atom, find its neighbouring atom. Keep them in a list.
 2. Go over the list, for each of them(element in the list) find their neighboring atom.
 3. Put the neighbours in the list if they are not added in the list before.
 4. Update the list so that all the elements are new in the list.
 5. Goto (2.) do agian until there is no new element in the list.
 --- One cluster will be detected.
 While doing so, we can mark the atom with some identity. Next scan for next atom, if it is not marked start from (1.) to find the next cluster and so on.
 
