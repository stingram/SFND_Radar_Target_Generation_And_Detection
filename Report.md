# Radar Target Generation and Detection

## Implementation steps for 2D CFAR

I implemented a 2D CFAR on the range doppler output using the following steps:
1. Determined number of training cells for each dimension Tr and Td. I similarly picked the number of guard cells Gr and Gd.
2. I slid the Cell Under Test (CUT) across the complete RDM.
3. I selected the grid that included the training, guard, and test cells.
4. I selected the smaller grid that only included the guard and test cells.
5. The values in these selected were converted from logarithmic to linear using the db2pow function.
6. I summed up all values from the larger grid and summed up all the values from the smaller grid.
7. I subtracted the sum of the smaller grid from the sum of the bigger grid.
8. The average was then computed by dividing the difference value by the total number of training cells around the CUT.
9. The value was then converted back to logarithmic using pow2db function.
10. I then added the offset to determine the threshold for this CUT.
11. I compare the signal under CUT against this threshold if CUT level > threshold then the RDM output was assigned a value of 1, else it was set to 0.

## Selection of Training, Guard cells and offset.

For selection of the number of Training cells, Guard cells, and the offset values I referred to past lessons and exercises. For each side of the CUT I chose: Tr = 12, Td= 12, Gr = 5, Gd =  3. For the offset I chose a value of 6. 

## Steps taken to suppress the non-thresholded cells at the edges.

Suppression of the non-thresholded cells at the edges was done by setting all cells in the edge region of the RDM to zero. This region was where the cell row index was a either <  (Tr+Gr+1) or > Nr-(Tr+Gr) while the column index was either < (Td+Gd+1) or n>Nd-(Td+Gd).
