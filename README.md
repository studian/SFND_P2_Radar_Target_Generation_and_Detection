# SFND_P2_Radar_Target_Generation_and_Detection
Udacity, Sensor Fusion, Project of Radar Target Generation and Detection

## Project Layout:
<img src="./ProjectLayout.JPG" width="700" />

## Implementation steps for the 2D CFAR process
* Loop over elements of RDM array each iteration selecting one cell to be the CUT (Cell Under Test)<br>
`for i = Tr+Gr+1 : (Nr/2)-(Gr+Tr)`<br>
`for j = Td+Gd+1 : Nd-(Gd+Td)`
* For each iteration loop over the training cells "excluding the guarding cells" to sum their values<br>
`for p = i-(Tr+Gr) : i+(Tr+Gr)`<br>
`for q = j-(Td+Gd) : j+(Td+Gd)`
* Calculate the average of the noise value<br>
`noise_level = noise_level + db2pow(RDM(p,q));`
* Convert using pow2db<br>
`th = pow2db(noise_level/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));`
* Add the offset value
* If the CUT is greater then the threshold replace it by `1`, else `0` <br>
and that’s all for the Implementation.

## Selection of Training, Guard cells and offset
* `Tr = 10, Td = 8` For both Range and Doppler Training Cells.
* `Gr = 4, Gd = 4` For both Range and Doppler Guard Cells.
* `offset = 1.4` the offset value.

## Steps taken to suppress the non-thresholded cells at the edges
This was done throught sclicing the output such that we have the surrounding rows and columns depending on the Training cells for both range and doppler.<br>
`RDM(union(1:(Tr+Gr),end-(Tr+Gr-1):end),:) = 0;  % Rows`<br>
`RDM(:,union(1:(Td+Gd),end-(Td+Gd-1):end)) = 0;  % Columns`

## Output:
<img src="results/Fig1_Range_from_First_FFT.bmp" width="700" />
<br>
<img src="results/Fig2_Range_and_Speed_From_FFT2.bmp" width="700" />
<br>
<img src="results/Fig3_CA-CFAR_Filtered_RDM.bmp" width="700" />
