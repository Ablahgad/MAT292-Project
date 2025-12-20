# MAT292 Project: CFD Solver - Ciara Browne, Nancy Fang, and Ablah Gadallah (Group 21)
A CFD solver for modelling the motion of a Zebrafish in water using a finite difference numerical model. Experimental data used for comparison sourced from ["Zebrafish swimming in the flow: a particle image velocimetry study" by Mwaffo et al](https://pmc.ncbi.nlm.nih.gov/articles/PMC5691796/).

## Instructions
1. Download repository and ensure essential libraries (scipy, numpy) are downloaded, change the working directory to the repository
2. Settings for the finite difference solver can be changed from the [Finite_difference.py file](Finite_difference.py)
3. Run solver or precalculated (flow speed equal to zero or 26mm/s) examples from the [visualizer.mlx file](visualizer.mlx)
4. Run percent difference calculation from the [Percent_Difference.py file](Percent_Difference.py) to compare calculated value with experimental PIV values. **Ensure divisions_x = 25 and divisions_y = 21**
5. Visualize percent difference in the respective section of the  [visualizer.mlx file](visualizer.mlx)

MATLAB file will also output GIF of outputted vector fields as [nodes.gif](nodes.gif). Please be warned that this solver is prone to blowing up when the time duration or time step are too large.
