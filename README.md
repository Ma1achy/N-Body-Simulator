
=====N-BODY SIMULATIONS=====

A collection of Python scripts to simulate the gravitational interactions between massive bodies using either the naive brute force algorithm
or the Barnes-Hut algorithm. This was written for a project for a scientific computing module, done as part of an undergraduate physics degree.

=====GETTING STARTED=====

To create a simulation open "generate-templates.py" and run the script. It will create 8 simulation templates by default, these are .csv files containing
the inital conditions to different simulations, the default directory is N-Body/sim templates.

To run a simulation open either "nbody-barnes-hut.py" or "nbody-brute-force.py" and change the path variable to contain the path of the simulation template 
you want to run. For ""nbody-brute-force.py" the path variable is "input_path".

For "nbody-barnes-hut.py" it is recommended to use the command line argument '--input' or '-i' to change the path variable, but you can also directly edit the "default"
variable if required. By default both alogrithms will use their test case templates. NOTE: The simulation templates are not compatiable between the two different algorithms.
(The barnes-hut file can be ran through command line arguments, I forgot to implement it with the rest of the scripts though whoops)

Once the simulation has finished processing to visualise it, find the associated simulation result "YYYY-MM-DD HH-MM-SS sim" folder in "N-Body/results". Copy the path to the 
"YYYY-MM-DD HH-MM-SS sim" folder and in "visualise.py" change the "path" variable to the path of the "YYYY-MM-DD HH-MM-SS sim" folder and run the script. 
A Pygame window should open and you will be able to see the simulation results animate. The default visualisation is set to a precomputed simulation of the test case. 
By default the visualisation saves every frame to a "YYYY-MM-DD HH-MM-SS render" folder in "N-Body/rendered"

To render a simulation to a .gif file copy the path to the associated "YYYY-MM-DD HH-MM-SS render" folder and in "render-to-gif.py" replace the "render_dir" variable with the
path to the render folder, then run the script. Inside the "YYYY-MM-DD HH-MM-SS render" folder a gif will be generated from the frames of the visualised simulation with the name
"YYYY-MM-DD HH-MM-SS render.gif".

There are two precomputed simulations, The test case, and a simulation of the solar system ready to be used with the Pygame visualisation.

=====DEPENDENCIES=====

Python v3.12+

numpy
scipy
matplotlib
plotly
astropy.time
astropy.coordinates
Pygame
imageio.v2

=====HELP=====

Why is my simulation taking a long time to compute?

Make sure have a reasonable maximum number of iterations.
Change the error tolerances in scipy.odeint
Reduce the number of bodies in the simulation
Disabling tree saving and relative energy error calculation in Barnes-Hut algorithm can improve performance
Increasing the minimum quad size in the simulation parameters for Barnes-Hut can improve performance

NOTE: The Barnes-Hut algorithm will take longer than should be expected for small N of particles whilist mainting the same accuracy as the brute force algorithm due to
the lower order integration scheme meaning more timesteps of higher temporal resoultion are required, that and the entire script for it is written in
Python whereas the brute force uses scipy and so it is faster at integrating.

How do I change the brute force algorithms error tolerances?

In the StateVector class change the rtol and atol arguments inside scipy.integrate.odeint in the scipy_odeint method

Where are there simulation templates?

Run "generate-templates.py" to generate 8 default simulation templates

How do I quit the Pygame visualisation?

Use the close window button in the top right corner of the Pygame window or
Use the "ESC" key to quit the visualisation

How do I control the camera in the Pygame visualisation?

Use the "W", "A", "S", "D" keys to pan the camera in the y+, x-, y-, x+ directions
Use the mouse scroll wheel to zoom in and out
Press the "E" key to reset the zoom of the camera back to the default zoom.
Press the "R" key to reset the position and zoom of the camera back to default (useful if you get lost)

How do I change what is rendered in the Pygame visualisation?

Press the "X" key to toggle the drawing of the centres of masses
Press the "V" key to toggle the drawing of the velocity and force vectors
Press the "Q" key to toggle the drawing of the Quadtree
Press the "F" key to toggle the drawing of body trails
Press the "C" key to toggle the debug culling hitboxes
