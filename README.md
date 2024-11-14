
<p align="center">
  <img src="https://github.com/user-attachments/assets/58e86a45-9eac-4d28-a979-79011e5ebe7a" />
</p>

# **About**
---
A Python based gravitational dynamics simulator, using a quadtree acceleration structure to represent clusters of massive bodies as a hierarchy of their centers of masses, the force calculation is accelerated by applying the Barnes-Hut algorithm. 
Used to study the interactions within and between star systems and galaxies, written as part of a scientific computing module in undergraduate Physics.

<p align="center">
  <img src="https://github.com/user-attachments/assets/ec2741fa-ff7b-45b0-929d-3dcfd5bd4528" width="400" height="400" />
  <img src="https://github.com/user-attachments/assets/cf0330df-a46d-4ce4-8f54-3c5982bcac37" width="400" height="400" />
</p>

# **Features**
---
<table style="width: 100%; table-layout: fixed;">
  <tr>
    <td colspan="2" style="text-align: center; width: 100%"><strong>A gravitational dynamics simulator that uses either:</strong></td>
  </tr>
  <tr>
    <td>
      <ul>
        <li>A Python-based implementation of the Barnes-Hut algorithm with a leapfrog integration scheme.</li>
      </ul>
      <p align="center">
        <img src="https://github.com/user-attachments/assets/21985872-9257-4264-bdf3-4bf309291936" width="400" height="400" />
      </p>
    </td>
    <td style="width: 50%; padding: 10px;">
      <ul>
        <li>A vectorized direct summation scheme using <code>numpy</code> & <code>scipy.integrate.odeint</code>.</li>
      </ul>
      <p align="center">
        <img src="https://github.com/user-attachments/assets/9e2fc945-bff8-4f4d-a877-8e94aaa7378b" width="512" height="200" />
      </p>
    </td>
  </tr>
</table>

<!-- Second Feature Section -->
  <table style="width: 80%; table-layout: fixed;">
    <tr>
      <td colspan="2" style="text-align: center; width: 100%"><strong>5 pre-defined simulation initial conditions:</strong></td>
    </tr>
    <tr>
      <td style="width: 50%; padding: 10px;">
        <ul>
          <li>Binary Star System</li>
        </ul>
        <p align="left">
          <img src="https://github.com/user-attachments/assets/0f59e275-049c-49f4-a9d6-a714ea22fbca" width="300" height="256" />
        </p>
        <ul>
          <li>Burrau's Problem</li>
        </ul>
        <p align="left">
          <img src="https://github.com/user-attachments/assets/f548dac5-0924-4a38-a148-afd5e1926280" width="300" height="256" />
        </p>
      </td>
      <td style="width: 50%; padding: 10px;">
        <ul>
          <li>The Solar System</li>
        </ul>
        <p align="left">
          <img src="https://github.com/user-attachments/assets/f6558622-4002-4235-acd0-c2cfadcbb67d" width="300" height="300" />
        </p>
        <ul>
          <li>Spiral Galaxy</li>
        </ul>
        <p align="left">
          <img src="https://github.com/user-attachments/assets/1121a530-549d-4898-9b16-1a4ae509e8d7" width="300" height="300" />
        </p>
        <ul>
          <li>Galaxy Merger</li>
        </ul>
        <p align="left">
          <img src="https://github.com/user-attachments/assets/cf0330df-a46d-4ce4-8f54-3c5982bcac37" width="300" height="300" />
        </p>
      </td>
    </tr>
  </table>
</div>

<table>
  <tr>
    <td colspan="2"><strong>Interactable simulation results via Pygame, supporting:</strong></td>
  </tr>
  <tr>
    <td>
      <ul>
        <li>Camera movement and zooming</li>
        <li>Frame saving</li>
        <li>GIF creation of simulations</li>
      </ul>
    </td>
    <td>
      <p align="center">
        <img src="https://github.com/user-attachments/assets/3d680856-2f6f-4c44-8a44-cdc9caef6ab2" width="600" height="512" />
      </p>
    </td>
  </tr>
</table>

# **Dependencies**
--- 
```
pip install -r requirements.txt
```
# **Usage**
--- 
<details>
  <summary>Creating a simulation:</summary>
   <ul style="margin-left: 20px;">
      <p>To create a simulation, run <code>generate-templates.py</code> to create 8 default simulation templates. These are .csv files containing the initial conditions for different simulations. The default directory is <code>sim templates/</code>. Or modify the script as desired to produce your own initial conditions.</p>
   </ul>
</details>

<details>
  <summary>Running a simulation:</summary>
  <ul style="margin-left: 20px;">
    <p>To run a simulation, open either <code>nbody-barnes-hut.py</code> or <code>nbody-brute-force.py</code> and change the variable <code>input_path</code> to contain the path of the simulation template you want to run, or use command line arguments <code>--input</code> or <code>-i</code> to change the path variable. By default, both algorithms will use their           test case templates.</p>
    <p><strong>NOTE:</strong> The simulation templates are not compatible between the two different algorithms!</p>
    <ul style="margin-left: 20px;">
</details>

<details>
  <summary>Visualising simulation results:</summary>
  <ul style="margin-left: 20px;">
    <p>Once the simulation has finished processing, to visualise it, find the associated simulation result folder, labeled <code>YYYY-MM-DD HH-MM-SS sim</code> in <code>results/</code>.</p>
    <p>In <code>visualise.py</code>, change the <code>path</code> variable to the path of the simulation result directory and run the script. A Pygame window should open, and you will be able to see the simulation results animate. The default visualisation is set to a precomputed simulation of the test case.</p>
    <p>By default, the visualisation saves every frame to a <code>YYYY-MM-DD HH-MM-SS render</code> folder in <code>rendered/</code>.</p>
    <p>To render a simulation to a .gif file in <code>render-to-gif.py</code>, replace the <code>render_dir</code> variable with the path of the associated <code>YYYY-MM-DD HH-MM-SS render</code> directory, then run the script. Inside the <code>YYYY-MM-DD HH-MM-SS render</code> folder, a gif will be generated from the frames of the visualised simulation with the           name <code>YYYY-MM-DD HH-MM-SS render.gif</code>.</p>
    <p>There are two precomputed simulations: The test case and a simulation of the solar system, ready to be used with the Pygame visualisation.</p>
    <ul style="margin-left: 20px;">
</details>

# **Help**
---

<details>
  <summary>Why is my simulation taking a long time to compute?</summary>
  <ul style="margin-left: 20px;">
    <li>Make sure to set a reasonable maximum number of iterations.</li>
    <li>Change the error tolerances in scipy.odeint for <code>nbody-brute-force.py</code>.</li>
    <li>Reduce the number of bodies in the simulation.</li>
    <li>Disable tree saving and relative energy error calculation in <code>nbody-barnes-hut.py</code> to improve performance.</li>
    <li>Increase the minimum quad size in the simulation parameters for <code>nbody-barnes-hut.py</code> to improve performance.</li>
    <li><strong>Note:</strong> The Barnes-Hut algorithm may take longer for small N of particles, given the lower order integration scheme and the script being fully in Python, whereas <code>nbody-brute-force.py</code> uses SciPy.</li>
  </ul>
</details>

<details>
  <summary>How do I change the <code>nbody-brute-force.py</code> error tolerances?</summary>
  <ul style="margin-left: 20px;">
    <li>In the <code>StateVector</code> class, change the <code>rtol</code> and <code>atol</code> arguments inside the <code>scipy.integrate.odeint</code> method.</li>
  </ul>
</details>

<details>
  <summary>Where are there simulation templates?</summary>
  <ul style="margin-left: 20px;">
    <li>Run <code>generate-templates.py</code> to generate 8 default simulation templates.</li>
  </ul>
</details>

<details>
  <summary>How do I quit the Pygame visualisation?</summary>
  <ul style="margin-left: 20px;">
    <li>Use the close window button in the top right corner of the Pygame window or press the <code>ESC</code> key to quit the visualisation.</li>
  </ul>
</details>

<details>
  <summary>How do I control the camera in the Pygame visualisation?</summary>
  <ul style="margin-left: 20px;">
    <li>Use <code>W</code>, <code>A</code>, <code>S</code>, <code>D</code> keys to pan the camera in the y+, x-, y-, x+ directions.</li>
    <li>Use the mouse <code>scroll wheel</code> to zoom in and out.</li>
    <li>Press <code>E</code> to reset the zoom of the camera back to the default zoom.</li>
    <li>Press <code>R</code> to reset the position and zoom of the camera back to default.</li>
  </ul>
</details>

<details>
  <summary>How do I change what is rendered in the Pygame visualisation?</summary>
  <ul style="margin-left: 20px;">
    <li>Press <code>X</code> to toggle the drawing of the centres of masses.</li>
    <li>Press <code>V</code> to toggle the drawing of velocity and force vectors.</li>
    <li>Press <code>Q</code> to toggle the drawing of the Quadtree.</li>
    <li>Press <code>F</code> to toggle the drawing of body trails.</li>
    <li>Press <code>C</code> to toggle the debug culling hitboxes.</li>
  </ul>
</details>


