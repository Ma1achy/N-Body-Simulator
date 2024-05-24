from astropy.time import Time
from astropy.coordinates import get_body_barycentric_posvel 
import csv
import os
from datetime import datetime, timezone
import plotly.graph_objects as go
import numpy as np
"""
This file generates templates for simulations of the solar system, test case, Burrau's problem for the Barnes Hut alogrithm and the brute force alogrithms.
It also generates templates for a multi-armed spiral galaxy simulation and for the merger of two galaxies for the Barnes Hut algorithm.

The templates should be generated before running the simulations. The templates are .csv files that contain the parameters for the simulations and the bodies to 
be simulated with their inital conditions.

The templates should be located in N-Body-Simulations/sim templates/ directory. The templates are read by the simulation scripts to generate the simulations.
So to reference a template in a simulation script, assign the path to the template to the required variable in the simulation script.

Simulation templates for the Barnes Hut algorithm are named with the suffix 'barnes-hut' and simulation templates for the brute force algorithm are named with the
suffix 'brute-force'. The simulation templates for the associated algorithms are identical except for the simulation parameters. However they are incompatible with 
each other, so a simulation template for the Barnes Hut algorithm cannot be used with the brute force algorithm and vice versa.
"""
#%%
#============================================SOLAR SYSTEM SIMULATION TEMPLATE FOR BARNES HUT===============================================================
"""
Generates a .csv containing the template to construct a simulation of the solar system planets with their positions,
velocities, masses, and radii at a given time.
"""
#simulation parameters
padding = 5*1.496e+11 #(m) padding around the universe to ensure all bodies remain within the universe

axes_size = 70*1.496e+11 + padding #(m) size of the universe in the simulation
min_size = 1E+3 #(m) minimum size of a quad in the simulation
theta = 0.25 #theta parameter for the Barnes-Hut algorithm, smaller values are more accurate but slower
max_iter = 1.49e+6*50 #maximum number of iterations to run the simulation for
dt = 3600/50 #(s) timestep of the simulation
write_step = 3600 #how many nth iterations to write to the csv file (rest are discarded)

#generate data for solar system bodies

time = Time("2023-11-19") #(YYYY-MM-DD) for positions and velocities of the planets at that time

bodies = ['sun', 'mercury', 'venus', 'earth', 'moon', 'mars', 'jupiter', 
          'saturn', 'uranus', 'neptune']

mass_radius_dict = {
    'sun': {'mass': 1.989 * 10**30, 'radius': 696340000},
    'mercury': {'mass': 3.3011 * 10**23, 'radius': 2439700},
    'venus': {'mass': 4.8675 * 10**24, 'radius': 6051800},
    'earth': {'mass': 5.97237 * 10**24, 'radius': 6371000},
    'moon': {'mass': 7.342 * 10**22, 'radius': 1737100},
    'mars': {'mass': 6.4171 * 10**23, 'radius': 3389500},
    'jupiter': {'mass': 1.8982 * 10**27, 'radius': 69911000},
    'saturn': {'mass': 5.6834 * 10**26, 'radius': 58232000},
    'uranus': {'mass': 8.6810 * 10**25, 'radius': 25362000},
    'neptune': {'mass': 1.02413 * 10**26, 'radius': 24622000},
}

fieldnames = ['id', 'name', 'position_x (m)', 'position_y (m)', 'velocity_x (m/s)', 'velocity_y (m/s)', 'mass (kg)', 'radius (AU)']

#create the "sim construction" directory if it doesn't exist
output_directory = 'sim templates'
os.makedirs(output_directory, exist_ok=True)

#generate the file name based on the simulation start time
current_datetime = datetime.now(timezone.utc)
file_name = f'solarsystem-{time.strftime("%Y-%m-%d-%H-%M-%S")}-barnes-hut.csv'
file_path = os.path.join(output_directory, file_name)

with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    csvfile.write(f'{axes_size}\n')
    csvfile.write(f'{min_size}\n')
    csvfile.write(f'{theta}\n')
    csvfile.write(f'{max_iter}\n')
    csvfile.write(f'{dt}\n')
    csvfile.write(f'{write_step}\n')

    #write data for each body
    for idx, body in enumerate(bodies):
        posvel = get_body_barycentric_posvel(body, time) #returns position and velocity of body in AU and AU/day NEED TO CONVERT TO SI UNITS!!!
   
        if body in mass_radius_dict:
            mass = mass_radius_dict[body]['mass']  
            radius = mass_radius_dict[body]['radius']  
        else:
            mass = 1
            radius = 1

        #apply coordinate transformations to centre the solar system to the middle of the 'universe'
        transformed_position_x = axes_size/2 + posvel[0].x.value*1.49598e+11
        transformed_position_y = axes_size/2 + posvel[0].y.value*1.49598e+11
        writer.writerow({  
            'id': idx,
            'name': body,
            'position_x (m)': transformed_position_x,
            'position_y (m)': transformed_position_y,
            'velocity_x (m/s)': posvel[1].x.value*1731460,
            'velocity_y (m/s)': posvel[1].y.value*1731460,
            'mass (kg)': float(mass),
            'radius (AU)': float(radius)
        })
#----------------------------------------PLOTTING THE SOLAR SYSTEM SIMULATION TEMPLATE-------------------------------------------------------
color_map = {
    'Sun': 'yellow',
    'Mercury': 'gray',
    'Venus': 'orange',
    'Earth': 'blue',
    'Moon': 'gray',
    'Mars': 'red',
    'Jupiter': 'brown',
    'Saturn': 'lightyellow',
    'Uranus': 'cyan',
    'Neptune': 'blue'
}

x_positions = []
y_positions = []
velocities_x = []
velocities_y = []
colors = []
radii = []
body_names = []

csv_file_path = 'sim templates/solarsystem-2023-11-19-00-00-00-barnes-hut.csv'

with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)
    axes_size, min_size, theta, max_iter, dt = [float(next(csv_reader)[0]) for _ in range(5)]

with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)

    for row in csv_reader:
        if len(row) >= 4:
            
            body_names.append(row[1].lower())
            x_positions.append(float(row[2]))
            y_positions.append(float(row[3]))
            velocities_x.append(float(row[4]))
            velocities_y.append(float(row[5]))

            body_name = row[1].lower()  
            colors.append(color_map.get(body_name, 'black'))
            radii.append(float(row[7]))
            
fig = go.Figure()

for i in range(len(x_positions)):
    fig.add_trace(go.Scattergl(
        x=[x_positions[i]],
        y=[y_positions[i]],
        mode='markers+text',
        marker=dict(color=colors[i]),
        text=[f'{body_names[i]}'],
        textposition='bottom center',
        textfont=dict(color = colors[i], size=12),
        name='Bodies',
        showlegend=False
    ))

# Plot velocity vectors
    fig.add_trace(go.Scatter(
        x=[x_positions[i], x_positions[i] + velocities_x[i]],
        y=[y_positions[i], y_positions[i] + velocities_y[i]],
        mode='lines+text',  # Change mode to 'lines+text'
        line=dict(color='skyblue', width = 10),
        textposition='top right',
        showlegend=False
    ))
    
#create orbit circles
for i in range(len(x_positions)):
    radius = ((x_positions[i] - (axes_size/2)) ** 2 + (y_positions[i] - (axes_size/2)) ** 2) ** 0.5
    orbit_circle = go.Scatter(
        x=[],
        y=[],
        mode='markers',
        marker=dict(color='white'),
        marker_size=1,
        name='Orbits',
        showlegend=False
    )
    beta = np.linspace(0, 2*np.pi, 100)
    orbit_circle['x'] = ((axes_size/2)) + radius * np.cos(beta)
    orbit_circle['y'] = ((axes_size/2)) + radius * np.sin(beta)
    fig.add_trace(orbit_circle)

#draw x-y axes centred at axes_size/2
fig.add_shape(type="line",
    x0=(axes_size)/2, y0=0, x1=(axes_size)/2, y1=axes_size,
    line=dict(color="darkgray",width=1, dash='dash')
)
fig.add_shape(type="line",
    x0=0, y0=(axes_size)/2, x1=axes_size, y1=(axes_size)/2,
    line=dict(color="darkgray",width=1, dash='dash')
)
#draw a square around the universe
fig.add_shape(type="rect",
    x0=0, y0=0, x1=axes_size, y1=axes_size,
    line=dict(color="red",width=1)
)

fig.update_layout(
    xaxis=dict(range=[0, axes_size],showgrid=False,),
    yaxis=dict(range=[0, axes_size],showgrid=False,),
    paper_bgcolor='black',
    plot_bgcolor='black',
    xaxis_title='X Position (m)',
    yaxis_title='Y Position (m)',
    title='Solar System',
)

fig.show()

#%%
#=============================================TEST CASE SIMULATION TEMPLATE FOR BARNES HUT===============================================================
"""
Generates a .csv containing the template to construct a simulation of the test case:

x_1 = -0.5 AU, y_1 = 0,       v_x1 = 0, V_y1 = -15 km/s
x_2 = 0.5 AU, y_2 = 0,        v_x2 = 0, V_y2 = 15 km/s

1 AU = 1.496e+11 m , M_1 = M_1 = M_sun = 1.989E+30 kg

i/ plot orbits over peroid of 5 years.
ii/ also plot x_1 and x_2 over time.
iii/ also verify the total energy of the system is conserved over time by plotting the relative change in total energy, dE = [E(t=0) - E(t)] / E(t=0)
where E(t) = T + U, the kinetic and potential energy of the system at time t (should find dE < 1E-5).

contains the positions, velocities, masses, and radii of the bodies in the test case.
"""

#simulation parameters
padding = 1.496e+11 #(m) padding around the universe to ensure all bodies remain within the universe

axes_size_test_case = 4*1.496e+11 + padding #(m) size of the universe in the simulation
min_size_test_case = 696340/10 #(m) minimum size of a quad in the simulation
theta_test_case = 0.5 #theta parameter for the Barnes-Hut algorithm, smaller values are more accurate but slower
max_iter_test_case = 43800*50 #maximum number of iterations to run the simulation for
dt_test_case = 30 #(s) timestep of the simulation
test_case_write_step = 24*50 #how many nth iterations to write to the csv file (rest are discarded)

#generate test case bodies:

test_case_bodies = {
    'Star 1': {'id': 0, 'x_position': -0.5*1.496e+11 + (axes_size_test_case)/2, 'y_position': 0 + (axes_size_test_case)/2, 'x_velocity': 0, 'y_velocity': -15000, 'mass': 1.989E+30, 'radius': 696340E+3},
    'Star 2': {'id': 1, 'x_position': 0.5*1.496e+11 + (axes_size_test_case)/2, 'y_position': 0 + (axes_size_test_case)/2, 'x_velocity': 0, 'y_velocity': 15000, 'mass': 1.989E+30, 'radius': 696340E+3}
}

fieldnames_test_case = ['id', 'name', 'position_x (m)', 'position_y (m)', 'velocity_x (m/s)', 'velocity_y (m/s)', 'mass (kg)', 'radius (m)']

#create the "sim construction" directory if it doesn't exist
output_directory = 'sim templates'
os.makedirs(output_directory, exist_ok=True)

#generate the file name based on the simulation start time
current_datetime = datetime.now(timezone.utc)
file_name_test_case = 'testcase-barnes-hut.csv'
file_path = os.path.join(output_directory, file_name_test_case)

with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames_test_case)

    csvfile.write(f'{axes_size_test_case}\n')
    csvfile.write(f'{min_size_test_case}\n')
    csvfile.write(f'{theta_test_case}\n')
    csvfile.write(f'{max_iter_test_case}\n')
    csvfile.write(f'{dt_test_case}\n')
    csvfile.write(f'{test_case_write_step}\n')

    # Write data for each body
    for body_name, body_data in test_case_bodies.items():
        writer.writerow({  
            'id': body_data['id'],
            'name': body_name,
            'position_x (m)': body_data['x_position'],
            'position_y (m)': body_data['y_position'],
            'velocity_x (m/s)': body_data['x_velocity'],
            'velocity_y (m/s)': body_data['y_velocity'],
            'mass (kg)': body_data['mass'],
            'radius (m)': body_data['radius']
        })

#----------------------------------------PLOTTING THE TEST CASE SIMULATION TEMPLATE-------------------------------------------------------
x_positions = []
y_positions = []
x_velocities = []
y_velocities = []
radii = []
body_names = []

csv_file_path = 'sim templates/testcase-barnes-hut.csv'

with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)
    axes_size, min_size, theta, max_iter, dt = [float(next(csv_reader)[0]) for _ in range(5)]
    
with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)

    for row in csv_reader:
        if len(row) >= 4:
            
            body_names.append(row[1].lower())
            x_positions.append(float(row[2]))
            y_positions.append(float(row[3]))
            x_velocities.append(float(row[4]))
            y_velocities.append(float(row[5]))
             
            colors.append(color_map.get('blue'))
            radii.append(float(row[7]))

fig2 = go.Figure()

for i in range(len(x_positions)):
    fig2.add_trace(go.Scattergl(
        x=[x_positions[i]],
        y=[y_positions[i]],
        mode='markers+text',
        marker=dict(color= 'yellow'),
        marker_size=radii[i]/1E+8,	
        showlegend=False
    ))
    
    # Plot velocity vectors
    fig2.add_trace(go.Scatter(
        x=[x_positions[i], x_positions[i] + velocities_x[i]],
        y=[y_positions[i], y_positions[i] + velocities_y[i]],
        mode='lines+text',  
        line=dict(color='skyblue', width = 10),
        textposition='top right',
        showlegend=False
    ))
    
#draw x-y axes centred at (axes_size + padding)/2
fig2.add_shape(type="line",
    x0=(axes_size)/2, y0=0, x1=(axes_size)/2, y1=axes_size,
    line=dict(color="darkgray",width=1, dash='dash')
)
fig2.add_shape(type="line",
    x0=0, y0=(axes_size)/2, x1=axes_size, y1=(axes_size)/2,
    line=dict(color="darkgray",width=1, dash='dash')
)
#draw a square around the universe
fig2.add_shape(type="rect",
    x0=0, y0=0, x1=axes_size, y1=axes_size,
    line=dict(color="red",width=1)
)

fig2.update_layout(
    xaxis=dict(range=[0, axes_size],showgrid=False,),
    yaxis=dict(range=[0, axes_size],showgrid=False,),
    paper_bgcolor='black',
    plot_bgcolor='black',
    xaxis_title='X Position (m)',
    yaxis_title='Y Position (m)',
    title='Test Case',
)

fig2.show()
#%%
#=============================================BURRAU'S PROBLEM SIMULATION TEMPLATE FOR BARNS HUT===============================================================

#simulation parameters
padding = 1.496e+11 #(m) padding around the universe to ensure all bodies remain within the universe

axes_size = 100*1.496e+11 + padding #(m) size of the universe in the simulation
min_size = 696340/10 #(m) minimum size of a quad in the simulation
theta = 0.5 #theta parameter for the Barnes-Hut algorithm, smaller values are more accurate but slower
max_iter = 3E+6 #maximum number of iterations to run the simulation for
dt_test = 63.12 #(s) timestep of the simulation
write_step = 24*50 #how many nth iterations to write to the csv file (rest are discarded)

bodies = {
    'Star 1': {'id': 0, 'x_position': 0 + (axes_size)/2, 'y_position': 0 + (axes_size)/2, 'x_velocity': 0, 'y_velocity': 0, 'mass': 3*1.989E+30, 'radius': 696340E+3},
    'Star 2': {'id': 1, 'x_position': 3*1.496E+11 + (axes_size)/2, 'y_position': 0 + (axes_size)/2, 'x_velocity': 0, 'y_velocity': 0, 'mass': 4*1.989E+30, 'radius': 696340E+3},
    'Star 3': {'id': 2, 'x_position': 3*1.496E+11 + (axes_size)/2, 'y_position': 4*1.496E+11 + (axes_size)/2, 'x_velocity': 0, 'y_velocity': 0, 'mass': 5*1.989E+30, 'radius': 696340E+3}
}

fieldnames = ['id', 'name', 'position_x (m)', 'position_y (m)', 'velocity_x (m/s)', 'velocity_y (m/s)', 'mass (kg)', 'radius (m)']

#create the "sim construction" directory if it doesn't exist
output_directory = 'sim templates'
os.makedirs(output_directory, exist_ok=True)

#generate the file name based on the simulation start time
current_datetime = datetime.now(timezone.utc)
file_name = 'burraus-problem-barnes-hut.csv'
file_path = os.path.join(output_directory, file_name)

with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    csvfile.write(f'{axes_size}\n')
    csvfile.write(f'{min_size}\n')
    csvfile.write(f'{theta}\n')
    csvfile.write(f'{max_iter}\n')
    csvfile.write(f'{dt_test}\n')
    csvfile.write(f'{write_step}\n')

    # Write data for each body
    for body_name, body_data in bodies.items():
        writer.writerow({  
            'id': body_data['id'],
            'name': body_name,
            'position_x (m)': body_data['x_position'],
            'position_y (m)': body_data['y_position'],
            'velocity_x (m/s)': body_data['x_velocity'],
            'velocity_y (m/s)': body_data['y_velocity'],
            'mass (kg)': body_data['mass'],
            'radius (m)': body_data['radius']
        })
#%%
#=============================================TEST CASE SIMULATION TEMPLATE FOR BRUTE FORCE===============================================================
"""
Generates a .csv containing the template to construct a simulation of the test case:

x_1 = -0.5 AU, y_1 = 0,       v_x1 = 0, V_y1 = -15 km/s
x_2 = 0.5 AU, y_2 = 0,        v_x2 = 0, V_y2 = 15 km/s

1 AU = 1.496e+11 m , M_1 = M_1 = M_sun = 1.989E+30 kg

i/ plot orbits over peroid of 5 years.
ii/ also plot x_1 and x_2 over time.
iii/ also verify the total energy of the system is conserved over time by plotting the relative change in total energy, dE = [E(t=0) - E(t)] / E(t=0)
where E(t) = T + U, the kinetic and potential energy of the system at time t (should find dE < 1E-5).

contains the positions, velocities, masses, and radii of the bodies in the test case.
"""

#simulation parameters
max_iter_test_case = 157800000 #maximum number of iterations to run the simulation for
dt_test_case = 10000 #(s) timestep of the simulation
test_case_write_step = 0 #how many nth iterations to write to the csv file (rest are discarded)

#generate test case bodies:

test_case_bodies = {
    'Star 1': {'id': 0, 'x_position': -0.5*1.496e+11, 'y_position': 0, 'x_velocity': 0, 'y_velocity': -15000, 'mass': 1.989E+30, 'radius': 696340E+3},
    'Star 2': {'id': 1, 'x_position': 0.5*1.496e+11, 'y_position': 0, 'x_velocity': 0, 'y_velocity': 15000, 'mass': 1.989E+30, 'radius': 696340E+3}
}

fieldnames_test_case = ['id', 'name', 'position_x (m)', 'position_y (m)', 'velocity_x (m/s)', 'velocity_y (m/s)', 'mass (kg)', 'radius (m)']

#create the "sim construction" directory if it doesn't exist
output_directory = 'sim templates'
os.makedirs(output_directory, exist_ok=True)

#generate the file name based on the simulation start time
current_datetime = datetime.now(timezone.utc)
file_name_test_case = 'testcase-brute-force.csv'
file_path = os.path.join(output_directory, file_name_test_case)

with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames_test_case)
    
    csvfile.write(f'{max_iter_test_case}\n')
    csvfile.write(f'{dt_test_case}\n')
    csvfile.write(f'{test_case_write_step}\n')

    # Write data for each body
    for body_name, body_data in test_case_bodies.items():
        writer.writerow({  
            'id': body_data['id'],
            'name': body_name,
            'position_x (m)': body_data['x_position'],
            'position_y (m)': body_data['y_position'],
            'velocity_x (m/s)': body_data['x_velocity'],
            'velocity_y (m/s)': body_data['y_velocity'],
            'mass (kg)': body_data['mass'],
            'radius (m)': body_data['radius']
        })

#----------------------------------------PLOTTING THE TEST CASE SIMULATION TEMPLATE-------------------------------------------------------
x_positions = []
y_positions = []
x_velocities = []
y_velocities = []
radii = []
body_names = []

csv_file_path = 'sim templates/testcase-brute-force.csv'

with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)
    max_iter, dt = [float(next(csv_reader)[0]) for _ in range(2)]
    
with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)

    for row in csv_reader:
        if len(row) >= 4:
            
            body_names.append(row[1].lower())
            x_positions.append(float(row[2]))
            y_positions.append(float(row[3]))
            x_velocities.append(float(row[4]))
            y_velocities.append(float(row[5]))
             
            colors.append(color_map.get('blue'))
            radii.append(float(row[7]))
            
fig3 = go.Figure()

for i in range(len(x_positions)):
    fig3.add_trace(go.Scattergl(
        x=[x_positions[i]],
        y=[y_positions[i]],
        mode='markers+text',
        marker=dict(color= 'yellow'),
        marker_size=radii[i]/1E+8,	
        showlegend=False
    ))
    
    # Plot velocity vectors
    fig3.add_trace(go.Scatter(
        x=[x_positions[i], x_positions[i] + velocities_x[i]],
        y=[y_positions[i], y_positions[i] + velocities_y[i]],
        mode='lines+text',  
        line=dict(color='skyblue', width = 10),
        textposition='top right',
        showlegend=False
    ))

fig3.update_layout(
    xaxis=dict(range=[-0.75*1.496e+11, 0.75*1.496e+11],showgrid=False,),
    yaxis=dict(range=[-0.75*1.496e+11, 0.75*1.496e+11],showgrid=False,),
    paper_bgcolor='black',
    plot_bgcolor='black',
    xaxis_title='X Position (m)',
    yaxis_title='Y Position (m)',
    title='Test Case Brute Force',
)
fig3.show()
#%%
#=====================================================SOLAR SYSTEM TEMPLATE FOR BRUTE FORCE===============================================================

#simulation parameters
max_iter = 157800000*40 #maximum number of iterations to run the simulation for
dt = 10000 #(s) timestep of the simulation
write_step = 0 #how many nth iterations to write to the csv file (rest are discarded)

#generate data for solar system bodies

time = Time("2023-12-25") #(YYYY-MM-DD) for positions and velocities of the planets at that time

bodies = ['Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter',
            'Saturn', 'Uranus', 'Neptune']

mass_radius_dict = {
    'Sun': {'mass': 1.989 * 10**30, 'radius': 696340000},
    'Mercury': {'mass': 3.3011 * 10**23, 'radius': 2439700},
    'Venus': {'mass': 4.8675 * 10**24, 'radius': 6051800},
    'Earth': {'mass': 5.97237 * 10**24, 'radius': 6371000},
    'Moon': {'mass': 7.342 * 10**22, 'radius': 1737100},
    'Mars': {'mass': 6.4171 * 10**23, 'radius': 3389500},
    'Jupiter': {'mass': 1.8982 * 10**27, 'radius': 69911000},
    'Saturn': {'mass': 5.6834 * 10**26, 'radius': 58232000},
    'Uranus': {'mass': 8.6810 * 10**25, 'radius': 25362000},
    'Neptune': {'mass': 1.02413 * 10**26, 'radius': 24622000},
}

fieldnames = ['id', 'name', 'position_x (m)', 'position_y (m)', 'velocity_x (m/s)', 'velocity_y (m/s)', 'mass (kg)', 'radius (AU)']

#create the "sim construction" directory if it doesn't exist
output_directory = 'sim templates'
os.makedirs(output_directory, exist_ok=True)

#generate the file name based on the simulation start time
current_datetime = datetime.now(timezone.utc)
file_name = f'solarsystem-brute-force-{time.strftime("%Y-%m-%d-%H-%M-%S")}.csv'
file_path = os.path.join(output_directory, file_name)

with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    csvfile.write(f'{max_iter}\n')
    csvfile.write(f'{dt}\n')
    csvfile.write(f'{write_step}\n')
    
    #write data for each body
    
    for idx, body in enumerate(bodies):
        posvel = get_body_barycentric_posvel(body, time)
        
        if body in mass_radius_dict:
            mass = mass_radius_dict[body]['mass']  
            radius = mass_radius_dict[body]['radius']
        else:
            mass = 1
            radius = 1
        
        writer.writerow({
            'id': idx,
            'name': body,
            'position_x (m)': posvel[0].x.value*1.49598e+11,
            'position_y (m)': posvel[0].y.value*1.49598e+11,
            'velocity_x (m/s)': posvel[1].x.value*1731460,
            'velocity_y (m/s)': posvel[1].y.value*1731460,
            'mass (kg)': float(mass),
            'radius (AU)': float(radius)
        })
#%%
#=============================================BURRAU'S PROBLEM SIMULATION TEMPLATE FOR BRUTE FORCE=============================================================== 

#simulation parameters
max_iter = 5E+8 #maximum number of iterations to run the simulation for
dt = 0.6312*100 #(s) timestep of the simulation
write_step = 0 #how many nth iterations to write to the csv file (rest are discarded)

bodies = {
    'Star 1': {'id': 0, 'x_position': 0, 'y_position': 0, 'x_velocity': 0, 'y_velocity': 0, 'mass': 1.989E+30, 'radius': 696340E+3},
    'Star 2': {'id': 1, 'x_position': 3*1.496E+11, 'y_position': 0, 'x_velocity': 0, 'y_velocity': 0, 'mass': 1.989E+30, 'radius': 696340E+3},
    'Star 3': {'id': 2, 'x_position': 3*1.496E+11, 'y_position': 4*1.496E+11, 'x_velocity': 0, 'y_velocity': 0, 'mass': 1.989E+30, 'radius': 696340E+3}
}

fieldnames = ['id', 'name', 'position_x (m)', 'position_y (m)', 'velocity_x (m/s)', 'velocity_y (m/s)', 'mass (kg)', 'radius (m)']

#create the "sim construction" directory if it doesn't exist
output_directory = 'sim templates'
os.makedirs(output_directory, exist_ok=True)

#generate the file name based on the simulation start time
current_datetime = datetime.now(timezone.utc)
file_name = 'barraus-problem-brute-force.csv'
file_path = os.path.join(output_directory, file_name)

with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    csvfile.write(f'{max_iter}\n')
    csvfile.write(f'{dt}\n')
    csvfile.write(f'{write_step}\n')

    # Write data for each body
    for body_name, body_data in bodies.items():
        writer.writerow({  
            'id': body_data['id'],
            'name': body_name,
            'position_x (m)': body_data['x_position'],
            'position_y (m)': body_data['y_position'],
            'velocity_x (m/s)': body_data['x_velocity'],
            'velocity_y (m/s)': body_data['y_velocity'],
            'mass (kg)': body_data['mass'],
            'radius (m)': body_data['radius']
        })
#%%
#========================================================GALAXY SIMULATION TEMPLATE FOR BARNES HUT ALGORITHM===================================================================================   
# Simulation parameters
padding = 10*63241.1 * 1.496e+11  # (m) padding around the universe to ensure all bodies remain within the universe

axes_size = 100e+4*63241.1*1.496e+11 + padding  # (m) size of the universe in the simulation
min_size = 10*63241*1.496e+11  # (m) minimum size of a quad in the simulation
theta = 0.5  # theta parameter for the Barnes-Hut algorithm, smaller values are more accurate but slower
max_iter = 9460*3  # maximum number of iterations to run the simulation for
dt = 3.154e+9*8000*10 # (s) timestep of the simulation
write_step = 1  # how many nth iterations to write to the csv file (rest are discarded)

# Generate data for a multi-armed spiral galaxy using a similar method as the provided code
seed = 1  # seed for random number generator
num_bodies = 500  # number of bodies to generate for the galaxy simulation
a = 0.5
b = 0.6
sx_scale = 0.25
sy_scale = 0.25

galaxy_size = 100e+3*63241.1*1.496e+11
rotation_direction = True # True for anti-clockwise, False for clockwise
galaxy_position_x = axes_size // 2
galaxy_position_y = axes_size // 2

def calculate_positions_with_seed(num_bodies, seed, a, b, galaxy_size, sx_scale, sy_scale, galaxy_position_x, galaxy_position_y, num_branches):
    np.random.seed(seed)
    num_bodies = num_bodies // num_branches  # Divide by the number of branches
    positions = []

    for branch in range(num_branches):
        th = np.random.randn(num_bodies)
        x_spiral = a * np.exp(b * th) * np.cos(th + 2 * np.pi * branch / num_branches) * galaxy_size
        y_spiral = a * np.exp(b * th) * np.sin(th + 2 * np.pi * branch / num_branches) * galaxy_size

        sx = np.random.normal(0, sx_scale * a, num_bodies) * galaxy_size
        sy = np.random.normal(0, sy_scale * a, num_bodies) * galaxy_size

        for i in range(num_bodies):
            x_position = x_spiral[i] + sx[i] + galaxy_position_x
            y_position = y_spiral[i] + sy[i] + galaxy_position_y
            positions.append((x_position, y_position))

    return positions

#calculate star velocities based on their x,y positions and the galaxy's rotation curve, and its rotation direction
#rotation curve is a function of radius, the total mass of the galaxy, and G: so calculate the radius of each star from the center of the 
#galaxy and the total mass of the galaxy when randomly generating each star by saving their masses 
#need to consider both branches of the spiral galaxy as the same set of stars.
#use the radius from the centre of the galaxy to calculate the star's velocity magnitude, 
#then use the angle of the star's position to calculate the x and y components of the velocity along with the rotation direction

def rotation_curve(radius, total_mass):
    return np.sqrt(6.67430e-11 * total_mass / radius)

def calculate_velocities(x_position, y_position, total_mass, galaxy_position_x, galaxy_position_y, rotation_direction):
    radius = np.sqrt((x_position - galaxy_position_x)**2 + (y_position - galaxy_position_y)**2)
    velocity_magnitude = rotation_curve(radius, total_mass)
    
    angle = np.arctan2(y_position - galaxy_position_y, x_position - galaxy_position_x)
    
    if rotation_direction:
        velocity_x = velocity_magnitude * np.cos(angle + np.pi / 2)
        velocity_y = velocity_magnitude * np.sin(angle + np.pi / 2)
    else:
        velocity_x = velocity_magnitude * np.cos(angle - np.pi / 2)
        velocity_y = velocity_magnitude * np.sin(angle - np.pi / 2)
        
    return velocity_x , velocity_y

# Define properties of the supermassive black hole
black_hole_mass = 8.54E+36  # Adjust the mass as needed
black_hole_radius = 1  # Adjust the radius as needed

# Add the supermassive black hole to the simulation data
black_hole = {
    'id': 0,  
    'name': 'supermassive_black_hole',
    'position_x (m)': galaxy_position_x,
    'position_y (m)': galaxy_position_y,
    'velocity_x (m/s)': 0,
    'velocity_y (m/s)': 0,
    'mass (kg)': black_hole_mass,
    'radius (m)': black_hole_radius,
}

total_mass = black_hole_mass  # (kg) total mass of the galaxy

# Generate star data with unique IDs and names
for i in range(num_bodies):
    # Random mass generation
    mass = np.random.uniform(70*1.898E+27, 230*1E+30)  # Adjust the mass range as needed
    total_mass += mass
    
    density = 1408  # (kg/m^3) average density of a star
    radius = (3 * mass / (4 * np.pi * density)) ** (1 / 3)  # (m) radius of the star
    
    # Calculate positions for both branches with random seed
    spiral_galaxy_stars = []  # Initialize the list
    num_branches = 4  # Adjust the number of branches as needed
    positions = calculate_positions_with_seed(num_bodies, seed, a, b, galaxy_size, sx_scale, sy_scale, galaxy_position_x, galaxy_position_y, num_branches)

    # Update star data for both branches
    for j in range(len(positions)):
        id = j
        x_position = positions[j][0]
        y_position = positions[j][1]

        spiral_galaxy_stars.append({
            'id': id + 1,
            'name': f'star {id}',
            'position_x (m)': x_position,
            'position_y (m)': y_position,
            'velocity_x (m/s)': 0,
            'velocity_y (m/s)': 0,
            'mass (kg)': mass,
            'radius (m)': radius,
        })

spiral_galaxy_stars.append(black_hole)

for i in range(num_bodies):
    velocity_x , velocity_y = calculate_velocities(spiral_galaxy_stars[i]['position_x (m)'], spiral_galaxy_stars[i]['position_y (m)'], total_mass, galaxy_position_x, galaxy_position_y, rotation_direction)
    
    spiral_galaxy_stars[i]['velocity_x (m/s)'] = velocity_x
    spiral_galaxy_stars[i]['velocity_y (m/s)'] = velocity_y

# Save the generated data to a CSV file for further inspection
output_directory = 'sim templates'
os.makedirs(output_directory, exist_ok=True)
file_path = os.path.join(output_directory, 'spiral-galaxy-barnes-hut.csv')

fieldname = ['id', 'name', 'position_x (m)', 'position_y (m)', 'velocity_x (m/s)', 'velocity_y (m/s)', 'mass (kg)', 'radius (m)']

with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldname)
    
    csvfile.write(f'{axes_size}\n')
    csvfile.write(f'{min_size}\n')
    csvfile.write(f'{theta}\n')
    csvfile.write(f'{max_iter}\n')
    csvfile.write(f'{dt}\n')
    csvfile.write(f'{write_step}\n')

    # Write data for each body
    for body in spiral_galaxy_stars:
        writer.writerow({  
            'id': body['id'],
            'name': body['name'],
            'position_x (m)': body['position_x (m)'],
            'position_y (m)': body['position_y (m)'],
            'velocity_x (m/s)': body['velocity_x (m/s)'],
            'velocity_y (m/s)': body['velocity_y (m/s)'],
            'mass (kg)': body['mass (kg)'],
            'radius (m)': body['radius (m)']
        })

# Plotting
# Load simulation parameters and galaxy data
csv_file_path = "sim templates/spiral-galaxy-barnes-hut.csv"

body_names = []
x_positions = []
y_positions = []
x_velocities = []
y_velocities = []
radii = []

with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)
    axes_size, min_size, theta, max_iter, dt = [float(next(csv_reader)[0]) for _ in range(5)]

with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)

    for row in csv_reader:
        if len(row) >= 4:
            
            body_names.append(row[1].lower())
            x_positions.append(float(row[2]))
            y_positions.append(float(row[3]))
            x_velocities.append(float(row[4]))
            y_velocities.append(float(row[5]))

            body_name = row[1].lower()  
            radii.append(float(row[7]))

# Create the figure
fig = go.Figure()

# Add markers for stars
for i in range(len(x_positions)):
    fig.add_trace(go.Scattergl(
        x=[x_positions[i]],
        y=[y_positions[i]],
        mode='markers',
        marker=dict(color='yellow', size=2),
        showlegend=False
    ))

# Draw x-y axes centered at (axes_size + padding)/2
fig.add_shape(type="line",
    x0=(axes_size) / 2, y0=0, x1=(axes_size) / 2, y1=axes_size,
    line=dict(color="darkgray", width=1, dash='dash')
)
fig.add_shape(type="line",
    x0=0, y0=(axes_size) / 2, x1=axes_size, y1=(axes_size) / 2,
    line=dict(color="darkgray", width=1, dash='dash')
)

# Draw a square around the universe
fig.add_shape(type="rect",
    x0=0, y0=0, x1=axes_size, y1=axes_size,
    line=dict(color="red", width=1)
)

# Plot velocity vectors
fig.add_trace(go.Scatter(
    x=[x_positions[i], x_positions[i] + x_velocities[i]],
    y=[y_positions[i], y_positions[i] + y_velocities[i]],
    mode='lines+text',  # Change mode to 'lines+text'
    line=dict(color='skyblue', width = 10),
    textposition='top right',
    showlegend=False
))
    
# Layout settings
fig.update_layout(
    xaxis=dict(range=[0, axes_size],showgrid=False,),
    yaxis=dict(range=[0, axes_size],showgrid=False,),
    paper_bgcolor='black',
    plot_bgcolor='black',
    xaxis_title='X Position (m)',
    yaxis_title='Y Position (m)',
    title='Galaxy',
)

fig.show()
#%%     
#=====================================================BINARY GALAXY MERGER SIMULATION TEMPLATE FOR BARNES HUT ALOGRITHM========================================================================            
# Simulation parameters
padding = 10*63241.1 * 1.496e+11  # (m) padding around the universe to ensure all bodies remain within the universe

axes_size = (100e+4*63241.1*1.496e+11 + padding )/2# (m) size of the universe in the simulation
min_size = 10*63241*1.496e+11  # (m) minimum size of a quad in the simulation
theta = 0.5  # theta parameter for the Barnes-Hut algorithm, smaller values are more accurate but slower
max_iter = 9460*3 # maximum number of iterations to run the simulation for
dt = 3.154e+9*8000*10  # (s) timestep of the simulation
write_step = 2 # how many nth iterations to write to the csv file (rest are discarded)

# Generate data for a multi-armed spiral galaxy using a similar method as the provided code
seed = 1  # seed for random number generator
num_bodies = 400  # number of bodies to generate for the galaxy simulation
num_branches = 2
a = 0.5
b = 0.6
sx_scale = 0.25
sy_scale = 0.25

galaxy_size = 100e+3*63241.1*1.496e+11*0.5
net_velocity_x, net_velocity_y = 0, 0
rotation_direction = True # True for anti-clockwise, False for clockwise
galaxy_position_x = axes_size // 2 - 2 * galaxy_size
galaxy_position_y = axes_size // 2 - 2 * galaxy_size

def calculate_positions_with_seed(num_bodies, seed, a, b, galaxy_size, sx_scale, sy_scale, galaxy_position_x, galaxy_position_y, num_branches):
    np.random.seed(seed)
    num_bodies = num_bodies // num_branches  # Divide by the number of branches
    positions = []

    for branch in range(num_branches):
        th = np.random.randn(num_bodies)
        x_spiral = a * np.exp(b * th) * np.cos(th + 2 * np.pi * branch / num_branches) * galaxy_size
        y_spiral = a * np.exp(b * th) * np.sin(th + 2 * np.pi * branch / num_branches) * galaxy_size

        sx = np.random.normal(0, sx_scale * a, num_bodies) * galaxy_size
        sy = np.random.normal(0, sy_scale * a, num_bodies) * galaxy_size

        for i in range(num_bodies):
            x_position = x_spiral[i] + sx[i] + galaxy_position_x
            y_position = y_spiral[i] + sy[i] + galaxy_position_y
            positions.append((x_position, y_position))

    return positions

#data for the second galaxy

num_bodies2 = 600  # number of bodies to generate for the galaxy simulation
num_branches2 = 4
galaxy_size2 = 100e+3*63241.1*1.496e+11*0.5
net_velocity_x2, net_velocity_y2 = 0, 0
rotation_direction2 = False # True for anti-clockwise, False for clockwise
galaxy2_position_x = axes_size // 2 + 2 * galaxy_size
galaxy2_position_y = axes_size // 2 + 2 * galaxy_size

#calculate star velocities based on their x,y positions and the galaxy's rotation curve, and its rotation direction
#rotation curve is a function of radius, the total mass of the galaxy, and G: so calculate the radius of each star from the center of the 
#galaxy and the total mass of the galaxy when randomly generating each star by saving their masses 
#need to consider both branches of the spiral galaxy as the same set of stars.
#use the radius from the centre of the galaxy to calculate the star's velocity magnitude, 
#then use the angle of the star's position to calculate the x and y components of the velocity along with the rotation direction

def calculate_velocities(x_position, y_position, total_mass, galaxy_position_x, galaxy_position_y, rotation_direction):
    radius = np.sqrt((x_position - galaxy_position_x)**2 + (y_position - galaxy_position_y)**2)
    velocity_magnitude = rotation_curve(radius, total_mass)
    
    angle = np.arctan2(y_position - galaxy_position_y, x_position - galaxy_position_x)
    
    if rotation_direction:
        velocity_x = velocity_magnitude * np.cos(angle + np.pi / 2)
        velocity_y = velocity_magnitude * np.sin(angle + np.pi / 2)
    else:
        velocity_x = velocity_magnitude * np.cos(angle - np.pi / 2)
        velocity_y = velocity_magnitude * np.sin(angle - np.pi / 2)
        
    return velocity_x , velocity_y

# Define properties of the supermassive black hole
black_hole_mass = 8.54E+36  # Adjust the mass as needed
black_hole_radius = 1  # Adjust the radius as needed

# Add the supermassive black hole to the simulation data
black_hole = {
    'id': 0,  
    'name': 'supermassive_black_hole',
    'position_x (m)': galaxy_position_x,
    'position_y (m)': galaxy_position_y,
    'velocity_x (m/s)': net_velocity_x,
    'velocity_y (m/s)': net_velocity_y,
    'mass (kg)': black_hole_mass,
    'radius (m)': black_hole_radius,
}

black_hole2 = {
    'id': 999,
    'name': 'supermassive_black_hole2',
    'position_x (m)': galaxy2_position_x,
    'position_y (m)': galaxy2_position_y,
    'velocity_x (m/s)': net_velocity_x2,
    'velocity_y (m/s)': net_velocity_y2,
    'mass (kg)': black_hole_mass,
    'radius (m)': black_hole_radius,
}

total_mass = black_hole_mass  # (kg) total mass of the galaxy
total_mass2 = black_hole_mass

# Generate star data with unique IDs and names
for i in range(num_bodies):
    # Random mass generation
    mass = np.random.uniform(70*1.898E+27, 230*1E+30)  # Adjust the mass range as needed
    #total_mass += mass
    
    density = 1408  # (kg/m^3) average density of a star
    radius = (3 * mass / (4 * np.pi * density)) ** (1 / 3)  # (m) radius of the star
    
    # Calculate positions for both branches with random seed
    spiral_galaxy_stars = []  # Initialize the list
    positions = calculate_positions_with_seed(num_bodies, seed, a, b, galaxy_size, sx_scale, sy_scale, galaxy_position_x, galaxy_position_y, num_branches)
    
    # Update star data for both branches
    for j in range(len(positions)):
        id = j
        x_position = positions[j][0]
        y_position = positions[j][1]

        spiral_galaxy_stars.append({
            'id': id + 1,
            'name': f'star {id}',
            'position_x (m)': x_position,
            'position_y (m)': y_position,
            'velocity_x (m/s)': 0,
            'velocity_y (m/s)': 0,
            'mass (kg)': mass,
            'radius (m)': radius,
        })
        
for i in range(num_bodies2):
    mass = np.random.uniform(70*1.898E+27, 230*1E+30) # Adjust the mass range as needed
    #total_mass2 += mass
    
    density = 1408  # (kg/m^3) average density of a star
    radius = (3 * mass / (4 * np.pi * density)) ** (1 / 3)  # (m) radius of the star
    
    spiral_galaxy_stars2 = []  # Initialize the list
    positions2 = calculate_positions_with_seed(num_bodies2, seed, a, b, galaxy_size2, sx_scale, sy_scale, galaxy2_position_x, galaxy2_position_y, num_branches2)
    
    for j in range(len(positions2)):
        id = j
        x_position = positions2[j][0]
        y_position = positions2[j][1]

        spiral_galaxy_stars2.append({
            'id': id + 1,
            'name': f'star {id}',
            'position_x (m)': x_position,
            'position_y (m)': y_position,
            'velocity_x (m/s)': 0,
            'velocity_y (m/s)': 0,
            'mass (kg)': mass,
            'radius (m)': radius,
        })
    
spiral_galaxy_stars.append(black_hole)
spiral_galaxy_stars.append(black_hole2)

for i in range(num_bodies):
    velocity_x , velocity_y = calculate_velocities(spiral_galaxy_stars[i]['position_x (m)'], spiral_galaxy_stars[i]['position_y (m)'], total_mass, galaxy_position_x, galaxy_position_y, rotation_direction)
    
    spiral_galaxy_stars[i]['velocity_x (m/s)'] = velocity_x + net_velocity_x
    spiral_galaxy_stars[i]['velocity_y (m/s)'] = velocity_y + net_velocity_y

for i in range(num_bodies2):
    velocity_x , velocity_y = calculate_velocities(spiral_galaxy_stars2[i]['position_x (m)'], spiral_galaxy_stars2[i]['position_y (m)'], total_mass2, galaxy2_position_x, galaxy2_position_y, rotation_direction2)
    
    spiral_galaxy_stars2[i]['velocity_x (m/s)'] = velocity_x + net_velocity_x2
    spiral_galaxy_stars2[i]['velocity_y (m/s)'] = velocity_y + net_velocity_y2
    
spiral_galaxy_stars.extend(spiral_galaxy_stars2)

# Save the generated data to a CSV file for further inspection
output_directory = 'sim templates'
os.makedirs(output_directory, exist_ok=True)
file_path = os.path.join(output_directory, 'galaxy-merge-barnes-hut.csv')

fieldname = ['id', 'name', 'position_x (m)', 'position_y (m)', 'velocity_x (m/s)', 'velocity_y (m/s)', 'mass (kg)', 'radius (m)']

with open(file_path, 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldname)
    
    csvfile.write(f'{axes_size}\n')
    csvfile.write(f'{min_size}\n')
    csvfile.write(f'{theta}\n')
    csvfile.write(f'{max_iter}\n')
    csvfile.write(f'{dt}\n')
    csvfile.write(f'{write_step}\n')

    # Write data for each body
    for body in spiral_galaxy_stars:
        writer.writerow({  
            'id': body['id'],
            'name': body['name'],
            'position_x (m)': body['position_x (m)'],
            'position_y (m)': body['position_y (m)'],
            'velocity_x (m/s)': body['velocity_x (m/s)'],
            'velocity_y (m/s)': body['velocity_y (m/s)'],
            'mass (kg)': body['mass (kg)'],
            'radius (m)': body['radius (m)']
        })

# Plotting
# Load simulation parameters and galaxy data
csv_file_path = "sim templates/galaxy-merge-barnes-hut.csv"

body_names = []
x_positions = []
y_positions = []
x_velocities = []
y_velocities = []
radii = []

with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)
    axes_size, min_size, theta, max_iter, dt = [float(next(csv_reader)[0]) for _ in range(5)]

with open(csv_file_path, newline='') as csvfile:
    csv_reader = csv.reader(csvfile)

    for row in csv_reader:
        if len(row) >= 4:
            
            body_names.append(row[1].lower())
            x_positions.append(float(row[2]))
            y_positions.append(float(row[3]))
            x_velocities.append(float(row[4]))
            y_velocities.append(float(row[5]))

            body_name = row[1].lower()  
            radii.append(float(row[7]))

# Create the figure
fig = go.Figure()

# Add markers for stars
for i in range(len(x_positions)):
    fig.add_trace(go.Scattergl(
        x=[x_positions[i]],
        y=[y_positions[i]],
        mode='markers',
        marker=dict(color='yellow', size=2),
        showlegend=False
    ))

# Draw x-y axes centered at (axes_size + padding)/2
fig.add_shape(type="line",
    x0=(axes_size) / 2, y0=0, x1=(axes_size) / 2, y1=axes_size,
    line=dict(color="darkgray", width=1, dash='dash')
)
fig.add_shape(type="line",
    x0=0, y0=(axes_size) / 2, x1=axes_size, y1=(axes_size) / 2,
    line=dict(color="darkgray", width=1, dash='dash')
)

# Draw a square around the universe
fig.add_shape(type="rect",
    x0=0, y0=0, x1=axes_size, y1=axes_size,
    line=dict(color="red", width=1)
)

# Plot velocity vectors
fig.add_trace(go.Scatter(
    x=[x_positions[i], x_positions[i] + x_velocities[i]],
    y=[y_positions[i], y_positions[i] + y_velocities[i]],
    mode='lines+text',  # Change mode to 'lines+text'
    line=dict(color='skyblue', width = 10),
    textposition='top right',
    showlegend=False
))
    
# Layout settings
fig.update_layout(
    xaxis=dict(range=[0, axes_size],showgrid=False,),
    yaxis=dict(range=[0, axes_size],showgrid=False,),
    paper_bgcolor='black',
    plot_bgcolor='black',
    xaxis_title='X Position (m)',
    yaxis_title='Y Position (m)',
    title='Galaxy',
)

fig.show()