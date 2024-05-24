import numpy as np
import matplotlib.pyplot as plt
import scipy as scipy
import csv

"""
This script simulates the motion of N bodies under the influence of gravity using a brute force approach.

To change the simulation template, change the input_path variable loctaed at the bottom of this script in the INITIALISATION section
to the path of the .csv file containing the simulation parameters and initial conditions of your choosing
"""
#=====================================================================CLASSES AND FUNCTIONS==============================================================================
class Vec2():
    def __init__(self, x, y):
        """
        Construct a 2D vector (x , y)
        
        args:
        x (float): x component of the vector
        y (float): y component of the vector
        """
        self.x = x
        self.y = y
    
    def __str__(self):
        return f"<Vec2 | x={self.x}, y={self.y}>" #define how the vector is printed
        
    def __truediv__(self, scalar): #define division by a scalar
        return Vec2(self.x/scalar , self.y/scalar) 
    
    def __add__(self, vec): #define vector addition
        return Vec2(self.x + vec.x , self.y + vec.y)
    
    def __sub__(self, vec): #define subtraction between vectors
        return Vec2((self.x - vec.x) , (self.y - vec.y))
    
    def magnitude(self): #define the magnitude of a vector
        return np.sqrt(self.x**2 + self.y**2)
    
    def normalise(a): #normalise a vector
        return a / a.magnitude()
    
    def __mul__(self, scalar): #define multipication by a scalar
        return Vec2(self.x * scalar, self.y * scalar)
    
    def distance(a, b): #define the distance between two vectors
        return np.sqrt((a.x - b.x)**2 + (a.y - b.y)**2)

def gravitational_force(body1, body2, G = 6.67E-11):
    """
    Find the gravitational force between two bodies
    
    args:
    body1 (Body)
    body2 (Body)
    G (float): gravitational constant
    
    returns:
    force (Vec2): gravitational force between body1 and body2
    """
    distance = Vec2.distance(body1.position, body2.position)
    direction = (body1.position - body2.position).normalise() 
    return direction * (-G * body1.mass * body2.mass) / (distance**2) #F = r_hat * (-G * m1 * m2) / r^2
    
class Body():
    def __init__(self, position, velocity, mass, radius, id):
        """
        Construct a body
        
        args:
        position (Vec2)
        velocity (Vec2)
        mass (float)
        raius (float)
        id (int)
        """
        self.id = id
        
        self.position = position
        self.velocity = velocity	
        self.mass = mass
        self.radius = radius
        
        self.acceleration = Vec2(0,0) #initialise the acceleration to zero
        self.net_force = Vec2(0,0) #initialise the net force to zero
        
    def __str__(self):
        return f"\nBody = (id: {self.id}, pos: {self.position} m, vel: {self.velocity} m/s, accel: {self.acceleration} m/s^2, net force: {self.net_force} N, mass: {self.mass} kg, radius: {self.radius} m)\n"
        
class StateVector():
    def __init__(self, bodies):
        """
        Construct a state vector
        
        x = [x1 y1 vx1 vy1
             x2 y2 vx2 vy2 
             ... 
             xn yn vxn vyn]
             
        args:
        bodies (list): contains all the bodies in the simulation
        """
        self.bodies = bodies
        self.positions_velocities = np.array([[b.position.x, b.position.y, b.velocity.x, b.velocity.y] for b in bodies]).flatten() #state vector = [x1 y1 vx1 vy1 x2 y2 vx2 vy2 ... xn yn vxn vyn]
        self.solution = None
        self.total_energy = []
        self.dE = []
        
    def __str__(self):
        return f"\nState Vector = \n{self.positions_velocities}\n" 
    #---------------------------------------------------------------------PHYSICS-------------------------------------------------------------------------------------------
    def __calculate_forces(self):
        """
        Calculate the net force on each body from all other bodies
        """
        for body in self.bodies:
            body.net_force = Vec2(0, 0) #reset the net force to zero for each body, otherwise forces will add up each iteration resulting in apocalyptic forces
        
        for body in self.bodies:
            for other_body in self.bodies:
                if body.id == other_body.id:
                    pass
                else:
                    body.net_force += gravitational_force(body, other_body)   
                    
    def get_state_vector_derivative(self, positions_velocities):
        """
        Computes the derivative of the state vector at time t
        
        args:
        positions_velocities (array): state vector
        
        returns:
        state_vector_derivative (array): the derivative of the state vector at time t
        """
        self.__calculate_forces()
        accelerations = [body.net_force / body.mass for body in self.bodies] #doesnt update accelerations of the body objects need to do that elsewhere
        self.state_vector_derivative = np.empty((0, 4))
        
        for i, body in enumerate(self.bodies):
            velocity_acceleration = np.array([positions_velocities[i * 4 + 2], positions_velocities[i * 4 + 3], accelerations[i].x, accelerations[i].y]) #state_vector = [x1 y1 vx1 vy1 x2 y2 vx2 vy2 ... xn yn vxn vyn]
            self.state_vector_derivative = np.vstack((self.state_vector_derivative, velocity_acceleration))
   
        print("Simulating...", end = "\r", flush = True)
        return self.state_vector_derivative
    
    def calculate_energy(self, G = 6.67E-11):
        """
        Calculate the total energy of the system, E = T + V, at each time step
        
        args:
        G (float): gravitational constant
        """
        checked_ids = []
        
        for i, body in enumerate(self.bodies):
            checked_ids.append(body.id)
            kinetic_energy = 0.5 * body.mass * (Vec2(self.solution[:, i * 4 + 2], self.solution[:, i * 4 + 3]).magnitude())**2
            potential_energy = 0
            
            for other_body in self.bodies:
                if body.id == other_body.id and body.id in checked_ids:
                    pass
                else:
                    body_pos = [Vec2(self.solution[:,i * 4 + 0], self.solution[:,i * 4 + 1]) for i, body in enumerate(self.bodies)]
                    other_body_pos = [Vec2(self.solution[:, j * 4 + 0], self.solution[:, j * 4 + 1]) for j, other_body in enumerate(self.bodies)]
                    seperation = [Vec2.distance(body_pos[i], other_body_pos[j]) for i, body in enumerate(self.bodies) for j, other_body in enumerate(self.bodies) if i != j]
                    
                    potential_energy += (-G * body.mass * other_body.mass / seperation[i] )* 0.5 #extra factor of 0.5 since we are double counting the potential energy between two bodies when the arrays are combined    
            
            self.total_energy.append(kinetic_energy + potential_energy)
        
        self.total_energy = np.array(self.total_energy)
        self.total_energy = np.sum(self.total_energy, axis=0)

    def calculate_dE(self):
        """
        Calculate the change in total energy of the system, dE = (E(t) - E(0))/E(0), at each time step
        """
        self.calculate_energy()
        initial_energy = self.total_energy[0]
        print(f"\n\nInitial Energy = {initial_energy} J\n")
        for i in range(len(self.total_energy)):
            curr_energy = self.total_energy[i]
            self.dE.append(np.abs((curr_energy - initial_energy)/initial_energy))  
    #-------------------------------------------------------------------------INTEGRATION------------------------------------------------------------------------------------------       
    def func_wrapper(self, y, t):
        """
        Wrapper function for scipy odeint to solve the state vector derivative
        
        args:
        y (array): state vector
        t (array): NOTE: redundant but odient requires this to be passed as an arg, stupid as f**k.
        
        returns:
        state_vector_derivative (array): the derivative of the state vector at time t
        """
        for i, body in enumerate(self.bodies): #update the positions and velocities of the body objects
            body.position = Vec2(y[i * 4 + 0], y[i * 4 + 1])
            body.velocity = Vec2(y[i * 4 + 2], y[i * 4 + 3])
        
        state_vector_derivative = self.get_state_vector_derivative(y) # y = [x1 y1 vx1 vy1 x2 y2 vx2 vy2 ... xn yn vxn vyn] new state vector from odeint
                                                                                                     
        return state_vector_derivative.flatten()
    
    def scipy_odeint(self, t):
        """
        Integrate the state vector using scipy odeint

        scipy.integrate.odeint(func = callable function i.e. f(y,t), y0 = array, t = array, args=(), ...)
        returns array of (positions and velocities), each row is the positions and velocities of every body at time t
        
        args:
        t (array): time array 0 to max_iter step of dt in seconds
        
        returns:
        solution (array): array of (positions and velocities), each row is the positions and velocities of every body at time t
        """
        func = self.func_wrapper  
        y0 = self.positions_velocities
        self.solution = scipy.integrate.odeint(func, y0, t, rtol = 1e-8, atol = 1e-8) 
        return self.solution           
    #---------------------------------------------------------------PLOTTING-----------------------------------------------------------------------------------------------------
    def draw(self):
        """
        Draw the positions of the bodies in the simulation
        """
        fig = plt.figure(constrained_layout=True, figsize = (16,5), dpi = 200) #create a figure

        axes = fig.subplot_mosaic( #define the subplots of the figure
            """
            AAAAAABBBBB
            AAAAAABBBBB
            AAAAAACCCCC
            AAAAAACCCCC
            """
        )
        
        for i, body in enumerate(self.bodies): #plot the positions of the bodies
            body_x_pos = self.solution[:, i * 4 + 0]
            body_y_pos = self.solution[:, i * 4 + 1]
            color_dict = {
                0: 'red',
                1: 'blue',
                2: 'lime',
                3: 'yellow',
                4: 'purple',
                5: 'cyan',
                6: 'pink',
                7: 'brown',
                8: 'orange',
                9: 'magenta',
                10: 'green',
            }
            
            #plot the positions of the bodies
            axes["A"].plot(body_x_pos, body_y_pos, label=body_names.get(body.id), color=color_dict.get(body.id, 'black'))
            axes["A"].set_xlabel('x (m)', fontsize=14)
            axes["A"].set_ylabel('y (m)', fontsize=14)
            axes["A"].set_title("Test Case Trajectories", fontsize=14)
            axes["A"].legend(loc='upper right', fontsize=8)
            
            #plot the x positions of the bodies over time
            axes["B"].plot(t, body_x_pos, label=body_names.get(body.id), color=color_dict.get(body.id, 'black'))
            axes["B"].axhline(y=0, color='grey', linestyle='dotted', alpha=0.5)
            axes["B"].set_xlabel('t (s)', fontsize=14)
            axes["B"].set_ylabel('x (m)', fontsize=14)
            axes["B"].set_title('X-Positions Over Time', fontsize=14)
            axes["B"].legend(loc='upper right', fontsize=8)
        
        #plot the total energy of the system over time  
        axes["C"].plot(t, self.dE, color='orange')
        axes["C"].set_xlabel('t (s)', fontsize=14)
        axes["C"].set_ylabel('|$\Delta E$| (J)', fontsize=14)
        axes["C"].set_title('Absolute Error in Total Energy', fontsize=14)
                
        plt.show()
#=====================================================================INITIALISATION=======================================================================================     
input_path = 'sim templates/testcase-brute-force.csv' #path to the csv file containing the simulation parameters and initial conditions

with open(input_path , newline='') as csvfile:
    csv_reader = csv.reader(csvfile) 
    bodies = []
    body_names = {}
    #read the simulation parameters from the first 5 rows of the csv file
    max_iter, dt, save_rate = [float(next(csv_reader)[0]) for _ in range(3)]

    for row in csv_reader:
        #create a body object from each row in the csv file
        data = Body(
            Vec2(float(row[2]), float(row[3])),  # position
            Vec2(float(row[4]), float(row[5])),  # velocity
            float(row[6]),  # mass
            float(row[7]),  # radius
            int(row[0]),  # id 
        )
        body_names[int(row[0])] = str(row[1]) 
        
        print(f'{row[1]}\t|  pos: {data.position} m \t | vel: {data.velocity} m/s \n\t| mass: {data.mass:.3g} kg\t | radius: {(data.radius):.3g} m\t\t | mag vel: {Vec2(data.velocity.x, data.velocity.y).magnitude()} m/s\t')
        print('\n')
        bodies.append(data)
        
print(f"name dict = {body_names}\n")
STATE_VECTOR = StateVector(bodies)        
t = np.arange(0, max_iter, dt)
STATE_VECTOR.scipy_odeint(t)
STATE_VECTOR.calculate_dE()
STATE_VECTOR.draw()