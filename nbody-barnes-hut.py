import os
import numpy as np
import plotly.graph_objects as go
import scipy as scipy
import time
import csv
import argparse
from datetime import datetime
"""
This script will simulate the motion of N bodies under gravity in 2D space using the Barnes-Hut algorithm.

To change the simulation template, change the csv file in the input directory. By either changing the input path variable in the main function, 
or by passing the path to the csv file as a command line argument when running the script via -i or --input.
"""
#=======================================================BACKGROUND ON ALGORITHM=====================================================================
"""
barnes-hut:
    
internal node: doesn't directly contain a body only children that contain bodies

external node: only contains one body

Recursive algorithm:
    
To build the tree (starting at the top, the root):
    1/ (base case) if the current node x is empty: put the new body b here.
    2/ If the current node x is an internal node : update the centre of mass and mass of x, 
    then place the body b in the correct quadrant
    3/ If the current node x is an external node: subdivide into four smaller quadrants. 
    Put the two bodies in their correct quadrants (there may be several subdivisions during 
    a single insertion), then update the cente of masses and mass of x afterwards (IMPORTANT TO DO AFTERWARDS!)

i.e.
 __________________
|       |  B |    |                        0(1,2,3,4)                          depth: 0
|   A   |____|____|              _______________|_________________________
|       |    | C  |             |         |             |                | 
|_______|____|____|            1(A)  2(i,ii,iii,iv)    3(i,ii,iii,iv)  4(H)    depth: 1
|   |   |         |                      |                  |     
|___|_D_|         |                {i(B) , iv(C)}   {ii(D) , iii(E), iv(F)}    depth: 2
| E |   |     H   |                 
|___|_F_|_________|

each node is going to contain its, mass, centre of mass, it's side length (s), and its absoulte
positon in space (top left corner), and its 4 children (for an internal node)

computing the force:

once we have our tree, computing the force for each body is approximated by doing the following on each
of the bodies (starting at the top of the tree doing a breadth first traversal):
    
    1/ (base case) If the current node is an external node (and not the body b) calculate the force exterted by the
    current node on b directly, and add this amount to b's net force.
    2/ otherwise compute s/d <= theta, treat this internal node as a single body, and calculate the force
    it exerts on the b, and add this the to b's net force
    3/ otherwise, run the prodedure recursively on each of the current nodes children
   
internal nodes of the same level are traversed breadth wise, the rest are traversed depth wise.
for the last node it can't exert a force on itself

condition for traversing internal nodes breadth wise:
    compute ratio s/d = theta, where s is the side length of the internal node, and d is the distance between
    the body and the nodes centre of mass. 
    theta ~= 0.5 - 0.7
    larger theta = less accurate, but faster . smaller theta = more accurate , slower
    theta = 0 -> direct sum O(N^2) (BAD!)

complexity: building tree takes O(nlogn) time , computing the force generally takes O(nlogn) (although this is
            dependent on the value of theta)

Generally more efficent to just delete the tree and rebuild it every iteration, rather than updating it, only starts to become more efficient
to update around n = 1M bodies and python is already too slow to simulate that many bodies anyway.
"""
#==============================================================CLASSES AND FUNCTIONS===========================================================

class Vec2():
    def __init__(self, x, y):
        """
        Construct a 2D vector (x , y)
        
        args:
        x (float): the x component of the vector
        y (float): the y component of the vector
        """
        self.x = x
        self.y = y
    
    def __str__(self):
        return f"<Vec2 | x={self.x} y={self.y}>" #define how the vector is printed
        
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
    
#------------------------------------------------------------------------------------------------------------------------------------------------

def gravitational_force(mass_1 , mass_2, position_1 , position_2, G = 6.67E-11): #define the gravitational force between two bodies
    """
    Calculates the gravitational force between two bodies
    
    args:
    mass_1 (float): the mass of the first body
    mass_2 (float): the mass of the second body
    position_1 (Vec2): the position of the first body
    position_2 (Vec2): the position of the second body
    G (float): the gravitational constant
    
    returns:
    (Vec2): the gravitational force between the two bodies
    """
    #F = G * m1m2/r^2 * r_hat , r_hat = r/|r|
    r = position_1 - position_2
    return Vec2.normalise(r) * - G * mass_1 * mass_2 / r.magnitude()**2

#------------------------------------------------------------------------------------------------------------------------------------------------  
class Body():
    def __init__(self, position, velocity, mass, radius, id):
        """
        Construct a body
        
        args:
        positon (Vec2): the x , y coordinate of the body
        velocity (Vec2): the x and y components of the bodies velocity
        mass (float): the mass of the object 
        radius (float): the radius of the object
        id (int): the unique id of the body
        """
        
        #body attributes
        self.id = id 
        
        #body physics
        self.position = position
        self.velocity = velocity
        self.mass = mass
        self.radius = radius
        self.net_force = Vec2(0,0) #intaialise the net force on the body to zero
        self.acceleration = Vec2(0,0) #initialise the acceleration of the body to zero
        self.energy = 0 #initialise the energy of the body to zero

#------------------------------------------------------------------------------------------------------------------------------------------------

class Quad():
    def __init__(self, top_left, bot_right, data = None):
        """
        Construct a Quad
        
        args:
        top_left (Vec2)
        bot_right (Vec2)
        data (Body)
        """
        
        self.top_left = top_left
        self.bot_right = bot_right
        
        self.data = data

        self.mass_total = None
        self.com = None
        
        self.tl_quad = None
        self.tr_quad = None
        self.bl_quad = None
        self.br_quad = None
        
    def is_leaf(self):
        """
        Checks if the quad is a leaf
        
        returns:
        (bool): True if the quad is a leaf, False otherwise
        """
        return not(self.tl_quad or self.tr_quad or self.bl_quad or self.br_quad)
    
    def side_length(self):
        return np.abs(self.top_left.x - self.bot_right.x)

#------------------------------------------------------------------------------------------------------------------------------------------------

class QuadTree():
    """
    Construct a Quad Tree of Quads
    """
    def __init__(self, top_left, bot_right, min_size):
        
        #quad tree attributes
        self.root = Quad(top_left, bot_right)
        self.min_size = min_size #minium size of quad
        self.size = 0
        self.seen_ids = []
        
        #physics attributes
        self.total_energy = 0
        self.max_accel = 0
        self.max_vel = 0
        
    def insert(self, data):
        """
        Inserts a quad into the quad tree
        
        args:
        data (Body): the body to insert into the quad tree
        """
        self.__quad_insert(data, self.root)
        
    def __quad_insert(self, data, quad):
        """
        Find the correct quad to insert the data into
        
        args:
        data (Body): the body to insert into the quad tree
        quad (Quad): the current quad that we are inserting the body into
        """
        if quad is None :
            return
        
        if not self.__in_bounds(data.position, quad):
            return
        
        if np.abs(quad.top_left.x - quad.bot_right.x) <= self.min_size:
            if quad.data is None:
                quad.data = data
                self.size += 1
            return
        
        if data.position.y <= (quad.top_left.y + quad.bot_right.y) / 2 : #check if in top or bottom half
        
            #top left quadrant
            if data.position.x <= (quad.top_left.x + quad.bot_right.x) / 2 : #check if in left or right half
                if quad.tl_quad is None:
                    quad.tl_quad = Quad(
                        quad.top_left,
                        (quad.top_left + quad.bot_right) / 2 
                        )
                self.__quad_insert(data, quad.tl_quad)
                
            #top right quadrant
            else:
                if quad.tr_quad is None:
                    quad.tr_quad = Quad(
                        Vec2((quad.top_left.x + quad.bot_right.x) / 2,
                             quad.top_left.y),
                        Vec2(quad.bot_right.x,
                            (quad.top_left.y + quad.bot_right.y) / 2)
                        )
                self.__quad_insert(data, quad.tr_quad)
        else:
            #bot left quadrant
            if data.position.x <= (quad.top_left.x + quad.bot_right.x) / 2 : #check if in left or right half
                if quad.bl_quad is None:
                    quad.bl_quad = Quad(
                        Vec2(quad.top_left.x,
                            (quad.top_left.y + quad.bot_right.y) / 2),
                        Vec2((quad.top_left.x + quad.bot_right.x) / 2,
                             quad.bot_right.y)
                        )
                self.__quad_insert(data, quad.bl_quad)
                
            #bot right quadrant
            else:
                if quad.br_quad is None:
                    quad.br_quad = Quad(
                        Vec2((quad.top_left.x + quad.bot_right.x)/2,
                             (quad.top_left.y + quad.bot_right.y)/2),
                            quad.bot_right
                        )
                self.__quad_insert(data, quad.br_quad)
            
    def __in_bounds(self, position, quad):
        """
        Finds if a position is within a quad
        
        args:
        position (Vec2): the position to check
        quad (Quad): the quad to check if the position is within
        
        returns:
        (bool): True if the position is within the quad, False otherwise
        """
        return (position.x >= quad.top_left.x and position.x <= quad.bot_right.x 
                and position.y >= quad.top_left.y and position.y <= quad.bot_right.y)

    #---------------------------------------------------PHYSICS------------------------------------------------------
    
    def calculate_com_approx(self, frame, save_rate, quaddatafilename):
        """
        Finds the centre of mass of each quad, and the total mass of each quad
        
        args:
        frame (int): the current frame
        save_rate (int): the rate at which to save the quad data to the csv file
        quaddatafilename (str): the path to the quad data csv file
        """
        self.__com_bfs(self.root, 0, frame, save_rate, quaddatafilename)
        
    def __com_bfs(self, quad, d, frame, save_rate, quaddatafilename):
        """
        Conducts a breadth first traversal of the quad tree, and calculates the total_mass and centre of mass of each quad
        
        args:
        quad (Quad): the current quad that we are calculating the centre of mass on
        d (int): the depth of the current quad
        frame (int): the current frame
        save_rate (int): the rate at which to save the quad data to the csv file
        quaddatafilename (str): the path to the quad data csv file
        """
        save_quads = True #set to true to save the quads to a csv file
        
        if quad is None:
            return 0, Vec2(0,0)
        
        tl_mass, tl_posvec = self.__com_bfs(quad.tl_quad, d+1, frame, save_rate, quaddatafilename) 
        tr_mass, tr_posvec = self.__com_bfs(quad.tr_quad, d+1, frame, save_rate, quaddatafilename) 
        bl_mass, bl_posvec = self.__com_bfs(quad.bl_quad, d+1, frame, save_rate, quaddatafilename) 
        br_mass, br_posvec = self.__com_bfs(quad.br_quad, d+1, frame, save_rate, quaddatafilename) 
        
        mass_total = tl_mass + tr_mass + bl_mass + br_mass 
        if mass_total == 0: #NOTE: if the system is give zero-masses then this check will break (along with the rest of the physics so whatever).
            mass_total = quad.data.mass
            com = quad.data.position 
        else:
            com = (tl_posvec*tl_mass + tr_posvec*tr_mass + bl_posvec*bl_mass + br_posvec*br_mass) / mass_total #calculate the centre of mass of the quad
    
        quad.mass_total = mass_total
        quad.com = com
        
        if save_quads:
            #save quad to csv file
            if frame % save_rate == 0:
                with open(quaddatafilename, 'a', newline='') as file:
                    writer = csv.writer(file)
                    quad_data = [frame, quad.top_left.x, quad.top_left.y, quad.bot_right.x, quad.bot_right.y, quad.mass_total, quad.com.x, quad.com.y]
                    writer.writerow(quad_data)
        
        return quad.mass_total, quad.com
    
    def calculate_forces(self, theta):
        """
        Finds the net force on each body in the system
        
        args:
        theta (float): the ratio s/d, where s is the side length of the quad, and d is the distance between the body and the quad's centre of mass
        """
        self.__force_per_body_bfs(self.root, theta)
    
    def __force_per_body_bfs(self, quad, theta):
        """
        Finds the net force on each body in the system by doing a breadth first traversal of the quad tree
        
        args:
        quad (Quad): the current quad that we are calculating the net force on
        theta (float): the ratio s/d, where s is the side length of the quad, and d is the distance between the body and the quad's centre of mass
        """	
        if quad is None:
            return
        
        if quad.is_leaf() and quad.data:
            
            quad.data.net_force = Vec2(0,0)
            self.__calculate_net_force_bfs(quad, self.root, theta)
        
        self.__force_per_body_bfs(quad.tl_quad, theta)
        self.__force_per_body_bfs(quad.tr_quad, theta)
        self.__force_per_body_bfs(quad.bl_quad, theta)
        self.__force_per_body_bfs(quad.br_quad, theta)
        
    def __calculate_net_force_bfs(self, quad_ref, quad_curr, theta):
        """
        for every body, we have to calculate the net force on that body due to all other bodies in the system.
        this is done by traversing the quad tree, and deciding whever to approximate the force exerted by a quad
        or to recurse down to its children and calculate the force exerted by each child:
        
        ==================== for any leaf node =====================

        We just have to calculate the pair-wise force directly and add it to a net force vector of the current body
        
        ==================== for any brach node ====================

        Calculate the ratio "s/d" by dividing the quad's side length by the distance between the current body and the quad's center of mass, defining theta.
        Recurse down the tree, comparing "s/d" to theta at each level. 

        If "s/d" is less than or equal to theta, approximate the forces using the center of mass and total mass of the branch node.

        If "s/d" is greater than theta:
        Treat the quad and its children as individual pair-wise interactions.
        Recurse down to leaf nodes, calculate the forces exerted by each leaf node, and add them to the net force vector of the current body.

        args:
        quad_ref (Quad): the current quad that we are calculating the net force on
        quad_curr (Quad): the quad that we are currently recursing down finding the force it applies to quad_ref
        theta (float): the ratio s/d, where s is the side length of the quad, and d is the distance between the body and the quad's centre of mass
        """   
        if quad_curr is None:
            return
        
        if quad_curr.data is quad_ref.data:
            return
        
        if quad_curr.is_leaf() and quad_curr.data:
            quad_ref.data.net_force += gravitational_force(quad_ref.data.mass, quad_curr.data.mass, quad_ref.data.position, quad_curr.data.position) #pair-wise force between two leaf quads
            
        s = quad_curr.side_length()
        d = Vec2.distance(quad_ref.data.position, quad_curr.com)
        
        if d == 0: #d can be zero because the parent quad can share the same cenre of mass as the bodies position
            return   
        
        s_d = s/d
        
        if s_d <= theta: #use force approx from com and masss_total
            quad_ref.data.net_force += gravitational_force(quad_ref.data.mass, quad_curr.mass_total, quad_ref.data.position, quad_curr.com)
            return 
        
        else:
            self.__calculate_net_force_bfs(quad_ref, quad_curr.tl_quad, theta)
            self.__calculate_net_force_bfs(quad_ref, quad_curr.tr_quad, theta)
            self.__calculate_net_force_bfs(quad_ref, quad_curr.bl_quad, theta)
            self.__calculate_net_force_bfs(quad_ref, quad_curr.br_quad, theta)
        
    def integrate(self, dt, theta, frame, save_rate, new_tree, bodydatafilename="sim.csv"):
        """
        Finds the velocity and new position of each body in the system
        """
        self.__leapfrog_intergation_bfs(self.root, dt, theta, frame, save_rate, new_tree, bodydatafilename)
        
    def __leapfrog_intergation_bfs(self, quad, dt, theta, frame, save_rate, new_tree, bodydatafilename):
        """
        Calculates the velocity and new position of each body in the system by doing a breadth first traversal of the quad tree and using the leapfrog integration method
        (this is a symplectic integrator, which means it conserves energy).
        
        args:
        quad (Quad): the current quad that we are calculating the velocity and new position on
        dt (float): the time step
        theta (float): the ratio s/d, where s is the side length of the quad, and d is the distance between the body and the quad's centre of mass
        frame (int): the current frame
        save_rate (int): the rate at which to save the body data to the csv file
        new_tree (QuadTree): the new tree to insert the bodies into
        bodydatafilename (str): the path to the body data csv file
        """
        
        if quad is None:
            return
        
        if quad.is_leaf() and quad.data: #if the quad is a leaf then integrate the body using leapfrog integration
            
            quad.data.acceleration = quad.data.net_force / quad.data.mass
            quad.data.position += quad.data.velocity * dt + quad.data.acceleration * dt**2 / 2
            quad.data.velocity += quad.data.acceleration * dt / 2
            quad.acceleration = quad.data.velocity / dt
            quad.data.velocity += quad.data.acceleration * dt / 2
            
            new_tree.insert(Body(quad.data.position, quad.data.velocity, quad.data.mass, quad.data.radius, quad.data.id)) #insert the body into the new tree with its new position and velocity
            
            if frame % save_rate == 0: #save the body to csv file
                with open(bodydatafilename, 'a', newline='') as file:
                    writer = csv.writer(file)
                    body_data = [frame, quad.data.id, quad.data.position.x, quad.data.position.y, quad.data.velocity.x, quad.data.velocity.y, quad.data.net_force.x, quad.data.net_force.y, quad.data.mass, quad.data.radius]
                    writer.writerow(body_data)                     
        else:
            self.__leapfrog_intergation_bfs(quad.tl_quad, dt, theta, frame, save_rate, new_tree, bodydatafilename)
            self.__leapfrog_intergation_bfs(quad.tr_quad, dt, theta, frame, save_rate, new_tree, bodydatafilename)
            self.__leapfrog_intergation_bfs(quad.bl_quad, dt, theta, frame, save_rate, new_tree, bodydatafilename)
            self.__leapfrog_intergation_bfs(quad.br_quad, dt, theta, frame, save_rate, new_tree, bodydatafilename)
    
    def get_net_energy(self):
        """
        Returns the total energy of the system
        """
        return self.total_energy
    
    def calculate_energy(self, G = 6.67E-11):
        """
        Finds the total energy of the system, E = T + V
        
        args:
        G (float): the gravitational constant
        """
        self.__kinetic_energy_bfs(self.root, G)   
          
    def __kinetic_energy_bfs(self, quad, G):
        """
        Calculates the kinetic energy and then the total energy of the system by doing a breadth first traversal of the quad tree
        
        args:
        quad (Quad): the current quad that we are calculating the kinetic energy on
        G (float): the gravitational constant
        """
        if quad is None:
            return 
         
        if quad.is_leaf() and quad.data:
            quad.data.energy += 0.5 * quad.data.mass * quad.data.velocity.magnitude()**2 #kinetic energy of a leaf quad
            self.__potential_energy_bfs(quad, self.root, G) 
            self.total_energy += quad.data.energy 
           
        self.__kinetic_energy_bfs(quad.tl_quad, G)
        self.__kinetic_energy_bfs(quad.tr_quad, G)
        self.__kinetic_energy_bfs(quad.bl_quad, G)
        self.__kinetic_energy_bfs(quad.br_quad, G)
    
    def __potential_energy_bfs(self, quad_ref, quad_curr, G):
        """
        Calculates the potential energy of the system by doing a breadth first traversal of the quad tree
        
        args:
        quad_ref (Quad): the current quad that we are calculating the potential energy on
        quad_curr (Quad): the quad that we are currently recursing down finding the potential energy it applies to quad_ref
        G (float): the gravitational constant
        """
        if quad_curr is None:
            return
                
        if quad_curr.is_leaf() and quad_curr.data and quad_curr.data is not quad_ref.data:
            if quad_curr.data.id in [ids for _, ids in enumerate(self.seen_ids)]: #if the body has already had its potential energy calculated then return
                return
            else:
                self.seen_ids.append(quad_ref.data.id)
                pe = -G * quad_ref.data.mass * quad_curr.data.mass / Vec2.distance(quad_ref.data.position, quad_curr.data.position) #gravitational potential energy between two leaf quads
                quad_ref.data.energy += pe
            
        else:
            self.__potential_energy_bfs(quad_ref, quad_curr.tl_quad, G)
            self.__potential_energy_bfs(quad_ref, quad_curr.tr_quad, G)
            self.__potential_energy_bfs(quad_ref, quad_curr.bl_quad, G)
            self.__potential_energy_bfs(quad_ref, quad_curr.br_quad, G)
            
    #---------------------------------------------------PLOTTING------------------------------------------------------  
    def draw(self, fig, draw_bodies, draw_quads, draw_com, draw_forces, draw_velocity):
        """
        Draws the plot of the quad tree for the inital conditions of the simulation.
        
        args:
        fig (go.Figure): the figure to add the traces to
        draw_bodies (bool): whether or not to draw the bodies
        draw_quads (bool): whether or not to draw the quads
        draw_com (bool): whether or not to draw the centre of mass of each quad
        draw_forces (bool): whether or not to draw the net force on each body
        draw_velocity (bool): whether or not to draw the velocity of each body
        """
        traces = []
        self.__draw_bfs(self.root, 0, traces, draw_bodies, draw_quads, draw_com, draw_forces, draw_velocity)
        fig.add_traces(traces)  
        
    def __draw_bfs(self, quad, d, traces, draw_bodies ,draw_quads, draw_com, draw_forces, draw_velocity):
        """
        Draws the quads and bodies in the quad tree along with the centre of mass of each quad for the inital conditions of 
        the simulation by conducting a breadth first traversal of the quad tree.
        
        args:
        quad (Quad): the current quad that we are drawing
        d (int): the depth of the current quad
        traces (list): the list of traces that we are appending to
        draw_bodies (bool): whether or not to draw the bodies
        draw_quads (bool): whether or not to draw the quads
        draw_com (bool): whether or not to draw the centre of mass of each quad
        draw_forces (bool): whether or not to draw the net force on each body
        draw_velocity (bool): whether or not to draw the velocity of each body
        """
        if quad is None:
            return
          
        if draw_bodies and quad.data: #draw the bodies
            traces.append(go.Scattergl(
                x=[quad.data.position.x],
                y=[quad.data.position.y],
                mode='markers',
                marker=dict(color='yellow', size = 5),
                showlegend=False
                
            ))
        
        if draw_com: #draw the centre of mass of each quad
            if quad.com:
                
                traces.append(go.Scatter(
                    x=[quad.com.x],
                    y=[quad.com.y],
                    mode='markers',
                    marker=dict(color='red', size=1E+10),
                    showlegend=False
                    
                    
                )) 
    
        if draw_quads: #draw the quads
            traces.append(go.Scattergl(
                x = [quad.top_left.x, quad.bot_right.x, quad.bot_right.x, quad.top_left.x, quad.top_left.x],
                y = [quad.top_left.y, quad.top_left.y, quad.bot_right.y, quad.bot_right.y, quad.top_left.y],
                
                mode='lines',
                line=dict(
                    color='lime'
                    ),
                showlegend=False
            ))
        
        if draw_forces and quad.data is not None: #draw the net force on each body
           
            traces.append(go.Scattergl(
                x= [quad.data.position.x, quad.data.position.x + Vec2.normalise(quad.data.net_force).x],
                y= [quad.data.position.y, quad.data.position.y + Vec2.normalise(quad.data.net_force).y],
                mode='lines',
                marker=dict(color='yellow', size=10),
                showlegend=False
                ))
        
        if draw_velocity and quad.data is not None: #draw the velocity of each body
 
            traces.append(go.Scattergl(
                x= [quad.data.position.x, quad.data.position.x + Vec2.normalise(quad.data.velocity).x],
                y= [quad.data.position.y, quad.data.position.y + Vec2.normalise(quad.data.velocity).y],
                mode='lines',
                marker=dict(color='cyan', size=10),
                showlegend=False
                ))
    
        self.__draw_bfs(quad.tl_quad, d+1, traces, draw_bodies, draw_quads, draw_com, draw_forces, draw_velocity)
        self.__draw_bfs(quad.tr_quad, d+1, traces, draw_bodies, draw_quads, draw_com, draw_forces, draw_velocity)
        self.__draw_bfs(quad.bl_quad, d+1, traces, draw_bodies, draw_quads, draw_com, draw_forces, draw_velocity)
        self.__draw_bfs(quad.br_quad, d+1, traces, draw_bodies, draw_quads, draw_com, draw_forces, draw_velocity)
#------------------------------------------------------------------------------------------------------------------------------------------------

#=======================================================================================================================

def main(input_path, dir):
    """
    Runs the simulation and outputs the results to a csv file in the results directory (default: results/)
    
    args:
    input_path (str): the path to the simulation template
    dir (str): the directory to save the simulation results to 
    """
    #----------------------------------------------SETUP SIM OUTPUT FILE----------------------------------------------------
    if not os.path.exists(dir):
        os.makedirs(dir)
    
    run_timestamp = datetime.now().strftime("%Y-%m-%d %H-%M-%S")
    
    #create a .csv file for the simulation parameters
    directory = os.path.join(dir, f"{run_timestamp} sim")
    os.makedirs(directory, exist_ok=True)
    simparamsfilename = "parameters.csv"
    simparamsfilename = os.path.join(directory, simparamsfilename)
    
    #create a .csv file for the body data 
    bodydatafilename = "bodies.csv"
    bodydatafilename = os.path.join(directory, bodydatafilename)
    
    #create a .csv file for the quad data
    quaddatafilename = "tree.csv"
    quaddatafilename = os.path.join(directory, quaddatafilename)
    
    #-------------------------------CONSTRUCT THE 0TH ITERATION FROM A SIMULATION TEMPLATE----------------------------------

    consider_energy = True #set to true to calculate the total energy of the system (this is computationally expensive)
    
    print('\nInserted bodies with parameters: \n')
    with open(input_path, newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        axes_size, min_size, theta, max_iter, dt, save_rate = [float(next(csv_reader)[0]) for _ in range(6)] #read the simulation parameters from the first 5 rows of the csv file
        body_names = {}
        curr_tree = QuadTree(Vec2(0, 0), Vec2(axes_size, axes_size), min_size)
        new_tree = QuadTree(Vec2(0, 0), Vec2(axes_size, axes_size), min_size)

        for row in csv_reader:
            body = Body(
                Vec2(float(row[2]), float(row[3])),  # position
                Vec2(float(row[4]), float(row[5])),  # velocity 
                float(row[6]),  # mass
                float(row[7]),  # radius
                int(row[0])  # id
            )
            body_names[int(row[0])] = str(row[1]) 
            curr_tree.insert(body)
            print(f'{row[1]}\t|  pos: {body.position} m \t | vel: {body.velocity} m/s \n\t| mass: {body.mass:.3g} kg\t | radius: {(body.radius):.3g} m\t\t | mag vel: {Vec2(body.velocity.x, body.velocity.y).magnitude()} m/s\t\n')
    
    #write the simulation parameters to the parameters.csv file
    with open(simparamsfilename, 'w', newline='') as file:
        writer = csv.writer(file)
        sim_params = [axes_size, min_size, theta, max_iter, dt, save_rate]
        writer.writerow(sim_params)
    
    if consider_energy:
        curr_tree.calculate_energy()
        initial_net_energy = curr_tree.get_net_energy()
        print(f'\nInitial Total Energy = {initial_net_energy:.5e} J\n')
        
    print('Generating plot of inital conditions...\n')
    print('\nSimulating...\n')

    fig = go.Figure()
    curr_tree.draw(fig, draw_bodies=True, draw_quads=True, draw_com=False, draw_forces=False, draw_velocity=False)
    
    fig.update_layout(
    width=1200,  
    height=1200,  
    xaxis=dict(range=[0, axes_size],showgrid=False,),
    yaxis=dict(range=[0, axes_size],showgrid=False,),
    paper_bgcolor='black',
    plot_bgcolor='black',
    xaxis_title='X Position (m)',
    yaxis_title='Y Position (m)',
    title='Initial Conditions',
    )
    fig.show()
    
    #=========================================================SIMULATE=========================================================

    early_exit_reason = "All bodies ejected from system"

    iter = 0 #starting frame
    mean_runtime = 0
    t_start = time.time()
   
    try:
        while iter < max_iter:
            if curr_tree.size < 1: #if all bodies have been ejected from the system stop the simulation
                break

            t0 = time.time()
            
            curr_tree.calculate_com_approx(iter, save_rate, quaddatafilename)
            curr_tree.calculate_forces(theta)
            curr_tree.integrate(dt, theta, iter, save_rate, new_tree, bodydatafilename)
            
            if consider_energy:
                curr_tree.calculate_energy()
                curr_tree_net_energy = curr_tree.get_net_energy()
            else:
                curr_tree_net_energy = 1
                initial_net_energy = 1
            
            curr_tree = new_tree 
            new_tree = QuadTree(Vec2(0, 0), Vec2(axes_size, axes_size), min_size)
            t1 = time.time()
            
            runtime = t1-t0
            mean_runtime += runtime
            
            bar_length = 100
            progress_percentage = ((iter + 1) / max_iter) * 100
            num_chars = int(bar_length * (progress_percentage / 100))
            
            if iter % save_rate == 0:
                print(
                    f'\033[92m[SIM] (i = {iter}) | n = {curr_tree.size} | t: {runtime:.5f}s | E: {curr_tree_net_energy:.5e} J | dE: {(curr_tree_net_energy - initial_net_energy) / initial_net_energy:.5e} | dt: {dt} s | ['
                    f'{"■" * num_chars}{"⠀" * (bar_length - num_chars)}] {progress_percentage:.3f} %\033[0m',
                    end='\r',
                    flush=True
                )
            else:
                print(
                    f'[SIM] (i = {iter}) | n = {curr_tree.size} | t: {runtime:.5f}s | E: {curr_tree_net_energy:.5e} J | dE: {(curr_tree_net_energy - initial_net_energy) / initial_net_energy:.5e} | dt: {dt} s | ['
                    f'{"■" * num_chars}{"⠀" * (bar_length - num_chars)}] {progress_percentage:.3f} %',
                    end='\r',
                    flush=True
                )
            
            iter += 1
     
    except KeyboardInterrupt:
        early_exit_reason = "Forcefully terminated"

    t_end = time.time()
    
    mean_runtime /= iter
    print(f"\nSimulation Finished | Mean iter runtime: {mean_runtime:.5f}s | Total runtime: {t_end-t_start:.5f}s")
    if iter < max_iter:
        print(f"\nExited Early at iteration {iter} | Reason: {early_exit_reason}.")

if __name__ == "__main__":
    """Run `python nbody --help` for argument information"""

    parser = argparse.ArgumentParser(prog="nbody.py")
    parser.add_argument('-i',   '--input',  type=str,   default="sim templates/testcase-barnes-hut.csv",       help="path to simulation template")
    parser.add_argument('-d',   '--dir',    type=str,   default="results/", help="output directory for simulation results")

    args = parser.parse_args()

    main(
        args.input,
        args.dir
    )       

