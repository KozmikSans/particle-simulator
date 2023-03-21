from dataclasses import dataclass
import csv
import math
import random
fi = open("simlog.tsv","w")
simulation_log = csv.writer(fi,delimiter="\t")
K = 8.99 * 10**9
DT = 0.0001
SM_DIST = 6.677 * 10**-9
time = 0

def vector_sum(vectors:list): # return mathematical vector sum of a vector list
    [x,y,z] = 0,0,0
    for v in vectors:
        x+=v[0]
        y+=v[1]
        z+=v[2]
    return [x,y,z]
def vector_through_two_coordinates(pos1,pos2): # analytical geometry basic, calculating vector given two points
    return pos2[0]-pos1[0],pos2[1]-pos1[1],pos2[2]-pos1[2]
def vect_to_distance(vect): # use of 3d pythagoran theorem to return the length of a vector
    total = vect[0]**2+vect[1]**2+vect[2]**2
    return math.sqrt(total)
def vect_divide(vect,num): # return mathematical division of 2 vectors
    return [vect[0]/num,vect[1]/num,vect[2]/num]
def vect_multiply(vect,num): # returns mathematical multiplication of 2 vectors
    return [vect[0]*num,vect[1]*num,vect[2]*num]

def compute_accel_for_particles(particle_main,particle): # returns accelleration vector of a particle-particle system
    particle_vector = vector_through_two_coordinates(particle_main.position,particle.position)
    v_len = vect_to_distance(particle_vector)
    coulomb_force = K*(particle_main.charge*particle.charge)/v_len**2
    unit_vector = vect_divide(particle_vector,v_len)
    force_coulomb_3d = vect_multiply(unit_vector,coulomb_force)
    accel = vect_divide(force_coulomb_3d,particle_main.mass)
    return accel

particles = []
@dataclass
class particle:
    id: int
    mass: int
    charge: float
    acceleration: list
    position: list
    def compute_acceleration(self): # cumputes accelleration
        all_vect = [self.acceleration]
        for part in particles:
            if part==self:
                continue
            vect = compute_accel_for_particles(self,part)
            all_vect.append(vect)
        self.acceleration = vector_sum(all_vect) # sum all accelleration vectors gathered to derive final accelleration of particle
        return self.acceleration
    def compute_position(self): # gain position using the accelleration vector
        toadd = [(self.acceleration[0]*(DT**2))/2,(self.acceleration[1]*(DT**2))/2,(self.acceleration[2]*(DT**2))/2]
        self.position = vector_sum([self.position,toadd])
        return self.position
    def __str__(self):
        return f"particle {self.id} => mass: {self.mass}, charge: {self.charge}, acceleration: {self.acceleration}, position: {self.position}"
    
def write_down(): # save simulation data of current iteration
    towrite = [time]
    for x in particles:
        towrite.append(str(x.position)[1:-1].replace(", ",","))
    simulation_log.writerow(towrite)

particles.append(particle(id = 0,mass=1.673 * 10**-27,charge=1,acceleration=[0,0,0],position=[0,0,0]))
particles.append(particle(id =1,mass=1.673 * 10**-27,charge=random.randint(-3,3),acceleration=[0,0,0],position=[0+random.randint(-100,100)*SM_DIST,0+random.randint(-100,100)*SM_DIST,0+random.randint(-100,100)*SM_DIST]))
particles.append(particle(id =2,mass=1.673 * 10**-27,charge=random.randint(-3,3),acceleration=[0,0,0],position=[0+random.randint(-100,100)*SM_DIST,0+random.randint(-100,100)*SM_DIST,0+random.randint(-100,100)*SM_DIST]))
particles.append(particle(id =3,mass=1.673 * 10**-27,charge=random.randint(-3,3),acceleration=[0,0,0],position=[0+random.randint(-100,100)*SM_DIST,0+random.randint(-100,100)*SM_DIST,0+random.randint(-100,100)*SM_DIST]))
particles.append(particle(id =4,mass=1.673 * 10**-27,charge=random.randint(-3,3),acceleration=[0,0,0],position=[0+random.randint(-100,100)*SM_DIST,0+random.randint(-100,100)*SM_DIST,0+random.randint(-100,100)*SM_DIST]))

for x in particles: # setting up particle information in the log
    simulation_log.writerow([str(x)])

headers = [f"part_{x.id}" for x in particles]
simulation_log.writerow(headers)

while time<0.005: # main simulation loop
    for x in particles: # first compute the accellerations of all particles
        x.compute_acceleration()
    for x in particles: # then their positions
        x.compute_position()
    write_down()
    time+=DT
