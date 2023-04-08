import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import shutil
from dataclasses import dataclass
import random
from vector_operations import *
import csv

# https://www.ucl.ac.uk/~ucfbasc/Theory/lenjon.html - leonard jones potential
#simulating helium ions - sigma,epsilon and mass are helium parameters
SIGMA = 258 * 10**-9
EPSILON = 10.22
HE_MASS = 6.6464731 * 10**-27

K = 8.99 * 10**9
DT = 2e-30
SM_DIST = 6.677 * 10**-9

# why does the LJ potential launch them like 10 kilometers away?
# TODO: this SHIT has to be fixed lol
def compute_accel_for_particles(particle_main, particle):
    """returns accelleration vector of a particle-particle system"""
    particle_vector = vector_through_two_coordinates(
        particle_main.position, particle.position)
    v_len = vect_to_distance(particle_vector)
    coulomb_force = K*-(particle_main.charge*particle.charge)/v_len**2
    LJ_potential = (24*EPSILON/v_len**2)*((2*(SIGMA/v_len)**12)-((SIGMA/v_len)**6)) # so called leonard jones potential
    total_force = coulomb_force - LJ_potential
    unit_vector = vect_divide(particle_vector, v_len)
    force_coulomb_3d = vect_multiply(unit_vector, total_force)
    accel = vect_divide(force_coulomb_3d, particle_main.mass)
    return accel


class Simulation():
    def __init__(self, filename: str = "simlog", n_frames: int = 100) -> None:
        self.filename = filename
        self.particles = []
        self.time = 0
        self.to_time = (n_frames-1) * DT    # there is a frame zero

    def create_particle(self, id: int, mass: int, charge: float, acceleration: list, position: list[float]):
        p = Particle(id, mass, charge, acceleration, position, self.particles)
        self.particles.append(p)

    def generate_particle(self):
        self.create_particle(id=len(self.particles), mass=HE_MASS, charge=random.choice([-1,1,2,-2]), acceleration=[0, 0, 0], position=[
                             0+random.randint(-100, 100)*SM_DIST, 0+random.randint(-100, 100)*SM_DIST, 0+random.randint(-100, 100)*SM_DIST]) 


    def simulate(self):
        if not self.particles:
            raise IndexError("Cannot simulate with no particles")

        def write_down():
            """Writes a row (=one frame) containing time and coordinates of every particle into the simulation log"""
            towrite = [self.time]
            for x in self.particles:
                towrite.append(str(x.position)[1:-1].replace(", ", ","))
            simulation_log.writerow(towrite)

        file = open(f"{self.filename}.tsv", "w")
        simulation_log = csv.writer(file, delimiter="\t")

        for x in self.particles:
            simulation_log.writerow([str(x)])

        headers = ["time"] + [f"particle_{x.id}" for x in self.particles]
        simulation_log.writerow(headers)

        write_down()
        # main simulation loop:
        while self.time < self.to_time:
            for x in self.particles:
                x.compute_acceleration()    # first compute the accellerations of all particles
            for x in self.particles:
                x.compute_position()    # then their positions
            self.time += DT
            write_down()

        file.close()

    def show(self):
        self._animate(path=None, show=True,
                      save_video=False, save_frames=False)

    def save(self, path: str = "saved/other", video: bool = True, frames: bool = True, log: bool = True):
        # create specified folder if it doesn't exists
        if not os.path.isdir(path):
            os.mkdir(path=path)
        if frames:
            os.mkdir(f"{path}/frames")
        # save frames and animation through the _animate function
        self._animate(path=path, show=False,
                      save_video=video, save_frames=frames)
        # if log, copy log to path
        if log:
            shutil.copy(f"{self.filename}.tsv", f"{path}/{self.filename}.tsv")

    def _crunch(self) -> list:
        """takes simlog.tsv (or specified_file.tsv) and puts the data into a nested list"""
        crunched = []

        myfile = open(f"{self.filename}.tsv", "r")
        # exclude the headers
        data = myfile.readlines()[len(self.particles)+1:]

        for line in data:
            frame_info = []
            line = line.strip()  # get rid of the newline char
            line = line.split("\t")
            # frame info in the following format: [time, (x1, y1, 1), ..., (xn, yn, zn)]
            time = float(line.pop(0))
            frame_info.append(time)
            for particle in line:
                particle = particle.split(",")
                particle = [float(coord) for coord in particle]
                frame_info.append(particle)
            crunched.append(frame_info)
        myfile.close()

        return crunched

    def _animate(self, path: str, show=True, save_video=True, save_frames=True) -> None:
        """Works with mpl animation API to create an animation object"""
        frames = self._crunch()

        # Create a figure and 3D axis object
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Add axis labels
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_zlabel('Z-axis')

        def update(frame):
            ax.clear()

            # Set the x, y, and z limits of the axis
            ax.set_xlim(-10e-7, 10e-7)
            ax.set_ylim(-10e-7, 10e-7)
            ax.set_zlim(-10e-7, 10e-7)
            # Draw the x, y, and z axis lines
            ax.plot([-10e-7, 10e-7], [0, 0], [0, 0], color='black')
            ax.plot([0, 0], [-10e-7, 10e-7], [0, 0], color='black')
            ax.plot([0, 0], [0, 0], [-10e-7, 10e-7], color='black')

            frame_info = frames[frame]
            time = frame_info[0]
            # "rotate" the matrix 90 degrees
            frame_info = np.rot90(frame_info[1:], k=1, axes=(1, 0))
            x_data, y_data, z_data = frame_info[0], frame_info[1], frame_info[2]
            # Add a title
            plt.title(f"T: {time}")

            scatter = ax.scatter(x_data, y_data, z_data, marker="o")
            if save_frames:
                fig.savefig(f"{path}/frames/frame_{frame}.png")
            return scatter

        # speed of the video, saved for calculationg the fps
        interval = 25
        animation = FuncAnimation(fig, update, frames=len(frames), interval=interval, blit=False, repeat=False)
        if save_video:
            animation.save(f"{path}/animation.mp4",
                           writer="ffmpeg", fps=1000/interval)
        if show:
            plt.show()


@dataclass
class Particle:
    id: int
    mass: int
    charge: float
    acceleration: list[float]
    position: list[float]
    particles: list["Particle"]

    def compute_acceleration(self):
        """computes acceleration of itself from influences of other particles"""
        all_vect = [self.acceleration]
        for part in self.particles:
            if part == self:
                continue
            vect = compute_accel_for_particles(self, part)
            all_vect.append(vect)
        # sum all accelleration vectors gathered to derive final accelleration of particle:
        self.acceleration = vector_sum(all_vect)
        return self.acceleration

    def compute_position(self):
        """gain position using the accelleration vector"""
        toadd = [(self.acceleration[0]*(DT**2))/2, (self.acceleration[1]
                                                    * (DT**2))/2, (self.acceleration[2]*(DT**2))/2]
        self.position = vector_sum([self.position, toadd])
        return self.position

    def __str__(self):
        return f"particle {self.id} => mass: {self.mass}, charge: {self.charge}, acceleration: {self.acceleration}, position: {self.position}"
