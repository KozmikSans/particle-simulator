from simulation import Simulation


if __name__ == "__main__":
    s = Simulation("simlog", 300)

    for i in range(10):
        s.generate_particle()

    s.simulate()
    s.show()
    s.save(path="saved/10_particles_test", video=True, frames=True, log=True)
