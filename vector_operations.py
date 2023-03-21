import math


def vector_sum(vectors):
    """return mathematical vector sum of a vector list"""
    [x, y, z] = 0, 0, 0
    for v in vectors:
        x += v[0]
        y += v[1]
        z += v[2]
    return [x, y, z]


def vector_through_two_coordinates(pos1, pos2):
    """analytical geometry basic, calculating vector given two points"""
    return pos2[0]-pos1[0], pos2[1]-pos1[1], pos2[2]-pos1[2]


def vect_to_distance(vect):
    """use of 3d pythagoran theorem to return the length of a vector"""
    total = vect[0]**2+vect[1]**2+vect[2]**2
    return math.sqrt(total)


def vect_divide(vect, num):
    """return mathematical division of 2 vectors"""
    return [vect[0]/num, vect[1]/num, vect[2]/num]


def vect_multiply(vect, num):
    """returns mathematical multiplication of 2 vectors"""
    return [vect[0]*num, vect[1]*num, vect[2]*num]
