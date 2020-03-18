from math import sqrt, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates


    def plus(self, v):
        assert self.dimension == v.dimension
        new_coordinates = [x+y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)


    def minus(self, v):
        assert self.dimension == v.dimension
        new_coordinates = [x-y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)


    def times_scaler(self, c):
        new_coordinates = [Decimal(c)*x for x in self.coordinates]
        return Vector(new_coordinates)


    def magnitude(self):
        return sqrt(sum([Decimal(x)**2 for x in self.coordinates]))


    def normalized(self):
        try:
            magnitude = self.magnitude()
            return self.times_scaler(1./magnitude)
        except ZeroDivisonError:
            raise Except(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)


    def dot(self, v):
        assert self.dimension == v.dimension
        return sum([x*y for x, y in zip(self.coordinates, v.coordinates)])


    def angle_with(self, v, in_degrees=False):
        try:
            u1 = self.normalized()
            u2 = v.normalized()
            angle_in_radians = acos(u1.dot(u2))

            if in_degrees:
                degrees_per_radian = 180./pi
                return angle_in_radians * degrees_per_radian
            else:
                return angle_in_radians
        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception("Cannot compute an angle with the zero vector")
            else:
                raise e

    def parallel_with(self, v):
        u1 = self.normalized()
        u2 = v.normalized()
        if u1.coordinates == u2.coordinates:
        	return True
        else:
        	return False

    def orthogonal_with(self, v):
        pass
