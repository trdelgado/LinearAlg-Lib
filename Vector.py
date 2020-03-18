import math

class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(coordinates)
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


    def time_scaler(self, s):
        new_coordinates = [s*e for e in self.dimension]
        return Vector(new_coordinates)


    def magnitude(self):
        magnitude = 0.
        for e in self.coordinates:
            magnitude += e**2
        return magnitude**(1/2)


    def normalized(self):
        try:
            mag = self.magnitude()
            return self.time_scaler(1./mag)
        except ZeroDivisonError:
            raise Except('Cannot normalize the zero vector')
        except TypeError:
            raise TypeError('The coordinates must be iterable')


    def dot_product(self, v):
        assert self.dimension == v.dimension
        return sum([x*y for x, y in zip(self.coordinates, v.coordinates)])


    def radians(self, v):
        assert self.dimension == v.dimension
        try:
            mag1 = self.magnitude()
            mag2 = v.magnitude()
            mag_product = mag1*mag2
            dot_product = self.dot_product(v)
            return math.acos(dot_product/mag_product)
        except ZeroDivisionError:
        	raise TypeError('Vectors must be nonzero')

    def degrees(self, v):
        rad = self.radians(v)
        return (rad*180.)/math.pi