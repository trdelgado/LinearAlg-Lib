from math import sqrt, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([round(Decimal(x),10) for x in coordinates])
            self.dimension = len(self.coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates


    def plus(self, v):
        new_coordinates = [Decimal(x+y) for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)


    def minus(self, v):
        new_coordinates = [Decimal(x-y) for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)


    def times_scalar(self, c):
        new_coordinates = [Decimal(c)*x for x in self.coordinates]
        return Vector(new_coordinates)


    def magnitude(self):
        return sqrt(sum([x**2 for x in self.coordinates]))


    def normalized(self):
        try:
            magnitude = self.magnitude()
            return self.times_scalar(1./magnitude)

        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)


    def dot(self, v):
        return Decimal(sum([x*y for x,y in zip(self.coordinates, v.coordinates)]))


    def angle_with(self, v, in_degrees=False):
        try:
            u1 = self.normalized()
            u2 = v.normalized()
            val = u1.dot(u2)
            if val > Decimal(1.0):
                val = Decimal(1.0)
            elif val < Decimal(-1.0):
                val = Decimal(-1.0)
            angle_in_radians = acos(val)

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


    def is_zero(self, tolerance=1e-10):
        return self.magnitude() < tolerance


    def is_parallel_to(self, v):
        return (self.is_zero()
            or v.is_zero()
            or (self.angle_with(v) == 0)
            or (self.angle_with(v) == pi))


    def is_orthogonal_to(self, v, tolerance=1e-10):
        return abs(self.dot(v)) < tolerance


def test_sets():
    # Test 1: Addition
    v = Vector([8.218,-9.341])
    w = Vector([-1.129,2.111])
    a = Vector([7.089,-7.230])
    assert(v.plus(w) == a)
    print("Passed test 1: Addition")

    # Test 2: Subtraction
    v = Vector([7.119,8.215])
    w = Vector([-8.223,0.878])
    a = Vector([15.342,7.337])
    assert(v.minus(w) == a)
    print("Passed test 2: Subtraction")

    # Test 3: Scaler Multiplication
    v = Vector([1.671,-1.012,-0.318])
    c = 7.41
    a = Vector([12.382,-7.499,-2.356])
    result = v.times_scalar(c)
    result = Vector([round(x,3) for x in result.coordinates])
    assert(result == a)
    print("Passed test 3: Scalar Multiplication")

    # Test 4: Magnitude
    v = Vector([-0.221,7.437])
    m = 7.440
    assert(round(v.magnitude(),3) == m)

    v = Vector([8.813,-1.331,-6.247])
    m = 10.884
    assert(round(v.magnitude(),3) == m)
    print("Passed test 4: Vector Magnitude")

    # Test 5: Direction
    v = Vector([5.581,-2.136])
    u = Vector([0.934,-0.357])
    result = v.normalized()
    result = Vector([round(x,3) for x in result.coordinates])
    assert(result == u)

    v = Vector([1.996,3.108,-4.554])
    u = Vector([0.340,0.530,-0.777])
    result = v.normalized()
    result = Vector([round(x,3) for x in result.coordinates])
    assert(result == u)
    print("Passed test 5: Vector Direction")

    # Test 6: Dot Products
    v = Vector([7.887,4.138])
    w = Vector([-8.802,6.776])
    s = round(Decimal(-41.382),3)
    assert(round(v.dot(w),3) == s)

    v = Vector([-5.955,-4.904,-1.874])
    w = Vector([-4.496,-8.755,7.103])
    s = round(Decimal(56.397),3)
    assert(round(v.dot(w),3) == s)
    print("Passed test 6: Dot Products")

    # Test 7: Angles
    v = Vector([3.183,-7.627])
    w = Vector([-2.668,5.319])
    r = 3.072
    assert(round(v.angle_with(w),3) == r)

    v = Vector([7.35,0.221,5.188])
    w = Vector([2.751,8.259,3.985])
    d = 60.276
    assert(round(v.angle_with(w,in_degrees=True),3) == d)
    print("Passed test 7: Angles")

    # Test 8: Parallel and Orthogonal
    v = Vector([-7.579,-7.88])
    w = Vector([22.737,23.64])
    assert(v.is_parallel_to(w) == True)
    assert(v.is_orthogonal_to(w) == False)

    v = Vector([-2.029,9.97,4.172])
    w = Vector([-9.231,-6.639,-7.245])
    assert(v.is_parallel_to(w) == False)
    assert(v.is_orthogonal_to(w) == False)

    v = Vector([-2.328,-7.284,-1.214])
    w = Vector([-1.821,1.072,-2.94])
    assert(v.is_parallel_to(w) == False)
    assert(v.is_orthogonal_to(w) == True)

    v = Vector([2.118,4.827])
    w = Vector([0,0])
    assert(v.is_parallel_to(w) == True)
    assert(v.is_orthogonal_to(w) == True)
    print("Passed test 8: Parallel and Orthogonal Vectors")