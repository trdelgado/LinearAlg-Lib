from math import sqrt, sin, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'No unique parallel component vector'
    ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG = "Only defined in two three dimensions"

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


    def is_parallel_to(self, v, tolerance=1e-5):
        return (self.is_zero()
            or v.is_zero()
            or (abs(self.angle_with(v)) < tolerance)
            or (abs(self.angle_with(v) - pi) < tolerance))


    def is_orthogonal_to(self, v, tolerance=1e-10):
        return abs(self.dot(v)) < tolerance


    def component_parallel_to(self, basis):
        try:
          u = basis.normalized()
          weight = self.dot(u)
          return u.times_scalar(weight)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)
            else:
                raise e


    def component_orthogonal_to(self, basis):
        try:
            projection = self.component_parallel_to(basis)
            return self.minus(projection)

        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e


    def cross(self, v):
        try:
            x = (self.coordinates[1]*v.coordinates[2] - v.coordinates[1]*self.coordinates[2])
            y = -(self.coordinates[0]*v.coordinates[2] - v.coordinates[0]*self.coordinates[2])
            z = (self.coordinates[0]*v.coordinates[1] - v.coordinates[0]*self.coordinates[1])
            return Vector([x,y,z])

        except ValueError as e:
            msg = str(e)
            if msg == 'need more than 2 values to unpack':
                self_embedded_in_R3 = Vector(self.coordinates + ('0',))
                v_embeeded_in_R3 = Vector(v.coordinates + ('0',))
                return self_embedded_in_R3.cross(v_embeeded_in_R3)
            elif (msg == 'too many values to unpack' or 
                  msg == 'need more than 1 value to unpack'):
                raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)
            else:
                raise e


    def area_of_parallelogram_with(self, v):
        cross_product = self.cross(v)
        return cross_product.magnitude()


    def area_of_triangle_with(self, v):
        return (self.area_of_parallelogram_with(v)/2.0)


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

    # Test 9: Component Vectors
    v = Vector([3.039,1.879])
    w = Vector([0.825,2.036])
    u = Vector([1.083,2.672])
    result = v.component_parallel_to(w)
    result = Vector([round(x,3) for x in result.coordinates])
    assert(result == u)

    v = Vector([-9.88,-3.264,-8.159])
    w = Vector([-2.155,-9.353,-9.473])
    u = Vector([-8.35,3.376,-1.434])
    result = v.component_orthogonal_to(w)
    result = Vector([round(x,3) for x in result.coordinates])
    assert(result == u)

    v = Vector([3.009,-6.172,3.692,-2.51])
    w = Vector([6.404,-9.144,2.759,8.718])
    u1 = Vector([1.969,-2.811,0.848,2.680])
    u2 = Vector([1.040,-3.361,2.844,-5.190])
    result1 = v.component_parallel_to(w)
    result1 = Vector([round(x,3) for x in result1.coordinates])
    result2 = v.component_orthogonal_to(w)
    result2 = Vector([round(x,3) for x in result2.coordinates])
    assert(result1 == u1)
    assert(result2 == u2)

    print("Passed test 9: Component Vectors")

    # Test 10: Cross Products, Areas of Parallelograms and Triangles 
    v = Vector([8.462,7.893,-8.187])
    w = Vector([6.984,-5.975,4.778])
    result = v.cross(w)
    result = Vector([round(x,3) for x in result.coordinates])
    u = Vector([-11.205,-97.609,-105.685])
    assert(result == u)

    v = Vector([-8.987,-9.838,5.031])
    w = Vector([-4.268,-1.861,-8.866])
    result = round(v.area_of_parallelogram_with(w),3)
    u = 142.122
    assert(result == u)

    v = Vector([1.5,9.547,3.691])
    w = Vector([-6.007,0.124,5.772])
    result = round(v.area_of_triangle_with(w),3)
    u = 42.565
    assert(result == u)

    print("Passed test 10: Cross Products and Areas")