from decimal import Decimal, getcontext

from Vector import Vector

getcontext().prec = 30


class Line(object):

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = all_zeros
        self.normal_vector = Vector(normal_vector)

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()


    def set_basepoint(self):
        try:
            n = self.normal_vector.coordinates
            c = self.constant_term
            basepoint_coords = ['0']*self.dimension

            initial_index = Line.first_nonzero_index(n)
            initial_coefficient = n[initial_index]

            basepoint_coords[initial_index] = c/initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e


    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector.coordinates

        try:
            initial_index = Line.first_nonzero_index(n)
            terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output


    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)


    def is_parallel_to(self, l):
        n1 = self.normal_vector.normalized()
        n2 = l.normal_vector.normalized()
        if n1 == n2:
            return True
        else:
            return False


    def is_equal_to(self, l):
        if self.is_parallel_to(l):
            k1 = self.constant_term
            b1 = self.normal_vector.coordinates[1]
            p1 = (0,k1/b1)
            k2 = l.constant_term
            b2 = l.normal_vector.coordinates[1]
            p2 = (0,k2/b2)
            if p1 == p2:
                return True
        return False


    def intersection_to(self, l):
        if self.is_parallel_to(l):
            if self.is_equal_to(l):
                print('L1 is equal to L2')
            else:
                print('L1 and L2 are parallel and have no intersections')
        else:
            A = self.normal_vector.coordinates[0]
            B = self.normal_vector.coordinates[1]
            C = l.normal_vector.coordinates[0]
            D = l.normal_vector.coordinates[1]
            k1 = self.constant_term
            k2 = l.constant_term
            x = (D*k1 - B*k2)/(A*D - B*C)
            y = ((-C*k1) + (A*k2))/(A*D - B*C)
            print("Intersection: ({},{})".format(x, y))


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


def test_sets():
    pass
    # Test 1: L1 is parallel to L2

    # Test 2: L1 is equal to L2

    # Test 3: Points of intersection of L1 and L2