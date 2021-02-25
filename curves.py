class Point:
    def __init__(self, x: float, y: float, curve):
        self.x = x
        self.y = y

        # if not curve.check_point(x, y):
        #     raise ValueError("point not on the curve!")

        if not isinstance(curve, WeierstrassNormalCurve):
            raise NotImplementedError()

        self.curve = curve

    def __str__(self) -> str:
        return "({}, {})".format(self.x, self.y)

    def __eq__(self, other) -> bool:
        if self.x == other.x and self.y == other.y:
            return True
        else:
            return False

    def __neg__(self):
        # Reflect over x-axis
        return Point(-self.x, self.y, self.curve)

    def __add__(self, other):
        R = self.curve.compute_R(self, other)
        return -R

    def __sub__(self, other):
        return self + -other

    def __mul__(self, n: int):
        if not isinstance(n, int):
            raise Exception("can only multiply by int")

        if n == 0:
            return IdentityPoint(self.curve)
        elif n < 0:
            return -self * -n
        else:
            Q = self
            R = self if n & 1 == 1 else IdentityPoint(self.curve)

            i = 2
            while i <= n:
                Q = Q + Q

                if n & i == i:
                    R = Q + R

                i = i << 1

    def __rmul__(self, n: int):
        return self * n

class IdentityPoint(Point):
    def __init__(self, curve):
        self.x = None
        self.y = None
        self.curve = curve

    def __neg__(self):
        return self

    def __str__(self) -> str:
        return "IdentityPoint"

    def __add__(self, other):  # since it is the identity element in the group
        return other

    def __mul__(self, n: int):
        return self


class WeierstrassNormalCurve:
    def __init__(self, a: float, b: float):
        self.a = a
        self.b = b

        if 4 * a**3 + 27 * b**2 == 0:
            raise ValueError("err: Curve is singular!")

    def check_point(self, x: float, y: float) -> bool:
        lhs = y ** 2
        rhs = x ** 3 + self.a * x + self.b
        return lhs == rhs

    def compute_R(self, P: Point, Q: Point) -> Point:
        """Find point R intersecting line formed by points P and Q"""

        R = Point(0, 0, self)

        if P.x != Q.x:
            m = (P.y - Q.y) / (P.x - Q.x)
        elif P.x == Q.x and P.y == Q.y:
            m = (3 * P.x**2 + self.a) / (2 * P.y)
        elif P.x == Q.x:
            return IdentityPoint(self)

        R.x = m**2 - P.x - Q.x
        R.y = P.y + m * (R.x - P.x)

        return R


if __name__=="__main__":
    # parts from https://jeremykun.com/2014/02/24/elliptic-curves-as-python-objects/
    curve = WeierstrassNormalCurve(-7, 10)
    assert Point(3, -2, curve) == Point(1, 2, curve) + Point(3, 4, curve)
    assert Point(-1, 2, curve) == Point(-1, 4, curve) + Point(1, 2, curve)
    P = Point(1, 2, curve)
    assert Point(1, 4, curve) == P + P
    assert P == IdentityPoint(curve) + P
    3 * P
    # the below doesn't pass yet
    # assert IdentityPoint(curve) == P + -P
