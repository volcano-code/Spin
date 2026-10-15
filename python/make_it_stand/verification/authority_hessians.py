"""Differentiate independent rational first derivatives; no native kernels."""
from ._spin2014_authority import _arithmetic, _independent_authority_physics_values
from fractions import Fraction


class _Dual:
    def __init__(self, value, derivative):
        self.value=Fraction(value) if type(value) is int else value
        self.derivative=tuple(derivative)

    def cast(self, other):
        return other if isinstance(other,_Dual) else _Dual(other,(0,)*10)

    def __add__(self, other):
        other=self.cast(other)
        return _Dual(self.value+other.value,[a+b for a,b in zip(self.derivative,other.derivative)])
    __radd__=__add__

    def __neg__(self):
        return _Dual(-self.value,[-x for x in self.derivative])

    def __sub__(self, other): return self+-self.cast(other)
    def __rsub__(self, other): return self.cast(other)+-self

    def __mul__(self, other):
        other=self.cast(other)
        return _Dual(self.value*other.value,[a*other.value+self.value*b for a,b in zip(self.derivative,other.derivative)])
    __rmul__=__mul__

    def __pow__(self, exponent):
        if type(exponent) is not int: raise TypeError('integer powers required')
        if exponent==0: return _Dual(1,(0,)*10)
        return _Dual(self.value**exponent,[exponent*self.value**(exponent-1)*d for d in self.derivative])

    def __truediv__(self, other): return self*self.cast(other)**-1
    def __rtruediv__(self, other): return self.cast(other)*self**-1
    def __abs__(self): return -self if self.value<0 else self
    def __le__(self, other): return self.value<=self.cast(other).value


def independent_authority_hessians(raw_moments, settings, *, arithmetic_backend='FRACTION_V1'):
    """Exact Hessians of existing independently written 2014 formulas.

    Also verifies first derivatives against automatic differentiation of the
    original scalar objectives/constraints. Not a material/derivative gate.
    Returned rationals stay internal, without decimal text round trips.
    """
    rational,read,write=_arithmetic(arithmetic_backend)
    if len(raw_moments)!=10: raise ValueError('ten moments required')
    variables=[_Dual(read(value),[rational(int(i==j)) for j in range(10)])
               for i,value in enumerate(raw_moments)]
    values=_independent_authority_physics_values(variables,settings,
        (rational,lambda value:value,write))
    def as_dual(value):
        return value if isinstance(value,_Dual) else _Dual(value,(rational(0),)*10)
    gradient=[as_dual(x) for x in values['objective_gradient']]
    objective=as_dual(values['objective'])
    if objective.derivative!=tuple(x.value for x in gradient):
        raise ValueError('independent objective first derivatives inconsistent')
    constraint_hessians=[]
    for scalar,row in zip(values['constraints'],values['constraint_jacobian']):
        scalar=as_dual(scalar); row=[as_dual(x) for x in row]
        if scalar.derivative!=tuple(x.value for x in row):
            raise ValueError('independent constraint first derivatives inconsistent')
        constraint_hessians.append([list(x.derivative) for x in row])
    return dict(objective_hessian=[list(x.derivative) for x in gradient],
        constraint_hessians=constraint_hessians)
