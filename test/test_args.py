

def quadratic(x, a, b, c):
    return a*x*x + b*x + c



p0 = (1, 1.1, 0.07)
print(quadratic(1, p0[0], p0[1], p0[2]))
print(quadratic(1, *p0))
