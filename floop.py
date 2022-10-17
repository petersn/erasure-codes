#field_size = 257
#field_neg = lambda a: -a % field_size
#field_add = lambda a, b: (a + b) % field_size
#field_mul = lambda a, b: (a * b) % field_size
#field_inv = lambda a: pow(a, field_size - 2, field_size)

field_size = 256
field_neg = lambda a: a
field_add = lambda a, b: a ^ b

def field_mul(a, b):
    p = 0
    for _ in range(8):
        if b & 1: p ^= a
        a <<= 1
        if a & 0x100: a ^= 0x11b
        b >>= 1
    return p

print("A")
a = 1
s = set()
for _ in range(2560):
	a = field_mul(a, 2)
	s.add(a)
print(len(s))
exit()

# field_size = 4
# field_neg = lambda a: a
# field_add = lambda a, b: a ^ b

# def field_mul(a, b):
#     p = 0
#     for _ in range(2):
#         if b & 1: p ^= a
#         a <<= 1
#         if a & 0b100: a ^= 0b111
#         b >>= 1
#     return p




inv_table = {i: [j for j in range(field_size) if field_mul(i, j) == 1][0] for i in range(1, field_size)}
field_inv = lambda a: inv_table[a]

def poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i in range(len(p)):
        for j in range(len(q)):
            r[i + j] = field_add(r[i + j], field_mul(p[i], q[j]))
    return r

def eval_poly(poly, x):
    r, a = 0, 1
    for coef in poly:
        r = field_add(r, field_mul(coef, a))
        a = field_mul(a, x)
    return r

def lagrange_eval(pos, deg, eval_x):
    r = 1
    for i in range(deg):
        if i == pos:
            continue
        r = field_mul(r, eval_x ^ i)
        r = field_mul(r, field_inv(pos ^ i))
    return r

for x in range(5):
    for eval_x in range(5, 8):
        print("%02x" % lagrange_eval(x, 5, eval_x), end=" ")
    print()
print()
exit()

def lagrange_poly(pos, deg):
    poly = [1]
    for i in range(deg):
        if i == pos:
            continue
        poly = poly_mul(poly, [field_neg(i), 1])
        f = pos ^ i
        poly = poly_mul(poly, [field_inv(f)])
    #poly = poly_mul(poly, [field_inv(eval_poly(poly, pos))])
    return poly


def lagrange_poly(pos, all_pos):
    poly = [1]
    for other in all_pos:
        if other == pos:
            continue
        poly = poly_mul(poly, [field_neg(other), 1])
        f = pos ^ other
        poly = poly_mul(poly, [field_inv(f)])
    #poly = poly_mul(poly, [field_inv(eval_poly(poly, pos))])
    return poly

nums = [0, 1]#list(range(2, 4))
for pos in nums:
	p = lagrange_poly(pos, nums)
	for x in range(len(nums) + 2):
		print(eval_poly(p, x), end=" ")
	print()


exit()



for x in range(5):
    p = lagrange_poly(x, 5)
    for eval_x in range(5, 8):
        print("%02x" % eval_poly(p, eval_x), end=' ')
    print()

def reed_solomon_encode(data, redundancy):
    extra = [0] * redundancy
    for i, val in enumerate(data):
        poly = poly_mul([val], lagrange_poly(i, len(data)))
        for j in range(redundancy):
            extra[j] = field_add(extra[j], eval_poly(poly, len(data) + j))
    return extra

data = [0x00, 0x01, 0x02, 0x03, 0x04]
for i in range(5):
    data[i] += 1
rs = reed_solomon_encode(data, 3)
print()
print(" ".join("%02x" % x for x in rs))

