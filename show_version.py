import itertools

# Remainder of p1 / p2 in GF(2)[x], lowest order bit is the constant term
def poly_rem(p1: int, p2: int) -> bool:
    while p2.bit_length() <= p1.bit_length():
        p1 ^= p2 << (p1.bit_length() - p2.bit_length())
    return p1

def poly_gcd(p1: int, p2: int) -> int:
    while p1 and p2:
        p1, p2 = p2, poly_rem(p1, p2)
    return p1 | p2

def check_irreducible(poly: int) -> bool:
    deg = poly.bit_length() - 1
    return poly_gcd(poly, 1 << 2**deg ^ 2) == poly and all(
        poly_gcd(poly, 1 << 2**(deg // d) ^ 2) == 1
        for d in range(2, deg + 1) if deg % d == 0
    )

for d in [2, 3, 4, 5, 6, 7, 8]:
    for i in range(2**d):
        poly = i | (1 << d)
        if check_irreducible(poly):
            print(f"// {poly:02x} is irreducible of degree {d}")
            break

def poly_mul(a: int, b: int) -> int:
    p = 0
    for i in range(b.bit_length()):
        if b & (1 << i):
            p ^= a << i
    return p

def gf_mul(poly: int, a: int, b: int) -> int:
    return poly_rem(poly_mul(a, b), poly)

def gf_inv(poly: int, a: int) -> int:
    r = 1
    for _ in range(poly.bit_length() - 2):
        a = gf_mul(poly, a, a)
        r = gf_mul(poly, r, a)
    return r

g = 3
poly = 0x13
# Make a log table in GF(16)
log = [0] * 16
c = 1
for i in range(16):
    log[i] = c
    c = gf_mul(poly, c, g)
print(log)
exit()

for i in range(1, 16):
    assert check_irreducible(0x13)
    print(gf_inv(0x13, i), end=", ")
print()
exit()

def rs_mat(poly: int, input_x: list[int], output_x: list[int]):
    table = [[1] * len(output_x) for _ in range(len(input_x))]
    for i in range(len(input_x)):
        for j in range(len(output_x)):
            for p in range(len(input_x)):
                if p != i:
                    table[i][j] = gf_mul(poly, table[i][j], gf_mul(
                        poly,
                        output_x[j] ^ input_x[p],
                        gf_inv(poly, input_x[i] ^ input_x[p]),
                    ))
    return table

def encode(poly: int, table: list[list[int]], data: list[int]):
    output = [0] * len(table[0])
    for i in range(len(data)):
        for j in range(len(output)):
            output[j] ^= gf_mul(poly, table[i][j], data[i])
    return output



def synthesize(poly: int, input_x: list[int], output_x: list[int]):
    deg = poly.bit_length() - 1
    table = rs_mat(poly, input_x, output_x)
    cost = 0
    #print(output_x, "Table:", table)
    #print(f"// data has {len(input_x)} pointers, each pointing to {deg} uint64_ts")
    #print(f"// out has {len(output_x)} pointers, each pointing to {deg} uint64_ts")
    #print("void encode(uint64_t **data, uint64_t **out) {")
    for i, row in enumerate(table):
        for j, c in enumerate(row):
            #print(end="    ")
            for bit_in in range(deg):
                effect = gf_mul(poly, c, 1 << bit_in)
                for bit_out in range(deg):
                    if effect & (1 << bit_out):
                        cost += 1
                        #print(f"    out[{j}][{bit_out}] ^= data[{i}][{bit_in}];")
            #print()
    #print("}")
    #print("// cost:", cost)
    return cost

for have in itertools.combinations(range(8), 5):
    want = [i for i in range(8) if i not in have]
    assert check_irreducible(0b10011)
    cost = synthesize(0b10011, have, want)
    print(f"{have} -> {want}: {cost}")
exit()

# === Tests ===

for poly in range(1, 100):
    continue
    if not check_irreducible(poly):
        continue
    print("Testing inverses in GF(2)[x]/" + bin(poly))
    for a in range(1, 2**(poly.bit_length() - 1)):
        assert gf_mul(poly, a, gf_inv(poly, a)) == 1

for p in range(256):
    poly = p | 0x100
    if not check_irreducible(poly):
        continue
    print("Testing RS matrix in GF(2)[x]/" + bin(poly))
    mat = rs_mat(poly, [0, 1, 2, 3, 4], [5, 6, 7])
    data = [1, 2, 3, 4, 5]
    data += encode(poly, mat, data)
    print("Data:", data)
    for have in itertools.combinations(range(len(data)), 5):
        want = [i for i in range(len(data)) if i not in have]
        print("  Have:", have, "Want:", want)
        rec_mat = rs_mat(poly, list(have), list(want))
        rec_data = [data[i] for i in have]
        print("  Rec data:", rec_data)
        recovered = encode(poly, rec_mat, rec_data)
        for i in range(len(recovered)):
            assert recovered[i] == data[want[i]]
        print("  Recovered:", recovered)
