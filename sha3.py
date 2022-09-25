import numpy as np
from math import log2

def dbg_print_3d(A):
    # def conv(arr):
    #     n = 0
    #     for i in range(0, 64):
    #         n |= (int(arr[i]) << i)
    #     return n
    # arr = list([list((conv(y) for y in x)) for x in A])
    # print(arr)
    pass

def keccak_p(nr: int, S: np.zeros):
    b = S.size

    w = b // 25
    l = int(log2(w))
    A = np.zeros((5, 5, w), dtype=bool)
    for x in range(5):
        for y in range(5):
            for z in range(w):
                A[x,y,z] = S[w*(5*y + x) + z]
    
    # here goes rounds
    def theta(A: np.zeros):
        dbg_print_3d(A)
        C = np.zeros((5, w), bool)
        D = np.zeros((5, w), bool)
        ret = np.zeros((5, 5, w), bool)
        for x in range(5):
            for z in range(w):
                C[x,z] = A[x, 0, z] ^ A[x, 1, z] ^ A[x, 2, z] ^ A[x, 3, z] ^ A[x, 4, z]
        for x in range(5):
            for z in range(w):
                D[x, z] = C[(x-1) % 5, z] ^ C[(x+1) % 5, (z-1) % w]
        for x in range(5):
            for y in range(5):
                for z in range(w):
                    ret[x,y,z] = A[x,y,z] ^ D[x,z]
        return ret
    # to - включительно
    def rho(A: np.zeros):
        # print("theta -> A -> rho")
        dbg_print_3d(A)
        ret = np.zeros((5, 5, w), bool)
        x, y = 1, 0
        for z in range(w):
            ret[0,0,z] = A[0, 0, z]
        for t in range(24):
            for z in range(w):
                ret[x,y,z] = A[x, y, (z-(t+1)*(t+2)//2)%w]
            x, y, = y, (2*x+3*y) % 5
        return ret
    
    def pi(A: np.zeros):
        dbg_print_3d(A)
        ret = np.zeros((5, 5, w), bool)
        for x in range(5):
            for y in range(5):
                for z in range(w):
                    ret[x,y,z] = A[(x+3*y)%5,x,z]
        return ret

    def chi(A: np.zeros):
        dbg_print_3d(A)
        ret = np.zeros((5, 5, w), bool)
        for x in range(5):
            for y in range(5):
                for z in range(w):
                    ret[x,y,z] = A[x,y,z] ^ ((A[(x+1)%5, y, z] ^ 1) & A[(x+2)%5,y,z])
        return ret

    def rc(t: int) -> bool:
        if t % 255 == 0:
            return True
        r = [True] + [False]*7
        for i in range(1, t%255 + 1):
            r = [False] + r
            r[0] = r[0] ^ r[8]
            r[4] = r[4] ^ r[8]
            r[5] = r[5] ^ r[8]
            r[6] = r[6] ^ r[8]
            r = r[:8]
        return r[0]
    
    def iota(A: np.zeros, i: int):
        ret = np.zeros((5, 5, w), bool)
        for x in range(5):
            for y in range(5):
                for z in range(w):
                    ret[x,y,z] = A[x,y,z]
        RC = np.zeros((w), dtype=bool)
        for j in range (l+1):
            RC[2**j - 1] = rc(j + 7*i)
        for z in range(w):
            ret[0, 0, z] = ret[0, 0, z] ^ RC[z]
        return ret


    for i in range(12 + 2*l - nr, 12 + 2*l - 1 + 1):
        A = iota(chi(pi(rho(theta(A)))), i)
        dbg_print_3d(A)
        # exit(0)
        # print(''.join(map(lambda a: '1' if a else '0', A.flatten())) + '\n')

    S_ = np.zeros((b), dtype=bool)

    for x in range(5):
        for y in range(5):
            for z in range(w):
                S_[w*(5*y + x) + z] = A[x,y,z]
    return S_

def pad101(x: int, m: int):
    j = (-m - 2) % x
    return np.array([True] + [False] * j + [True], bool)



# SHA3-512(M) = KECCAK[1024] (M || 01, 512)
# KECCAK[c] (N, d) = SPONGE[KECCAK-p[1600, 24], pad10*1, 1600 – c] (N, d)

def sponge(r: int, N: np.array, d: int):
    # print(f'r = {r}')
    padding = pad101(r, N.size)
    P = np.concatenate((N, padding))
    n = P.size//r
    b=1600
    c = b - r
    PS = np.split(P, n)
    S = np.zeros((b), bool)
    for i in range(n-1+1):
        S = keccak_p(24, np.logical_xor(S, np.concatenate((PS[i], np.zeros(c, bool)))))
        # print(''.join(map(lambda a: '1' if a else '0', S)))
    Z = np.zeros((0), bool)
    while True:
        Z = np.concatenate((Z, S[:r]))
        if d <= Z.size:
            return Z[:d]
        S = keccak_p(24, S)

def sha3_512(b: bytes):
    M = np.zeros((len(b)*8 + 2), bool)
    for offset in range(len(b)):
        for i in range(8):
            M[offset * 8 + i] = ((b[offset] >> i) & 0b1 ) != 0
    M[-1] = True
    M[-2] = False
    c=1024
    R = sponge(1600 - c, M, 512)
    # print(R.size)
    res = bytearray(512//8)
    for offset in range(512//8):
        n = 0
        for i in range(8):
            n = n | ((1 if R[offset*8 + i] else 0) << i)
        res[offset] = n
    return(res)


import hashlib

def test_sha(s):
    my = sha3_512(s).hex()
    from_hashlib = hashlib.sha3_512(s).hexdigest()
    print()
    print(f'data    = {s}')
    print(f'my      = {my}')
    print(f'hashlib = {from_hashlib}')
    print(f'match   = {my == from_hashlib}')

test_sha(b'hello'*1000)



