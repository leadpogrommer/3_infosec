import numpy as np
import numba
import time



# rotate (positive = left)
@numba.jit
def rot64(a, n):
    a = np.uint64(a)
    n = np.uint64(n)
    return np.uint64((a >> np.uint64((64 - (n % 64)))) + (a << np.uint64(n % 64)))


@numba.jit
def keccak_p(nr: int, S: np.zeros):
    b = S.size * 64

    w = b // 25
    l = 6 # int(log2(w))
    A = np.zeros((5, 5), dtype=np.uint64)
    for x in range(5):
        for y in range(5):
            A[x,y] = S[(5*y + x)]
    
    # here goes rounds
    def theta(A: np.zeros):

        C = np.zeros((5), np.uint64)
        D = np.zeros((5), np.uint64)
        for x in range(5):
                C[x] = A[x, 0] ^ A[x, 1] ^ A[x, 2] ^ A[x, 3] ^ A[x, 4]
        for x in range(5):
            D[x] = C[(x-1) % 5] ^ np.uint64(rot64(C[(x+1) % 5], 1))
        for x in range(5):
            for y in range(5):
                    A[x,y] = A[x,y] ^ D[x]

        return A
    # to - включительно
    def rho(A: np.zeros):

        ret = np.zeros((5, 5), np.uint64)
        x, y = 1, 0
        ret[0,0] = A[0, 0]
        for t in range(24):
            ret[x,y] = rot64(A[x, y], (t+1)*(t+2)//2)
            x, y, = y, (2*x+3*y) % 5

        return ret
    
    def pi(A: np.zeros):

        ret = np.zeros((5, 5), np.uint64)
        for x in range(5):
            for y in range(5):
                ret[x,y] = A[(x+3*y)%5,x]
        return ret

    def chi(A: np.zeros):

        ret = np.zeros((5, 5), np.uint64)
        for x in range(5):
            for y in range(5):
                ret[x,y] = A[x,y] ^ ((~A[(x+1)%5, y]) & A[(x+2)%5,y])

        return ret

    def rc(t: int) :
        if t % 255 == 0:
            return 1
        r = [1] + [0]*7
        for i in range(1, t%255 + 1):
            r = [0] + r
            r[0] = r[0] ^ r[8]
            r[4] = r[4] ^ r[8]
            r[5] = r[5] ^ r[8]
            r[6] = r[6] ^ r[8]
            r = r[:8]
        return r[0]
    
    def iota(A: np.zeros, i: int):

        RC = np.uint64(0)
        for j in range (l+1):
            RC |= np.uint64((rc(j + 7*i) << (2**j - 1)))
        A[0, 0] = A[0, 0] ^ RC

        return A


    for i in range(12 + 2*l - nr, 12 + 2*l - 1 + 1):
        A = iota(chi(pi(rho(theta(A)))), i)

    S_ = np.zeros((b//64), dtype=np.uint64)

    for x in range(5):
        for y in range(5):
                S_[(5*y + x)] = A[x,y]
    return S_

@numba.jit
def pad101(x: int, m: int):
    j = (-m - 2) % x
    return np.array([1] + [0] * j + [1], np.uint8)



# SHA3-512(M) = KECCAK[1024] (M || 01, 512)
# KECCAK[c] (N, d) = SPONGE[KECCAK-p[1600, 24], pad10*1, 1600 – c] (N, d)

@numba.jit(nopython=True)
def sponge(r, N, d):
    # print(f'r = {r}')
    padding = pad101(r, N.size)
    P_bits = np.concatenate((N, padding))
    n = P_bits.size//r
    b=1600
    c = b - r

    P = np.zeros((len(P_bits)//64), np.uint64)

    for y in range(len(P)):
        tmp = 0
        for i in range(64):
            tmp |= (P_bits[y*64 + i]) << i
        P[y] = tmp

    # now we are using only uint64
    PS = np.split(P, n)
    S = np.zeros((b//64), np.uint64)
    print('started iterations')
    for i in range(n-1+1):
        S = keccak_p(24, np.bitwise_xor(S, np.concatenate((PS[i], np.zeros(c//64, np.uint64)))))
    Z = np.zeros((0), np.uint64)
    while True:
        Z = np.concatenate((Z, S[:r//64]))
        if d//64 <= Z.size:
            print('finished iterations')
            return Z[:d]
        S = keccak_p(24, S)


def sha3_512(b: bytes):
    M = np.zeros((len(b)*8 + 2), np.uint8)
    for offset in range(len(b)):
        for i in range(8):
            M[offset * 8 + i] = ((b[offset] >> i) & 0b1 )
    M[-1] = 1
    M[-2] = 0
    c=1024
    R = sponge(1600 - c, M, 512)
    # print(R.size)
    res = R[:512//64].tobytes()
    # for i in R[:512//64]:
        # res += bytearray(int(i).to_bytes(8, byteorder='little', signed=False))
    return(res)


import hashlib

def test_sha(s):
    start_time = time.time()
    my = sha3_512(s).hex()
    finish_time = time.time()
    
    from_hashlib = hashlib.sha3_512(s).hexdigest()
    print(f'time    = {finish_time - start_time}')
    print(f'my      = {my}')
    print(f'hashlib = {from_hashlib}')
    print(f'match   = {my == from_hashlib}')

test_sha(b'E'*1024*1024)

# TODO: run on 1 meg of input (should take approximately 85 minutes)