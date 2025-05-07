import os
import random
import time
from hashlib import shake_256

from montgomery_isogenies.kummer_line import KummerLine
from montgomery_isogenies.kummer_isogeny import KummerLineIsogeny
from montgomery_isogenies.isogenies_x_only import lift_image_to_curve

from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny

from utilities.discrete_log import BiDLP, discrete_log_pari
from utilities.supersingular import torsion_basis, torsion_basis_with_pairing

proof.all(False)    


#set_random_seed(1)
#random.seed(int(1))



def eval_dimtwo_isog(Phi, q, P, Q, ord, RS=None):

    R1 = P
    R2 = Phi.domain()[1](0)
    phiP = Phi(CouplePoint(R1, R2))
    imP = phiP[0]

    R1 = Q
    R2 = Phi.domain()[1](0)
    phiQ = Phi(CouplePoint(R1, R2))
    imQ = phiQ[0]

    R1 = P-Q
    R2 = Phi.domain()[1](0)
    phiPQ = Phi(CouplePoint(R1, R2))
    imPQ = phiPQ[0]

    if (imP - imQ)[0] != imPQ[0]:
        imQ = -imQ


    if RS == None:
        R, S, WP = torsion_basis_with_pairing(Phi.domain()[1], ord)
    else:
        R, S = RS
        WP = R.weil_pairing(S, ord)

    R1 = Phi.domain()[0](0)
    R2 = R
    phiR = Phi(CouplePoint(R1, R2))
    imR = phiR[0]

    R1 = Phi.domain()[0](0)
    R2 = S
    phiS = Phi(CouplePoint(R1, R2))
    imS = phiS[0]

    R1 = Phi.domain()[0](0)
    R2 = R-S
    phiRS = Phi(CouplePoint(R1, R2))
    imRS = phiRS[0]

    if (imR - imS)[0] != imRS[0]:
        imS = -imS

    wp = WP^q
    x, y = BiDLP(imP, imR, imS, ord, ePQ=wp)
    w, z = BiDLP(imQ, imR, imS, ord, ePQ=wp)


    imP = x*R + y*S
    imQ = w*R + z*S

    imP *= q
    imQ *= q

    return imP, imQ


def point_to_xonly(P, Q):
    L = KummerLine(P.curve())

    PQ = P - Q
    xP = L(P[0])
    xQ = L(Q[0])
    xPQ = L(PQ[0])

    return L, xP, xQ, xPQ




def make_prime(N):
    for f in range(1, 1000):
        if (N*f - 1).is_prime():
            return Integer(N*f - 1)
        
def random_unit(modulus):
    while True:
        alpha = ZZ.random_element(modulus)
        if gcd(alpha, modulus) == 1:
            break

    return alpha

def xof_kdf(X, Y):
    """
    A simple key-derivation function which initalised an XOF from two
    finite field elements. A key is intended to then be derived from
    this by requesting len(key) bytes, see meth:xof_encrypt
    """
    X_bytes = X.to_bytes()
    Y_bytes = Y.to_bytes() 
    return shake_256(X_bytes + Y_bytes)

def xof_encrypt(xof, msg):
    """
    Compute the XOR encryption of a message by sampling
    len(msg) bytes from a suitible XOF.
    """
    key = xof.digest(len(msg))
    return bytes(x ^^ y for (x, y) in zip(key, msg))

def random_matrix(modulus):
    while True:
        d1 = ZZ.random_element(modulus)
        d2 = ZZ.random_element(modulus)
        d3 = ZZ.random_element(modulus)
        d4 = ZZ.random_element(modulus)
        if gcd(d1*d4 - d2*d3, modulus) == 1:
            break

    return d1, d2, d3, d4


#bs = [137, 177, 239]
bs = [162, 243, 324]

for i, lambda_ in enumerate([128, 192, 256]):
    #if i == 0 or i == 1:
    #    continue
    a = lambda_ - 2
    b = bs[i]
    # b = ceil(2*lambda_ * log(2)/log(3) * 216 / 256)
    # print(b)
    # b = 137 #ceil(2*lambda_ * log(2)/log(3) * )
    c = floor(lambda_/3 * log(2)/log(5))
    A = Integer(2**a)
    B = Integer(3**b)
    assert B > 2**(2 * lambda_)
    C = 5**c #next_prime(2**(lambda_//2))
    p = make_prime(4*A*B*C)

    x = C

    a = factor(p+1)[0][1] - 2
    A = Integer(2**a)
    print(f"a = {a}, b = {b}, c = {c} ")
    print(factor(p+1))

    # print(a, b, c)

    # print(f"{p = }")
    print(f"{p.nbits() = }")
    # print(f"{p.nbits() / 64. = }")
    # print(f"{p.nbits() / float(lambda_) = }")


    print("pk", ceil(6 * p.nbits() / 8.))
    print("ct", 2 * ceil(6 * p.nbits() / 8.))
    print()
    print("pk_cmp", ceil(5 * p.nbits() / 8.))
    print("ct_cmp", ceil((4 * p.nbits()) /8.) + ceil(3*A.nbits() / 8.))

    FF.<xx> = GF(p)[]
    F.<i> = GF(p^2, modulus=xx^2+1)
    E0 = EllipticCurve(F, [1, 0])
    E0.set_order((p+1)**2)
    PA, QA = torsion_basis(E0, 4*A)
    PB, QB = torsion_basis(E0, B)
    X0, Y0 = torsion_basis(E0, C)


    PQA = PA - QA
    PQB = PB - QB
    XY0 = X0 - Y0

    _E0 = KummerLine(E0)
    _PB = _E0(PB[0])
    _QB = _E0(QB[0])
    _PQB = _E0(PQB[0])

    xPA = _E0(PA[0])
    xQA = _E0(QA[0])
    xPQA = _E0(PQA[0])

    xX_0 = _E0(X0[0])
    xY_0 = _E0(Y0[0])
    xXY_0 = _E0(XY0[0])

        
    def keygenA():
        for _ in range(1000):
            q = random.randint(0, 2^(a-1)-1)
            if q % 2 == 1 and q % 3 != 0 and q % 5 != 0 and (A-q) % 3 != 0 and (A-q) % 5 != 0 :
                break
        else:
            raise ValueError("Could not find suitable q.")

        deg = q
        rhs = deg * (2**a - deg) * B

        upper_bound = ZZ((rhs.nbits() - p.nbits()) // 2 - 2)
        alpha = random_unit(A)
        beta =  random_unit(A)
        gamma = random_unit(B)
        delta = random_unit(C)

        if upper_bound < 0:
            raise ValueError("Degree is too small.")

        P2, Q2, P3, Q3 = PA, QA, PB, QB

        QF = gp.Qfb(1, 0, 1)

        for _ in range(10_000):
            zz = randint(0, ZZ(2**upper_bound))
            tt = randint(0, ZZ(2**upper_bound))
            sq = rhs - p * (zz**2 + tt**2)

            if sq <= 0:
                continue

            if not sq.is_prime() or sq % 4 != 1:
                continue

            # Try and get a solution sq = x^2 + y^2
            try:
                xx, yy = QF.qfbsolve(sq)
                break
            except ValueError:
                continue
        else:
            raise ValueError("Could not find a suitable endomorphism.")


        i_end = lambda P: E0(-P[0], i*P[1])
        pi_end = lambda P: E0(P[0]**p, P[1]**p)
        θ = lambda P: xx*P + yy*i_end(P) + zz*pi_end(P) + tt*i_end(pi_end(P))
        
        P2_ = θ(P2)
        Q2_ = θ(Q2)
        P3_ = θ(P3)
        Q3_ = θ(Q3)

        try:
            R = Q3
            wp_ = P3_.weil_pairing(R, B)
            wp = Q3_.weil_pairing(R, B)
            ## If DLOG exists, then Q3 is guaranteed to be lin. indep with K3
            discrete_log_pari(wp_, wp, B)

            K3_dual = θ(Q3)
        except TypeError:
            R = P3
            K3_dual = θ(P3) # P3 is guaranteed to be lin. indep with K3


        K3_dual = _E0(K3_dual[0])
        phi3 = KummerLineIsogeny(_E0, K3_dual, B)

        xP2_ = _E0(P2_[0])
        xQ2_ = _E0(Q2_[0])
        xPQ2_ = _E0((P2_ - Q2_)[0])

        ximP2_ = phi3(xP2_)
        ximQ2_ = phi3(xQ2_)
        ximPQ2_ = phi3(xPQ2_)

        P2_ = ximP2_.curve_point()
        Q2_ = ximQ2_.curve_point()

        if (P2_ - Q2_)[0] != ximPQ2_.x():
            Q2_ = -Q2_

        inverse = inverse_mod(B, 4*A)
        P2_ = inverse * P2_
        Q2_ = inverse * Q2_
        
        # We compute the (2^a,2^a)-isogeny to map points through
        P, Q = CouplePoint(-deg * P2, P2_), CouplePoint(-deg * Q2, Q2_)
        kernel = (P, Q)

        Phi = EllipticProductIsogeny(kernel, a)

        P23x = PA + PB + X0
        Q23x = QA + QB + Y0

        imP23x, imQ23x = eval_dimtwo_isog(Phi, A-deg, P23x, Q23x, 4*A*B*x)

        X_A = inverse_mod(4*A*B, x) * (4*A*B * imP23x)
        Y_A = inverse_mod(4*A*B, x) * (4*A*B * imQ23x)

        P2_og = alpha * inverse_mod(B*x, 4*A) * (B*x * imP23x)
        Q2_og =  beta * inverse_mod(B*x, 4*A) * (B*x * imQ23x)

        P3_og = gamma * inverse_mod(4*A*x, B) * (4*A*x * imP23x)
        Q3_og = gamma * inverse_mod(4*A*x, B) * (4*A*x * imQ23x)

        _, xP2, xQ2, xPQ2 = point_to_xonly(P2_og, Q2_og)
        _, xP3, xQ3, xPQ3 = point_to_xonly(P3_og, Q3_og)


        return (deg, alpha, beta, delta), (xP2, xQ2, xPQ2, xP3, xQ3, xPQ3, delta * X_A, delta * Y_A)



    def encrypt(pkA, m):
        xP2, xQ2, xPQ2, xP3, xQ3, xPQ3, X_A, Y_A = pkA
        beta = 0 #random_unit(B) # determines isogeny kernel
        d1, d2, d3, d4 = random_matrix(C) # determines X mask
        omega = random_unit(A)
        omega_inv = inverse_mod(omega, A)

        _KB = _QB.ladder_3_pt(_PB, _PQB, beta)
        phiB = KummerLineIsogeny(_E0, _KB, B)

        EB = phiB.codomain()
        xP2_B = phiB(xPA)
        xQ2_B = phiB(xQA)
        X_B = phiB(xX_0).curve_point()
        Y_B = phiB(xY_0).curve_point()
        xXY_B = phiB(xXY_0)

        if (X_B - Y_B)[0] != xXY_B.x():
            Y_B = -Y_B

        X_B, Y_B = d1*X_B + d2*Y_B, d3*X_B + d4*Y_B
        
        P2_B, Q2_B = lift_image_to_curve(PA, QA, xP2_B, xQ2_B, 4 * A, B)
        P2_B =     omega * P2_B
        Q2_B = omega_inv * Q2_B
        
        # sanity checks everywhere

        EA = xP3.parent()
        xK = xQ3.ladder_3_pt(xP3, xPQ3, beta)
        phiB_ = KummerLineIsogeny(EA, xK, B)

        EAB = phiB_.codomain().curve()
        xP2_AB = phiB_(xP2)
        xQ2_AB = phiB_(xQ2)
        xPQ2_AB = phiB_(xPQ2)


        P2_AB = xP2_AB.curve_point()
        Q2_AB = xQ2_AB.curve_point()

        if (P2_AB - Q2_AB)[0] != xPQ2_AB.x():
            Q2_AB = -Q2_AB

        P2_AB *= omega
        Q2_AB *= omega_inv

        xX_AB = phiB_(EA(X_A[0]))
        xY_AB = phiB_(EA(Y_A[0]))
        xXY_AB = phiB_(EA((X_A - Y_A)[0]))

        X_AB = xX_AB.curve_point()
        Y_AB = xY_AB.curve_point()

        if (X_AB - Y_AB)[0] != xXY_AB.x():
            Y_AB = -Y_AB

        X_AB, Y_AB = d1*X_AB + d2*Y_AB, d3*X_AB + d4*Y_AB

        # An XOF is initiated from the x-coordinates of two elliptic curve points
        # and a message is encrypted by sampling len(m) bytes from the XOF and
        # using this as a XOR key
        xof = xof_kdf(X_AB[0], Y_AB[0])
        ct = xof_encrypt(xof, m)
            
        return EB.curve(), P2_B, Q2_B, X_B, Y_B, EAB, P2_AB, Q2_AB, ct



    def decrypt(skA, ct):
        deg, alpha, beta, delta = skA
        EB, P2_B, Q2_B, X_B, Y_B, EAB, P2_AB, Q2_AB, ct_bytes = ct
        
        P2_AB = inverse_mod(alpha, A) * P2_AB
        Q2_AB =  inverse_mod(beta, A) * Q2_AB


        UAB, VAB, wp = torsion_basis_with_pairing(EAB, x)
        
        P, Q = CouplePoint(-deg * P2_B, P2_AB), CouplePoint(-deg * Q2_B, Q2_AB)
        kernel = (P, Q)
        Phi = EllipticProductIsogeny(kernel, a)


        X_AB, Y_AB = eval_dimtwo_isog(Phi, A-deg, X_B, Y_B, x, RS=(UAB, VAB))

        X_AB *= delta
        Y_AB *= delta

        # An XOF is initiated from the x-coordinates of two elliptic curve points
        # and a message is recovered by sampling a key from the XOF and computing
        # the XOR with the ciphertext
        xof = xof_kdf(X_AB[0], Y_AB[0])
        pt = xof_encrypt(xof, ct_bytes)

        return pt

    N = 100 # number of iterations
    tt = [0, 0, 0]
    correct = 0

    for _ in range(N):
        # Sample messages of length [1, 128]
        k = random.randint(1, 128)
        m = os.urandom(k)
        
        t0 = time.time_ns()
        sk, pk = keygenA()
        tt[0] += time.time_ns() - t0

        t0 = time.time_ns()
        ct = encrypt(pk, m)
        tt[1] += time.time_ns() - t0

        t0 = time.time_ns()
        m_ = decrypt(sk, ct)
        tt[2] += time.time_ns() - t0


        assert m == m_

    tt = [float(t) / N / 10^6 for t in tt]

    if len(tt) > 50: 
        tt = tt[20:]
        N = N - 20

    print("\n"*2)
    print("="*60)
    print(" "*12, "Benchmarking %d iterations (λ = %d)" % (N, lambda_), " "*12)
    print(f"KeyGen: {(tt[0]):.1f} ms")
    print(f"Encrypt: {(tt[1]):.1f} ms")
    print(f"Decrypt: {(tt[2]):.1f} ms")
    print("="*60)
    print("\n"*2)
