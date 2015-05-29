reset()

n = 32

#generate q such that X^n+1 splits into n polynomials modulo q
q = 2^n+1
while not q in Primes():
    q = q+2^n

Zq = Integers(q)
Rtmp.<X> = PolynomialRing(Zq);

#verify that X^n+1 splits into n polynomials of degree 1 mod q
P = [l[0] for l in (Rtmp(X^n+1)).factor()]

b = True
for l in P:
    if (l.degree() != 1 or l[1] != 1):
        b = False

if b == False:
    print "attention, P = ", P

#As all polynomials have degree one, we are only interested into their constant coefficient
P = [-l[0] for l in P]

I = Rtmp.ideal([X^n+1])
R = Rtmp.quotient_ring(I)
Y = R.gen()

#convert an element of R into an element of Rtmp
def convert_R_to_Rtmp(P):
    list_coeff = list(P)
    Q = Rtmp(0)
    Y = 1
    for l in list_coeff:
        Q = Rtmp(Q+l*Y);
        Y = X*Y
    return Q
    
#generate a random invertible matrix S
S = random_matrix(R,2)
while (gcd(convert_R_to_Rtmp(det(S)), Rtmp(X^n+1)) != 1):
    S = random_matrix(R,2)
    
#S1 is the inverse of S modulo q
S1 = matrix([[S[1,1],-S[0,1]],[-S[1,0], S[0,0]]])*det(S)^(-1)

#m is the number of public matrices (triangular matrices hidden with S)
m = 2

#M is the list of the secret triangular matrices
M = [matrix([[R.random_element(), R.random_element()],[0,R.random_element()]]) for i in range(m)]

#N is the list of matrices hidden by S
N = [S1*M[i]*S for i in range(m)]


#w and z are the bottom coefficients of the matrix Q we want to recover, and W and Z are the list of the values of w and z in each of the subfields
W = []
Z = []

for k in P:
    var('w,z')
    mat = matrix([[convert_R_to_Rtmp(N[i][0,0])(k)-convert_R_to_Rtmp(N[i][1,1])(k), convert_R_to_Rtmp(N[i][1,0])(k), -convert_R_to_Rtmp(N[i][0,1])(k)] for i in range(m)])
    kernel = list((mat.right_kernel()).basis()[0])
    w1 = sqrt(kernel[1])
    if not w1 in Zq:
        kernel = [l*kernel[1] for l in kernel]
        w1 = sqrt(kernel[1])
    z1 = sqrt(kernel[2])
    if (w1 * z1 != kernel[0]):
	    w1 = -w1
    if (w1*z1 != kernel[0]):
	    print "\n\nattention, w1*z1 != kernel[0]\n\n"
    W = W+[(k,w1)]
    Z = Z+[(k,z1)]

#reconstruction w and Z from their values in each of the subfields
wf = R(Rtmp.lagrange_polynomial(W))
zf = R(Rtmp.lagrange_polynomial(Z))
	
#generating the top coefficients of Q so that Q is invertible
Q = matrix([[R.random_element(), R.random_element()],[zf,wf]])

while(gcd(convert_R_to_Rtmp(det(Q)), Rtmp(X^n+1)) != 1):
    Q = matrix([[R.random_element(), R.random_element()],[zf,wf]])

print "Q*S^(-1) = ", Q*S1
#test if Q*S^{-1} is indeed triangular
if (Q*S1)[1][0] == 0:
    print "\nQ*S^(-1) is upper triangular"
else:
    print "\nFailure: Q*S^(-1) is not upper triangular"
