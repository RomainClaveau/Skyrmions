# Based on "70 lines of Numpy" (Copyright (C) 2014 Claas Abert)
# Solving the sLLG equation : (1 + l^2)ds/dt = w x s + l * s x (w x s)


########################
# Importing librairies #
########################
import numpy as np
from math import asinh, atan, sqrt, pi

#######################################
# Setup mesh and material constraints #
#######################################
n     = (1, 1, 1)
dx    = (5e-9, 5e-9, 3e-9)
gamma = 2.211e5
#gamma = 1
eps   = 1e-25 # Small number to avoid division by 0
mu0   = 4e-7 * pi
ms    = 8e5
A     = 1.3e-11
D     = 1e9 # 1 rad.GHz
tau   = 1e-11 # Weak correlation time

############################################################
# Computing the magnetization dynamics using an integrator #
############################################################

def computing_W(): # Computing <w x s> with w the random Gaussian variable
    value = np.zeros(n + (3,3,))
    value[:,:,:,1,2] = -D
    value[:,:,:,2,1] = D
    return value


def F(s, h, W, S, tau):
    value = np.zeros(3)
    for i in range(3):
        for k in range(3):
            value[i] -= coeff * alpha * (h[k] * S[k][i] - h[i] * S[k][k])
            for j in range(3):
                value[i] += coeff * (epsilon[i][j][k] * h[j] * s[k] + epsilon[i][j][k] * W[j][k])
    return value

def G(s, h, W, S, tau):
    value = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            value[i][j] -= W[i][j] / tau
            for l in range(3):
                value[i][j] += coeff * D * epsilon[j][i][l] * s[l] / tau - coeff * alpha * (h[l] * (s[l] * W[i][j] + s[j] * W[i][l]) - 2 * h[j] * s[l] * W[i][l])
                for k in range(3):
                    value[i][j] += coeff * (epsilon[j][k][l] * h[k] * W[i][l])
    return value

def H(s, h, W, S, tau):
    value = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            for l in range(3):
                value[i][j] -= coeff * alpha * (h[l] * (s[i] * S[l][j] + s[l] * S[i][j] + s[j] * S[i][l] - 2 * s[i] * s[j] * s[l]) - h[j] * (s[i] * S[l][l] + 2 * s[l] * S[i][l] - 2 * s[i] * s[l] * s[l]) + h[l] * (s[j] * S[l][i] + s[l] * S[j][i] + s[i] * S[j][l] - 2 * s[j] * s[i] * s[l]) - h[i] * (s[j] * S[l][l] + 2 * s[l] * S[j][l] - 2 * s[j] * s[l] * s[l]))
                for k in range(3):
                    value[i][j] += coeff * (epsilon[j][k][l] * h[k] * S[i][l] + epsilon[j][k][l] * (s[i] * W[k][l] + s[l] * W[k][i]) + epsilon[i][k][l] * h[k] * S[j][l] + epsilon[i][k][l] * (s[j] * W[k][l] + s[l] * W[k][j]))
    return value

def RK4(min, max, s, S, W, gamma, alpha, step):
    t = 0
    dt = (max - min) / step
    c = 0
    Q = 1e-1

    # Computing W
    W = computing_W()

    with open('dllb_zeeman_rk4_s.dat', 'w') as f:
        while t < max:
            h = h_zee # Computing effective field depending on the magnetization s

            vec_s = map(lambda i: np.mean(s[:,:,:,i]), range(3))

            s_sqr = sqrt(vec_s[0]**2 + vec_s[1]**2 + vec_s[2]**2)
            
            print "%.8f\t%f\t%f\t%f\t%f" % (t*1e9, vec_s[0], vec_s[1], vec_s[2], s_sqr)
            f.write("%.8f\t%f\t%f\t%f\t%f\n" % (t*1e9, vec_s[0], vec_s[1], vec_s[2], s_sqr))

            for a in range(n[0]):
                for b in range(n[1]):
                    for c in range(n[2]):
                        a1 = dt * F(s[a][b][c], h[a][b][c], W[a][b][c], S[a][b][c], tau)
                        b1 = dt * G(s[a][b][c], h[a][b][c], W[a][b][c], S[a][b][c], tau)
                        c1 = dt * H(s[a][b][c], h[a][b][c], W[a][b][c], S[a][b][c], tau)

                        a2 = dt * F(s[a][b][c] + a1/2, h[a][b][c], W[a][b][c] + b1/2, S[a][b][c] + c1/2, tau)
                        b2 = dt * G(s[a][b][c] + a1/2, h[a][b][c], W[a][b][c] + b1/2, S[a][b][c] + c1/2, tau)
                        c2 = dt * H(s[a][b][c] + a1/2, h[a][b][c], W[a][b][c] + b1/2, S[a][b][c] + c1/2, tau)

                        a3 = dt * F(s[a][b][c] + a2/2, h[a][b][c], W[a][b][c] + b2/2, S[a][b][c] + c2/2, tau)
                        b3 = dt * G(s[a][b][c] + a2/2, h[a][b][c], W[a][b][c] + b2/2, S[a][b][c] + c2/2, tau)
                        c3 = dt * H(s[a][b][c] + a2/2, h[a][b][c], W[a][b][c] + b2/2, S[a][b][c] + c2/2, tau)

                        a4 = dt * F(s[a][b][c] + a3, h[a][b][c], W[a][b][c] + b3, S[a][b][c] + c3, tau)
                        b4 = dt * G(s[a][b][c] + a3, h[a][b][c], W[a][b][c] + b3, S[a][b][c] + c3, tau)
                        c4 = dt * H(s[a][b][c] + a3, h[a][b][c], W[a][b][c] + b3, S[a][b][c] + c3, tau)

                        s += (a1 + 2 * a2 + 2 * a3 + a4) / 6
                        W += (b1 + 2 * b2 + 2 * b3 + b4) / 6
                        S += (c1 + 2 * c2 + 2 * c3 + c4) / 6

            # Variable step
            vec_h = map(lambda i: np.mean(h[:,:,:,i]), range(3))
            norm_h = sqrt(vec_h[0] * vec_h[0] + vec_h[1] * vec_h[1] + vec_h[2] * vec_h[2])

            dt = Q/norm_h

            t += dt
            c += 1

    return s

# Initial values
s = np.zeros(n + (3,))
s[:,:,:,:] = 0
S = np.zeros(n + (3, 3,))
S[:,:,:,:,:] = 0
W = np.zeros(n + (3, 3,))
W[:,:,:,:,:] = 0
s[:,:,:,0]   = 1.0 # s = (1, 0, 0)
S[:,:,:,0,0] = 1.0

# Levi-Civita symbol
epsilon = np.zeros((3,3,3))

for a in range(3):
  for b in range(3):
    for c in range(3):
      if [a,b,c] == [0, 1, 2] or [a,b,c] == [1,2,0] or [a,b,c] == [2,0,1]: 
        epsilon[a][b][c] = 1
      elif a == b or a == c or b == c: 
        epsilon[a][b][c] = 0
      else: 
        epsilon[a][b][c] = -1

# Switching
print "Setting up Zeeman field"
step = 1e3
h_zee = np.tile([0.0, 0.0, 176 * 1e9], np.prod(n)).reshape(s.shape) # h_zee = (0, 0, 176 rad.GHz)
alpha = 0.1
coeff = 1 / (1 + alpha * alpha)
s = RK4(0.0, 3e-10, s, S, W, gamma, alpha, step)