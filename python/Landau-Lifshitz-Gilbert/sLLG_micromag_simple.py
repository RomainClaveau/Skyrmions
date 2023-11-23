# Based on "70 lines of Numpy" (Copyright (C) 2014 Claas Abert)
# Solving the sLLG equation : (1 + l^2)ds/dt = w x s + l * s x (w x s)


########################
# Importing librairies #
########################
import numpy as np
from math import asinh, atan, sqrt, pi, fmod, log10

#######################################
# Setup mesh and material constraints #
#######################################
n     = (100, 25, 1)
dx    = (5e-9, 5e-9, 3e-9)
gamma = 2.211e5
eps   = 1e-18 # Small number to avoid division by 0
mu0   = 4e-7 * pi
ms    = 8e5
A     = 1.3e-11

#############################
# Defining Newell functions #
#############################
def f(p):
  x, y, z = abs(p[0]), abs(p[1]), abs(p[2])
  return + y / 2.0 * (z**2 - x**2) * asinh(y / (sqrt(x**2 + z**2) + eps)) \
         + z / 2.0 * (y**2 - x**2) * asinh(z / (sqrt(x**2 + y**2) + eps)) \
         - x*y*z * atan(y*z / (x * sqrt(x**2 + y**2 + z**2) + eps))       \
         + 1.0 / 6.0 * (2*x**2 - y**2 - z**2) * sqrt(x**2 + y**2 + z**2)

def g(p):
  x, y, z = p[0], p[1], abs(p[2])
  return + x*y*z * asinh(z / (sqrt(x**2 + y**2) + eps))                         \
         + y / 6.0 * (3.0 * z**2 - y**2) * asinh(x / (sqrt(y**2 + z**2) + eps)) \
         + x / 6.0 * (3.0 * z**2 - x**2) * asinh(y / (sqrt(x**2 + z**2) + eps)) \
         - z**3 / 6.0 * atan(x*y / (z * sqrt(x**2 + y**2 + z**2) + eps))        \
         - z * y**2 / 2.0 * atan(x*z / (y * sqrt(x**2 + y**2 + z**2) + eps))    \
         - z * x**2 / 2.0 * atan(y*z / (x * sqrt(x**2 + y**2 + z**2) + eps))    \
         - x*y * sqrt(x**2 + y**2 + z**2) / 3.0

################################
# Setup demagnetization tensor #
################################

def set_n_demag(c, permute, func):
  print "Setup magnetization tensor"
  it = np.nditer(n_demag[:,:,:,c], flags=['multi_index'], op_flags=['writeonly'])
  while not it.finished:
    value = 0.0
    for i in np.rollaxis(np.indices((2,)*6), 0, 7).reshape(64, 6):
      idx = map(lambda k: (it.multi_index[k] + n[k]) % (2*n[k]) - n[k], range(3))
      value += (-1)**sum(i) * func(map(lambda j: (idx[j] + i[j] - i[j+3]) * dx[j], permute))
    it[0] = - value / (4 * pi * np.prod(dx))
    it.iternext()

n_demag = np.zeros([1 if i==1 else 2*i for i in n] + [6])
for i, t in enumerate(((f,0,1,2),(g,0,1,2),(g,0,2,1),(f,1,2,0),(g,1,2,0),(f,2,0,1))):
  set_n_demag(i, t[1:], t[0])

m_pad     = np.zeros([1 if i==1 else 2*i for i in n] + [3])
f_n_demag = np.fft.rfftn(n_demag, axes = filter(lambda i: n[i] > 1, range(3)))

######################################################################
# Computing effective field (Demagnetization field + Exchange field) #
######################################################################
def h_eff(m):
  # demag field
  m_pad[:n[0],:n[1],:n[2],:] = m
  f_m_pad = np.fft.rfftn(m_pad, axes = filter(lambda i: n[i] > 1, range(3)))
  f_h_demag_pad = np.zeros(f_m_pad.shape, dtype=f_m_pad.dtype)
  f_h_demag_pad[:,:,:,0] = (f_n_demag[:,:,:,(0, 1, 2)]*f_m_pad).sum(axis = 3)
  f_h_demag_pad[:,:,:,1] = (f_n_demag[:,:,:,(1, 3, 4)]*f_m_pad).sum(axis = 3)
  f_h_demag_pad[:,:,:,2] = (f_n_demag[:,:,:,(2, 4, 5)]*f_m_pad).sum(axis = 3)
  h_demag = np.fft.irfftn(f_h_demag_pad, axes = filter(lambda i: n[i] > 1, range(3)))[:n[0],:n[1],:n[2],:]

  # exchange field
  h_ex = - 2 * m * sum([1/x**2 for x in dx])
  for i in range(6):
    h_ex += np.repeat(m, 1 if n[i%3] == 1 else [i/3*2] + [1]*(n[i%3]-2) + [2-i/3*2], axis = i%3) / dx[i%3]**2

  h = ms * h_demag + 2 * A / (mu0 * ms) * h_ex

  return h

############################################################
# Computing the magnetization dynamics using an integrator #
############################################################

def F(m, h, gamma, alpha):
    return gamma / (1 + alpha**2) * np.cross(h, m) + gamma * alpha / (1 + alpha**2) * np.cross(m, np.cross(h, m))

def RK4(m, min, max, gamma, alpha, step, write):
    t = 0
    dt = (max - min) / step
    c = 0

    if write == "yes":
        with open('dynamics.dat', 'w') as f:
            while t < max:
                h = h_eff(m) + h_zee # Computing effective field depending on the magnetization m

                k1 = dt*F(m, h, gamma, alpha)
                k2 = dt*F(m + k1 / 2, h, gamma, alpha)
                k3 = dt*F(m + k2 / 2, h, gamma, alpha)
                k4 = dt*F(m + k3, h, gamma, alpha)
                m += (k1 + 2*k2 + 2*k3 + k4) / 6

                m /= np.repeat(np.sqrt((m*m).sum(axis=3)), 3).reshape(m.shape)

                if fmod(c, 100) == 0: # Printing 1 point for 100 points computed
                    m_vec = map(lambda i: np.mean(m[:,:,:,i]), range(3))
                    print "%.15f %f %f %f" % (t/max, m_vec[0], m_vec[1], m_vec[2])
                    f.write("%.15f %f %f %f\n" % (t/max, m_vec[0], m_vec[1], m_vec[2]))

                t += dt
                c += 1
    else:
        while t < max:
            h = h_eff(m) + h_zee # Computing effective field depending on the magnetization m

            k1 = dt*F(m, h, gamma, alpha)
            k2 = dt*F(m + k1 / 2, h, gamma, alpha)
            k3 = dt*F(m + k2 / 2, h, gamma, alpha)
            k4 = dt*F(m + k3, h, gamma, alpha)
            m += (k1 + 2*k2 + 2*k3 + k4) / 6

            m /= np.repeat(np.sqrt((m*m).sum(axis=3)), 3).reshape(m.shape)

            m_vec = map(lambda i: np.mean(m[:,:,:,i]), range(3))
            print "%.15f %f %f %f" % (t/max, m_vec[0], m_vec[1], m_vec[2])

            t += dt


    return m

def sLLG(m, gamma, alpha, step, write):
    min = 0
    max = 1e-9

    m = RK4(m, min, max, gamma, alpha, step, write)

    return m

# Initial value of m
m = np.zeros(n + (3,))
m[1:-1,:,:,0]   = 1.0
m[(-1,0),:,:,1] = 1.0

# Settings magnetization in s-state
print "Computing s-state"
step = 5000
h_zee = np.zeros(n + (3,))
alpha = 1
m = sLLG(m, gamma, alpha, step, "no")

# Switching
print "Setting up Zeeman field"
step = 50000
h_zee = np.tile([-24.6e-3/mu0, +4.3e-3/mu0, 0.0], np.prod(n)).reshape(m.shape)
alpha = 0.02
m = sLLG(m, gamma, alpha, step, "yes")