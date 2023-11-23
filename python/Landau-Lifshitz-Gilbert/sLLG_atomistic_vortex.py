# Based on "70 lines of Numpy" (Copyright (C) 2014 Claas Abert)
# Solving the sLLG equation : (1 + l^2)ds/dt = w x s + l * s x (w x s)

#!/usr/bin/python


########################
# Importing librairies #
########################
import numpy as np
from math import asinh, atan, sqrt, pi, sin, cos, floor, fmod
import random
import pprint as pp

hbar_ev = 6.582119514e-7

#######################################
# Setup mesh and material constraints #
#######################################
n     = (10, 10, 1)

def h_eff(s, l):
  h_ex_value = np.zeros(n + (3,))
  J = 10 * 1e-3
  
  for a in range(n[0]):
    for b in range(n[1]):
      for i in range(len(l[a][b][0])):
        h_ex_value[a][b][0] += s[l[a][b][0][i][0]][l[a][b][0][i][1]][l[a][b][0][i][2]]

  return 0.5 * 1e9 * J * h_ex_value / hbar_ev
  

def h_dmi(s, l):
  h_dmi_value = np.zeros(n + (3,))
  k = 8 * 1e-3
  
  for a in range(n[0]):
    for b in range(n[1]):
      for i in range(len(l[a][b][0])):
        aprime = l[a][b][0][i][0]
        bprime = l[a][b][0][i][1]
        r_ij = np.array([a - aprime, b - bprime, 0])
        r_ij_norm = sqrt(r_ij[0]**2 + r_ij[1]**2 + r_ij[2]**2)
        r_ij /= r_ij_norm
        h_dmi_value[a][b][0] += np.cross(s[aprime][bprime][0], r_ij)

  return 0.5 * 1e9 * k * h_dmi_value / hbar_ev

def h_dd(s):
  h_dd_value = np.zeros(n + (3,))
  D = 0.1 * 1e-3

  for a in range(n[0]):
    for b in range(n[1]):
      for aprime in range(n[0]):
        for bprime in range(n[1]):
          if(a != aprime and b != bprime):
            r_ij = np.array([a - aprime, b - bprime, 0])
            r_ij_norm = sqrt(r_ij[0]**2 + r_ij[1]**2 + r_ij[2]**2)
            r_ij /= r_ij_norm

            if r_ij_norm < 10:
              h_dd_value[a][b][0] += (3 * r_ij * np.dot(s[aprime][bprime][0], r_ij) - s[aprime][bprime][0]) / r_ij_norm**3

  return -1e9 * D * h_dd_value / hbar_ev

def PC(a, b):
  if a == n[0]: a = 0
  if a == -1: a = n[0]-1
  if b == n[1]: b = 0
  if b == -1: b = n[1]-1

  return [a,b]

def list_neighbors(condition, s):
  neighbors = [[ [0 for col in range(n[2])] for col in range(n[1])] for row in range(n[0])]

  if condition == "periodical":
    for a in range(n[0]):
      for b in range(n[1]):
        current_neighbors = []

        current_neighbors.append([PC(a-1,b)[0], PC(a-1,b)[1], 0])
        current_neighbors.append([PC(a+1,b)[0], PC(a+1,b)[1], 0])
        current_neighbors.append([PC(a,b-1)[0], PC(a,b-1)[1], 0])
        current_neighbors.append([PC(a,b+1)[0], PC(a,b+1)[1], 0])

        neighbors[a][b][0] = current_neighbors

  if condition == "strict":
    for a in range(n[0]):
      for b in range(n[1]):
        current_neighbors = []
        for i in range(-1, 2, 1):
          for j in range(-1, 2, 1):

            aprime = a + i
            bprime = b + j

            distance = sqrt((a - aprime)**2 + (b - bprime)**2)

            if distance < 1.1 and aprime >= 0 and aprime <= n[0]-1 and bprime >= 0 and bprime <= n[1]-1 and [a,b] != [aprime, bprime]:
              current_neighbors.append([aprime, bprime, 0])
        neighbors[a][b][0] = current_neighbors


  return neighbors

def h_ani(s):
  h_ani_value = np.zeros(n + (3,))
  d = np.array([0.0, 0.0, 1.0])

  K = 0.0 * 1e-3

  for a in range(n[0]):
    for b in range(n[1]):
      h_ani_value[a][b][0] = np.dot(d, s[a][b][0])

  return -0.5 * K * 1e9 * h_ani_value / hbar_ev

############################################################
# Computing the magnetization dynamics using an integrator #
############################################################

def F(m, h, gamma, alpha):
    value = np.zeros(3)
    for i in range(3):
        for k in range(3):
            value[i] += alpha * m[k] * (h[i] * m[k] - h[k] * m[i]) / (1 + alpha**2)
            for j in range(3):
                value[i] += epsilon[i][j][k] * m[k] * h[j] / (1 + alpha**2)
    return value

def RK4(m, new_m, min, max, gamma, step, l):
    t = 0
    dt = (max - min) / step
    counter = 0
    alpha = 1.0

    with open('skyrmions_10_10.dat', 'w') as f:
        while t < max:
            m = np.copy(new_m)
            h = h_eff(m, l) + h_dmi(m, l) + h_dd(m) + h_zee # Computing effective field depending on the magnetization m

            #if fmod(counter, 100) == 0: 
            m_vec = map(lambda i: np.mean(m[:,:,:,i]), range(3))
            m_norm = sqrt(m_vec[0]**2 + m_vec[1]**2 + m_vec[2]**2)

            print "%f %f %f %f %f" % (t/max, m_vec[0], m_vec[1], m_vec[2], m_norm)

            #print h
            #if fmod(counter, 100) == 0: 
            f.write("%i\n" % (n[0]*n[1]))
            f.write("Time=%f\n" % (t/max))
          
            for a in range(n[0]):
                for b in range(n[1]):
                      f.write("%i\t%i\t%i\t%.6f\t%.6f\t%.6f\n" % (a, b, 0, m[a][b][0][0], m[a][b][0][1], m[a][b][0][2]))

            v = np.zeros(3,)
            for a in range(n[0]):
              for b in range(n[1]):
                v += F(m[a][b][0], h[a][b][0], gamma, alpha)

            if sqrt(v[0]**2 + v[1]**2 + v[2]**2) < 1e10 and t > 0:
              break

            for a in range(n[0]):
                for b in range(n[1]):
                    '''
                    if a == n[0]-1 or b == n[1]-1 : alpha = 1.0
                    else: alpha = 0.02
                    '''

                    k1 = dt*F(m[a][b][0], h[a][b][0], gamma, alpha)
                    k2 = dt*F(m[a][b][0] + k1 / 2, h[a][b][0], gamma, alpha)
                    k3 = dt*F(m[a][b][0] + k2 / 2, h[a][b][0], gamma, alpha)
                    k4 = dt*F(m[a][b][0] + k3, h[a][b][0], gamma, alpha)

                    new_m[a][b][0] += (k1 + 2*k2 + 2*k3 + k4) / 6      

            t += dt
            counter += 1
    return m

def sLLG(m, m_new, gamma, step, l):
    min = 0
    max = 1e-9

    m = RK4(m, new_m, min, max, gamma, step, l)

    return m

def create_random_configuration():
    for a in range(n[0]):
      for b in range(n[1]):
        for c in range(n[2]):
          rd_nb = random.uniform(0, 2*pi)
          rd_nb_2 = random.uniform(0, 2*pi)

          m[a][b][c][0] = sin(rd_nb) * cos(rd_nb_2)
          m[a][b][c][1] = sin(rd_nb) * sin(rd_nb_2)
          m[a][b][c][2] = cos(rd_nb)

          m[a][b][c] /= sqrt(m[a][b][c][0]**2 + m[a][b][c][1]**2 + m[a][b][c][2]**2)

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

# Initial value of m
m = np.zeros(n + (3,))
m[:,:,:,2] = -1
#m[(0,1,2,3,4),:,:,1] = -1.0
#m[(5,6,7,8,9),:,:,1] = 1.0

#create_random_configuration()

#create_vortex(floor(n[0]/2), floor(n[1]/2))

new_m = np.copy(m)

l = list_neighbors("strict", m)
#pp = pp.PrettyPrinter(depth=6)
#pp.pprint(l)

# Switching
print "Setting up Zeeman field"
step = 1e5
h_zee = np.tile([0.0, 0.0, 176e9 * 25], np.prod(n)).reshape(m.shape)
m = sLLG(m, new_m, 0.0, step, l)