#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys 
import os
import numpy as np
from math import asinh, atan, sqrt, pi, sin, cos, floor, fmod, acos
import random
from scipy.spatial import Delaunay
#import numba

hbar_ev = 6.582119514 * 1e-7
hbar = 1.05457148 * 1e-34
kboltzmann = 1.380649 * 1e-23

sys.path.append(os.path.dirname(os.path.realpath("create_configuration.py")))
sys.path.append(os.path.dirname(os.path.realpath("load_configuration.py")))

import create_configuration
import load_configuration

# Creating configuration
filename_config = "configuration_15_15_1_strict_random.json"
filename_output = "influence_temperature_skyrmions_15_15_1.dat"
filename_output_charge = "charge_skyrmions_15_15_1.dat"

#create_configuration.create_configuration(500, 500, 1, "strict", "random", filename_config)
data = load_configuration.load_configuration(filename_config)

n     = (data[0], data[1], data[2])
m = np.asarray(data[3])
W = np.asarray(data[4])
S = np.asarray(data[5])

def compute_Q_2(s):
  Q = 0.0

  for a in range(n[0]):
    for b in range(n[1]):
      dx = 0.0
      dy = 0.0

      if a == 0:
        dx = (s[a+1][b][0] - s[a][b][c])/(1.0)
      elif a > 0 and a < n[0]-1:
        dx = (s[a+1][b][0] - s[a-1][b][0])/(2.0)
      else:
        dx = (s[a-1][b][0] - s[a][b][c])/(1.0)

      if b == 0:
        dy = (s[a][b+1][0] - s[a][b][c])/(1.0)
      elif b > 0 and b < n[1]-1:
        dy = (s[a][b+1][0] - s[a][b-1][0])/(2.0)
      else:
        dy = (s[a][b-1][0] - s[a][b][c])/(1.0)

      Q += np.dot(s[a][b][c], np.cross(dx, dy))

  return -Q/(4 * pi)/0.667

# Computing topological charge
def compute_Q(s):

  # Delaunay triangulation
  points = []
  A = 0.0

  for a in range(n[0]):
    for b in range(n[1]):
      points.append([a + 1e-15 * random.random(), b + 1e-15 * random.random()])

  points = np.array(points)

  tri = Delaunay(points).simplices

  for a in tri:
    counter_clockwise = False
    counter = 0
    permutations = [
      [a[0], a[1], a[2]],
      [a[0], a[2], a[1]],
      [a[1], a[0], a[2]],
      [a[1], a[2], a[0]],
      [a[2], a[0], a[1]],
      [a[2], a[1], a[0]]
    ]

    pt = permutations[counter]

    while(counter_clockwise == False and counter < 5):
      rij = np.array([points[pt[1]][0] - points[pt[0]][0],points[pt[1]][1] - points[pt[0]][1], 0])

      rik = np.array([points[pt[2]][0] - points[pt[0]][0],points[pt[2]][1] - points[pt[0]][1], 0])

      r_normal = np.array([0, 0, 1])

      if np.sign(np.dot(r_normal, np.cross(rij, rik))) == -1:
        counter_clockwise = True
        break

      counter += 1
      pt = permutations[counter]
    
    ni = s[points[pt[0]][0]][points[pt[0]][1]][0]
    nj = s[points[pt[1]][0]][points[pt[1]][1]][0]
    nk = s[points[pt[2]][0]][points[pt[2]][1]][0]

    A_tmp = ((1 + np.dot(ni, nj) + np.dot(ni, nk) + np.dot(nj, nk)) / sqrt(2 * (1 + np.dot(ni, nj)) * (1 + np.dot(ni, nk)) * (1 + np.dot(nj, nk))))

    #A_tmp = ((1 + np.dot(ni, nj) + np.dot(ni, nk) + np.dot(nj, nk)) / sqrt(2 * (1 + np.dot(ni, nj)) * (1 + np.linalg.norm(ni) * np.linalg.norm(nj)) * (1 + np.linalg.norm(ni) * np.linalg.norm(nk))))

    if A_tmp > 1: A_tmp = 1
    if A_tmp < -1: A_tmp = -1

    A_tmp = 2 * np.sign(np.dot(ni, np.cross(nj, nk))) * acos(A_tmp)

    A += A_tmp

  return A / (4 * pi)

def h_eff(s, l):
        h_ex_value = np.zeros(n + (3,))
        J = 10 * 1e-3
        
        for a in range(n[0]):
            for b in range(n[1]):
                for c in range(n[2]):
                    for i in range(len(l[a][b][c])):
                        h_ex_value[a][b][c] += s[l[a][b][c][i][0]][l[a][b][c][i][1]][l[a][b][c][i][2]]

        return 0.5 * 1e9 * J * h_ex_value / hbar_ev
  
    # Computing Dzyaloshinskii-Moriya field
def h_dmi(s, l):
    h_dmi_value = np.zeros(n + (3,))
    k = 8 * 1e-3
    
    for a in range(n[0]):
        for b in range(n[1]):
            for c in range(n[2]):
                for i in range(len(l[a][b][c])):
                    aprime = l[a][b][c][i][0]
                    bprime = l[a][b][c][i][1]
                    cprime = l[a][b][c][i][2]
                    r_ij = np.array([a - aprime, b - bprime, c - cprime])
                    r_ij_norm = sqrt(r_ij[0]**2 + r_ij[1]**2 + r_ij[2]**2)
                    r_ij = r_ij / r_ij_norm
                    h_dmi_value[a][b][c] += np.cross(s[aprime][bprime][cprime], r_ij)

    return 0.5 * 1e9 * k * h_dmi_value / hbar_ev

def h_dd(s):
  h_dd_value = np.zeros(n + (3,))
  D = 0.5 * 1e-3

  for a in range(n[0]):
    for b in range(n[1]):
      for aprime in range(n[0]):
        for bprime in range(n[1]):
          if(a != aprime and b != bprime):
            r_ij = np.array([a - aprime, b - bprime, 0])
            r_ij_norm = sqrt(r_ij[0]**2 + r_ij[1]**2 + r_ij[2]**2)
            r_ij = r_ij / r_ij_norm

            if r_ij_norm < 2:
              h_dd_value[a][b][c] += (3 * r_ij * np.dot(s[aprime][bprime][0], r_ij) - s[aprime][bprime][0]) / r_ij_norm**3

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

        neighbors[a][b][c] = current_neighbors

  if condition == "strict":
    for a in range(n[0]):
      for b in range(n[1]):
        for c in range(n[2]):
          current_neighbors = []
          for i in range(-1, 2, 1):
            for j in range(-1, 2, 1):
              for k in range(-1, 2, 1):

                aprime = a + i
                bprime = b + j
                cprime = c + k

                distance = sqrt((a - aprime)**2 + (b - bprime)**2 + + (c - cprime)**2)

                if distance < 1.1 and aprime >= 0 and aprime <= n[0]-1 and bprime >= 0 and bprime <= n[1]-1 and cprime >= 0 and cprime <= n[2]-1 and [a,b,c] != [aprime, bprime, cprime]:
                  current_neighbors.append([aprime, bprime, cprime])
          neighbors[a][b][c] = current_neighbors


  return neighbors

def h_ani(s):
  h_ani_value = np.zeros(n + (3,))
  d = np.array([0.0, 0.0, 1.0])

  K = 0.0 * 1e-3

  for a in range(n[0]):
    for b in range(n[1]):
      h_ani_value[a][b][c] = np.dot(d, s[a][b][c])

  return -0.5 * K * 1e9 * h_ani_value / hbar_ev

############################################################
# Computing the magnetization dynamics using an integrator #
############################################################

def F(s, h, W, S, tau):
    value = np.zeros(3)
    for i in range(3):
        for k in range(3):
            value[i] += coeff * alpha * (h[i] * S[k][k] - h[k] * S[k][i])
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

def RK4(m, new_m, min, max, gamma, step, l):
    global D

    t = 0
    dt = (max - min) / step
    counter = 0
    E = 1e-5
    
    # Printing initial configuration
    m_vec = map(lambda i: np.mean(m[:,:,:,i]), range(3))
    m_norm = sqrt(m_vec[0]**2 + m_vec[1]**2 + m_vec[2]**2)
    with open(filename_output_charge, 'w') as f2, open(filename_output, 'w') as f:
        #Q = compute_Q_2(m)

        #f2.write(("%f %f %.2f\n") % (t/max, D, Q))

        while t < max:
            m = np.copy(new_m)
            W = np.copy(new_W)
            S = np.copy(new_S)

            #Q = compute_Q_2(m)

            h = h_eff(m, l) + h_dmi(m, l) + h_zee # Computing effective field depending on the magnetization m

            m_vec = map(lambda i: np.mean(m[:,:,:,i]), range(3))
            m_norm = sqrt(m_vec[0]**2 + m_vec[1]**2 + m_vec[2]**2)

            # Printing
            print "%f %f %f %f %f %f" % (t/max, m_vec[0], m_vec[1], m_vec[2], m_norm, D)
            #print "%f %f %.2f" % (t/max, D, Q)
            
            
            f.write("%i\n" % (n[0]*n[1]*n[2]))
            f.write("Time=%f\n" % (t/max))
            

            #f2.write(("%f %f %.2f\n") % (t/max, D, Q))

            
            for a in range(n[0]):
                for b in range(n[1]):
                  for c in range(n[2]):
                    f.write("%i\t%i\t%i\t%f\t%f\t%f\t%f\n" % (a, b, c, m[a][b][c][0], m[a][b][c][1], m[a][b][c][2], D))
            

            for a in range(n[0]):
                D = 1e15 * a / (n[0]-1)
                for b in range(n[1]):
                    for c in range(n[2]):
                      a1 = dt * F(m[a][b][c], h[a][b][c], W[a][b][c], S[a][b][c], tau)
                      b1 = dt * G(m[a][b][c], h[a][b][c], W[a][b][c], S[a][b][c], tau)
                      c1 = dt * H(m[a][b][c], h[a][b][c], W[a][b][c], S[a][b][c], tau)

                      a2 = dt * F(m[a][b][c] + a1/2, h[a][b][c], W[a][b][c] + b1/2, S[a][b][c] + c1/2, tau)
                      b2 = dt * G(m[a][b][c] + a1/2, h[a][b][c], W[a][b][c] + b1/2, S[a][b][c] + c1/2, tau)
                      c2 = dt * H(m[a][b][c] + a1/2, h[a][b][c], W[a][b][c] + b1/2, S[a][b][c] + c1/2, tau)

                      a3 = dt * F(m[a][b][c] + a2/2, h[a][b][c], W[a][b][c] + b2/2, S[a][b][c] + c2/2, tau)
                      b3 = dt * G(m[a][b][c] + a2/2, h[a][b][c], W[a][b][c] + b2/2, S[a][b][c] + c2/2, tau)
                      c3 = dt * H(m[a][b][c] + a2/2, h[a][b][c], W[a][b][c] + b2/2, S[a][b][c] + c2/2, tau)

                      a4 = dt * F(m[a][b][c] + a3, h[a][b][c], W[a][b][c] + b3, S[a][b][c] + c3, tau)
                      b4 = dt * G(m[a][b][c] + a3, h[a][b][c], W[a][b][c] + b3, S[a][b][c] + c3, tau)
                      c4 = dt * H(m[a][b][c] + a3, h[a][b][c], W[a][b][c] + b3, S[a][b][c] + c3, tau)

                      new_m[a][b][c] += (a1 + 2 * a2 + 2 * a3 + a4) / 6
                      new_W[a][b][c] += (b1 + 2 * b2 + 2 * b3 + b4) / 6
                      new_S[a][b][c] += (c1 + 2 * c2 + 2 * c3 + c4) / 6

            m_vec_new = map(lambda i: np.mean(new_m[:,:,:,i]), range(3))
            m_norm_new = sqrt(m_vec_new[0]**2 + m_vec_new[1]**2 + m_vec_new[2]**2)

            print abs(m_norm_new - m_norm)
            print counter

            if counter == 200:
                break

            h_vec = map(lambda i: np.mean(h[:,:,:,i]), range(3))
            h_norm = sqrt(h_vec[0]**2 + h_vec[1]**2 + h_vec[2]**2)

            dt = 1e-14

            t += dt
            counter += 1
    return m

def dLLB(m, m_new, gamma, step, l):
    min = 0
    max = 1e-9

    m = RK4(m, new_m, min, max, gamma, step, l)

    return m

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

# Parameters
tau = 1e-11
step = 1e5
alpha = 1.0
coeff = 1 / (1 + alpha * alpha)
D = 1e14

# Copying dynamics at t-1
new_m = np.copy(m)
new_W = np.copy(W)
new_S = np.copy(S)

l = list_neighbors("strict", m)

# Switching
E = 1e-5
print "Setting up Zeeman field"
h_zee = np.tile([0.0, 0.0, 176e9 * 0], np.prod(n)).reshape(m.shape)
m = dLLB(m, new_m, 0.0, step, l)
