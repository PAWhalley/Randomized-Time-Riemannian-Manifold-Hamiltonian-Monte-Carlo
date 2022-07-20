#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:00:03 2022

@author: s2110992
"""

import numba
import numpy as np
from numpy.linalg import inv

import time
import random

import os
from os import sysplot

data_matrix = np.zeros((120,2000))


for k in range(1,2001):
    file_location = './cov_shrink_simulations/'
    if k < 10:
        #print(k)
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_01_01_000' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp1             =  np.append(xi_p,xi_m, axis = 0)
        
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_01_02_000' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp2             =  np.append(xi_p,xi_m, axis = 0)
        
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_02_02_000' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp3             =  np.append(xi_p,xi_m, axis = 0)   
    
    if 10 <= k < 100:
        #print(k)
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_01_01_00' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp1             =  np.append(xi_p,xi_m, axis = 0)
        
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_01_02_00' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp2              =  np.append(xi_p,xi_m, axis = 0)
        
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_02_02_00' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp3             =  np.append(xi_p,xi_m, axis = 0)   
        
    if 100 <= k < 1000:
        #print(k)
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_01_01_0' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp1             =  np.append(xi_p,xi_m, axis = 0)
        
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_01_02_0' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp2             =  np.append(xi_p,xi_m, axis = 0)
        
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_02_02_0' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp3             =  np.append(xi_p,xi_m, axis = 0)  
        
    if 1000 <= k:
        #print(k)
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_01_01_' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp1             =  np.append(xi_p,xi_m, axis = 0)
        
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_01_02_' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp2             =  np.append(xi_p,xi_m, axis = 0)
        
        file_name         =  'xi_athena_final_shrinkage_lognormal_2zbin_shapenoise_02_02_' + str(k) + '.dat'
        file              =  os.path.join(file_location,file_name)
        dat_file          =  np.genfromtxt(fname=file,  skip_header =1, dtype='unicode', invalid_raise = False)
        xi_p              =  dat_file[:,1]
        xi_m              =  dat_file[:,2]
        xi_p              =  np.reshape(xi_p,(len(xi_p),1))
        xi_m              =  np.reshape(xi_m,(len(xi_m),1))
        temp3             =  np.append(xi_p,xi_m, axis = 0)  
        
    temp12  = np.append(temp1, temp2, axis = 0)
    temp123 = np.append(temp12,temp3, axis = 0)
    #print(np.shape(temp123))
    data_matrix[:,k-1]  = np.reshape(temp123,(len(temp123),))
        
    
    
    
# Astronomical Data
#n = 2000
N = 120
N = 60

p =  N//6
d =  N
m =  p


n = N*2//3
Data = data_matrix[:N,:2000]

mean = np.mean(Data,axis = 1)
var = np.var(Data,axis = 1)

    
    
for i in range(2000):
    Data[:,i] -= mean
    Data[:,i] /= np.sqrt(var)


# sample variance
mean2 = np.zeros(d)
for i in range(80):
    mean2 += Data[:,i]/80
    
S_small = np.zeros((d,d))
for i in range(80):
    x = Data[:,i] - mean2
    S_small += np.outer(x,x)/80
    

D_2 = np.diag(np.diag(S_small))
Pmat = (S_small - D_2)


EIG = np.linalg.eig(Pmat)
A  = np.zeros((d,p))
for i in range(p):
    A[:,i] = EIG[1][:,i]



D_1 = np.diag(EIG[0][:p])

Approx = np.matmul(A,np.matmul(D_1,np.transpose(A)))
D_2 = np.diag(np.diag(S_small)) - np.diag(np.diag(Approx))
Approx += D_2



#p<=d



dimension = int(d*p+p+d)
num_of_constraints  = int(p*(p+1)/2)

σ_1 = 2
σ_2 = 2
#Matrix in distribution

@numba.jit(nopython=True)
def vec_to_matrix(q):
    X = np.zeros((d,p))
    for i in range(d):
        for j in range(p):
            X[i,j] = q[j*d+i]
    return X

@numba.jit(nopython=True)
def matrix_to_vec(X):
    #initialising filler array
    x = np.zeros(d*p)
    
    for i in range(d*p):
        i_index = i%d
        j_index =  int((i - i_index)/d)
        x[i] = X[i_index,j_index]
    return x

@numba.jit(nopython=True)
def dot_product(v1,v2):
    dot = 0
    for i in range(len(v1)):
        
        dot += v1[i]*v2[i]
        
    return dot

@numba.jit(nopython=True)
def matmul(matrix1,matrix2):
    a = matrix1.shape[0]
    b = matrix2.shape[1]
    c = matrix2.shape[0]
    rmatrix = np.zeros((a,b))
    for i in range(a):
        for j in range(b):
            for k in range(c):
                rmatrix[i,j] += matrix1[i,k] * matrix2[k,j]
    return rmatrix

@numba.jit(nopython=True)
def matrix_vec_multiplication(A,x):
    v = np.zeros(len(A))
    
    for i in range(len(A)):
        for j in range(len(x)):
                v[i] += A[i][j] * x[j]
    return v


@numba.jit(nopython=True)
def g_ij(q,i,j):
    
    q_mat = q[:d*p]
    
    X = vec_to_matrix(q_mat)
    
    if i==j:
        y = np.linalg.norm(X[:,i])**2 - 1
    else:
        y = dot_product(X[:,i],X[:,j])
        
    return y


@numba.jit(nopython=True)
def G(q): #considering i<j.
    
    q_mat = q[:d*p]
    
    X = vec_to_matrix(q_mat)
    
    z = np.zeros((dimension,num_of_constraints))
    
    for i in range(p): #block diagonals
        z[d*i:d*(i+1),int(p*i-0.5*i*(i-1)):int(p*(i+1) - 0.5*i*(i+1))] = X[:,i:]
    
        #vector diagonals
        for j in range(p-i):
            z[(j+i)*d:(j+i+1)*d,int(p*i-0.5*i*(i-1) + j)] += X[:,i]  
    z = z.T #could implement this above
    return z

@numba.jit(nopython=True)
def potential_derv_fast(q):
    #Can check with numerical differentiation.
    
    X = vec_to_matrix(q[:d*p])
    d_1 = q[d*p:d*p+p]
    d_2 = q[d*p+p:]
    
    D_1 = np.diag(d_1)
    D_2 = np.diag(d_2)
    
    Σ = matmul(matmul(X,D_1),np.transpose(X)) + D_2
    
    Σ_inv_T = np.transpose(np.linalg.inv(Σ))

    #Constructing M
    M_kl = 0.5*n*Σ_inv_T 
    
    for k in range(d):
        for l in range(d):
            for r in range(n):
                M_kl[k,l] -= 0.5*dot_product(Data[:,r]-mean,Σ_inv_T[k,:])*dot_product(Σ_inv_T[:,l],Data[:,r]-mean)
    
    #dUdX
    dUdX = np.zeros((d,p))
    
    for i in range(d):
        
        for j in range(p):
            
            for k in range(d):
                
                for l in range(d):
                    
                    if k == i and l == i:
                        
                        dΣ_kl_dX_ij = 2*X[i,j]*d_1[j]
                        
                    elif k == i:
                        
                        dΣ_kl_dX_ij = d_1[j]*X[l,j]
                        
                    elif l == i:
                        
                        dΣ_kl_dX_ij = X[k,j]*d_1[j]
                    
                    else:
                        continue
                    
                    dUdX[i,j] += M_kl[k,l]*dΣ_kl_dX_ij
                    
                
            
    
    #dUd1
    dUd1 = np.zeros(p)
    
    for j in range(p):
    
        for k in range(d):
                
            for l in range(d): 
                    
                dΣ_kl_dD1_jj = X[k,j]*X[l,j]

                dUd1[j] += M_kl[k,l]*dΣ_kl_dD1_jj

        #adding extra term
        dUd1[j] += d_1[j]/(σ_1)**2
    
   
    #dUd2
    dUd2 = np.zeros(d)
    
    for j in range(d):

        dUd2[j] += M_kl[j,j]*1.
        
        #adding extra term
        dUd2[j] += d_2[j]/(σ_2)**2
    
    pot_derv = np.zeros(int(d*p+p+d))
    pot_derv[:d*p] = matrix_to_vec(dUdX)
    pot_derv[d*p:d*p+p] = dUd1
    pot_derv[d*p+p:] = dUd2
    
    return pot_derv


q_initial = list(matrix_to_vec(np.eye(d,p)))
d_1 = abs(np.random.normal(0,σ_1,p))
d_2 = abs(np.random.normal(0,σ_2,d))
q_initial += list(d_1)
q_initial += list(d_2)
q_initial = np.array(q_initial)
x_init = q_initial
#x_init = q_initial


t1 = time.time()
potential_derv_fast(q_initial)
t2 = time.time()
print('time before compiling = ',t2-t1)
t1 = time.time()
potential_derv_fast(q_initial)
t2 = time.time()
print('time after compiling = ',t2-t1)


#RATTLE Hamiltonian Flow
@numba.jit(nopython=True)
def RATTLE_with_Potential(x0,v0,t,dt,max_elim_iters):
    n = np.floor(t/dt)
    vn = v0
    qn = x0
    vhalf = v0
    G_q = G(qn)
    
    #Gram Matrix is GG^T
    gram = matmul(G_q,G_q.T)
    gram_inv = np.linalg.inv(gram)
    
    
    
    pderv = potential_derv_fast(qn)#potential_derv(qn)
    
    residual_list = np.zeros(num_of_constraints)
    for i in range(int(n)):
        
        
        #solver for Lagrange position multipliers
        Q = qn + vhalf*dt - dt*dt*pderv
        
        #non-linear gaussian elimination
        for k in range(max_elim_iters): #i>j
            for i in range(p):
                for j in range(i,p):
                    g_Q = g_ij(Q,i,j)
                    index = int(i*p - 0.5*i*(i-1) + j-i)
                    
                    residual_list[index] = g_Q
                    if abs(g_Q) < 1e-8:
                        continue
                    G_Q = G(Q)
                    
                    #should be sum of i's and js in indexing below
                    dlambda = g_Q/dot_product(G_Q[index,:],G_q[index,:])
                    Q = Q - G_q[index,:]*dlambda
            #break condition
            if np.all(np.abs(residual_list)<1e-8):
                break

        
        #half step
        vhalf = (Q-qn)/dt
        qn = Q
        
        pderv = potential_derv_fast(qn) #potential_derv(qn)
        G_q = G(qn)
        
        gram = matmul(G_q,G_q.T)
        gram_inv = np.linalg.inv(gram)
        
        #linear solver Lagrange velocity multipliers
        b = matrix_vec_multiplication(G_q,2*vhalf/dt - pderv)
        coeffs_v = matrix_vec_multiplication(gram_inv,b)
        
        #full step
        vn = vhalf - 0.5*dt*pderv - 0.5*dt*matrix_vec_multiplication(G_q.T,coeffs_v)
       
    return qn,vn




#Sampling Event Times
@numba.jit(nopython=True)
def time_exp(lam):
    t = np.random.exponential(lam)
    return t

@numba.jit(nopython=True)
def tangent_space_gaussian(q):
    
    
    v = np.random.normal(0.,1.0,dimension).T
    
    
    G_q = G(q)
    
    
    gram = matmul(G_q,G_q.T)
    gram_inv = np.linalg.inv(gram)
    
    proj_matrix = np.eye(dimension) - matmul(G_q.T,matmul(gram_inv,G_q))
    
    #sample 3d gaussian and then project onto tangent space.
    v = matrix_vec_multiplication(proj_matrix,v)
    
    return v

Z = tangent_space_gaussian(x_init)
t1 = time.time()
RATTLE_with_Potential(q_initial,Z,0.01,0.001,50)
t2 = time.time()
print('time before compiling = ',t2-t1)
t1 = time.time()
RATTLE_with_Potential(q_initial,Z,0.01,0.001,50)
t2 = time.time()
print('time after compiling = ',t2-t1)

@numba.jit(nopython=True)
def U(q):
    
    X = vec_to_matrix(q[:d*p])
    
    d_1 = q[d*p:d*p+p]
    d_2 = q[d*p+p:]
    
    D_1 = np.diag(d_1)
    D_2 = np.diag(d_2)
    
    Σ = matmul(matmul(X,D_1),np.transpose(X)) + D_2
    Σ_inv = np.linalg.inv(Σ)
    
    #likelihood
    pot = 0.5*n*np.log(np.linalg.det(Σ)) + 0.5*n*p*np.log((2*np.pi))
    for i in range(n):
        Bx = matrix_vec_multiplication(Σ_inv,Data[:,i]-mean)
        pot += 0.5*dot_product(Data[:,i]-mean,Bx)
        
    #prior
    #don't need uniform prior because it's constant
    pot += 0.5*np.log(2*np.pi*(σ_1)**2)
    pot += 0.5*dot_product(d_1,d_1)/(σ_1)**2
    
    pot += 0.5*np.log(2*np.pi*(σ_2)**2)
    pot += 0.5*dot_product(d_2,d_2)/(σ_2)**2
    
    return pot


#@numba.jit(nopython=True)
#def f(q):
    
#    z = U(q)
    
#    return z

@numba.jit(nopython=True)
def hamiltonian(x,v):
    return U(x) + 0.5*dot_product(v,v)


@numba.jit(nopython=True)
def f(q):
    X = vec_to_matrix(q[:d*p])
    d_1 = q[d*p:d*p+p]
    d_2 = q[d*p+p:]
    
    D_1 = np.diag(d_1)
    D_2 = np.diag(d_2)
    
    Σ = matmul(matmul(X,D_1),np.transpose(X)) + D_2
    return Σ


def Potential_approx_deriv(q,ϵ):
    δV = np.zeros(len(q))
    for i in range(len(q)):
        h = np.zeros(len(q))
        h[i] = ϵ
        δV[i] = (U(q+h)-U(q-h))/(2*ϵ)
    return δV



#Testing 
q_initial = list(matrix_to_vec(A))

d_1 = np.diag(D_1) #abs(np.random.normal(0,σ_1,p))
d_2 = np.diag(D_2) #abs(np.random.normal(0,σ_2,d))

q_initial += list(d_1)
q_initial += list(d_2)
q_initial = np.array(q_initial)
x_init = q_initial
#x_init = [-0.1189457869072246, -0.1532052405446348, -0.055524795344386466, 0.1318029832358942, -0.06721621082965719, -0.04367800854594563, -0.0890051672947462, 0.041271950008578695, 0.08000937016093752, 0.01895517996621171, 0.29337367442077683, -0.21344638406224942, -0.0055184061026798585, 0.1530809932265405, -0.12823959258700082, -0.04568951668877275, 0.14524942318297487, -0.05935427730271483, 0.13499519900470927, -0.07701941773818943, -0.23413115441039592, 0.1565218835920915, -0.0598183965429737, -0.1021580076316431, 0.09492461552607588, -0.16273249257610414, -0.11671520788416634, -0.02172720442750551, -0.02281899388583262, -0.034701186414861335, -0.044676591037603325, 0.259101083966414, 0.11969591712316144, 0.031757848628513274, -0.21596764069790939, -0.06233738329526475, 0.0951124113044121, -0.2058754403740536, 0.087763605816172, 0.18805428087209272, -0.13284362625677545, 0.02538410991051019, 0.07060478286330155, -0.18345245114765946, 0.08215688997549841, -0.11671754776651086, 0.08617006689249714, -0.020246777792063596, 0.07139764514378885, 0.009273221198617756, 0.18354967409103085, -0.2195564960480616, -0.0377588342333215, 0.08691304082109619, -0.06366531563709497, -0.09037402959201545, 0.31821536894654723, 0.05410308138733849, -0.0030906146939607493, 0.1579931039159656, -0.11145255156516551, 0.012082298585083563, 0.02219342159984371, -0.01171723418748735, 0.034816193188120605, 0.05037213911356894, 0.05827011350249334, 0.14426930305340996, -0.18422205177236295, -0.16163228634013282, 0.24056579257217744, -0.04392755994764904, -0.34773382670577474, -0.1169452811710747, 0.175379895315866, -0.013747429798157098, -0.023363891560940892, 0.035874589605728376, -0.03841504658463743, -0.008646686660342405, 0.051386794801259794, -0.13422500857790012, -0.26698497113113306, 0.00380290842282765, -0.06810057521106717, 0.12373892312696262, -0.2150195818740738, -0.07538942857690274, 0.039322070109323363, 0.08680583727996508, -0.017423098116496175, -0.0027345310182991267, -0.116367323776691, 0.12250449826026595, 0.05191980245948846, -0.2527979851009291, -0.009802896999031896, 0.1142457488181253, 0.12094446710792155, 0.23254305073936374, 0.03578821603127601, -0.06333414344442012, 0.10729623996190643, 0.13574501445877182, -0.016712448807328122, 0.04612080116919837, -0.002193831973109665, 0.25943878432987283, -0.16782536698744793, -0.11631657099185504, 0.20764740818124858, -0.04503705916873586, -0.23604632736841752, -0.08958849466077427, 0.14559599289841313, -0.1103855445566403, -0.034554995932981926, -0.045338061778162714, -0.15313998793443184, -0.06496071657631498, 0.2137374788591269, -0.09597685690044527, -0.0726768376178967, -0.14884244646149639, 0.03407494412063018, 0.12810550944013718, -0.10869513551309705, 0.1517458621209612, 0.0012488406236690005, 0.1068589323695738, 0.2399373262302408, -0.043770209334820835, 0.16447815791030584, -0.3016392899658222, -0.08823529696857174, 0.1858469827755413, 0.155395499992338, 0.08506873070717504, 0.0778794889345394, 0.042455018143264536, 0.16374275578019776, 0.037395749158683586, 0.09557969662961749, 0.026351348138295927, 0.13292997127739606, 0.10976505196045058, 0.002913523757062229, -0.10265957663329911, -0.07842930609118177, 0.08182802522856467, -0.021237544140566295, 0.11186502118107898, 0.02779541108779882, 0.09398381772278333, 0.047351306852443834, 0.16311717978841106, 0.0035539257042423006, 0.0006957353124289155, -0.10233219371327955, -0.07956660838206381, 0.0844846686711103, 0.22186815131244583, 0.08821260293085377, 0.10968116147878121, -0.00933575372867809, 0.33467153965845015, -0.06517122012569318, 0.2186252281689461, 0.003286687191099869, 0.05816065949773853, 0.15256046587606947, -0.0546343059079189, 0.10563006371317148, -0.2567458527427682, -0.10352729805414576, 0.10748150885252959, 0.23821973549835274, 0.04765309253200741, 0.04505103661717917, -0.07060318356080529, -0.13714032476069196, 0.3359154773897516, -0.25909391748501737, 0.010360161666705435, -0.03131485106913592, 0.04647213656971988, 0.055953651493398995, 0.007015720106133663, -0.053581742686872, -0.3894590932107038, -0.009917099049146732, 0.06892157819671739, 0.1647436086358956, 0.2667365787113288, -0.06051824358021412, 0.08077294485279707, 0.04524605853243987, 0.08949583874353752, -0.025138763125361268, 0.011206581248394071, -0.15923578152125226, 0.04221166444378727, 0.14603173520621499, -0.027001439575413434, -0.016968480587828343, 0.13530978105564756, -0.10451125263256905, -0.20007924519005924, -0.01698952731547163, -0.0872359262863817, 0.009458791234824562, -0.09838120148258456, -0.15370963477265992, 0.16094721897320494, -0.04531471754561836, 0.0939978367696399, 0.16344134089409887, -0.002765046298641239, -0.17957327447669016, 0.0005163919882324526, -0.09281296977884759, 0.09640216723008473, 0.09842743954524923, 0.026232505322275123, 0.05609154892782628, 0.1720482321130052, 0.15297858623810798, 0.07967282489303947, -0.012921077329471758, -0.3248363828336502, -0.05105204473096207, -0.03815783244955159, 0.01759951789209686, 0.028282761116092063, -0.22754797609603258, -0.05799529354601056, -0.05635288353041896, 0.0016529348929668454, 0.017223725077309866, 0.007469440947977366, -0.0043167016586454715, 0.04892969227034329, -0.012457612939778304, -0.011054629405167148, 0.06682182329509741, 0.09368828019443366, 0.04962348634057525, 0.054769690830793936, 0.08489318405814095, 0.0977619335843102, 0.17588926924857134, 0.22810771702970295, 0.23916372736367067, 0.28009231847911775, 0.2181954954904803, 0.12247317328064346, -0.023242756104500172, -0.07140671020483949, -0.15268583801918476, -0.17168636080273386, -0.04939962242884995, -0.06913553598067242, -0.044117684637976905, -0.007490889657347759, 0.01902882184319406, -0.08775428492193522, -0.009503844114369001, 0.039306866462171276, 0.03998099201571125, -0.029950436549877016, -0.03002224377194808, -0.03973095539873898, -0.026258173436761442, 0.009462563813939189, 0.05042225354343654, 0.08369654848594256, 0.004945287507275056, 0.11534029882793786, 0.22565759155906737, 0.17574552236408159, 0.052840168547081484, 0.07007601215555888, 0.11742956574431401, 0.020240863701376496, 0.12911702826322682, 0.12791517098847704, 0.11844470240574072, 0.0901570095618336, 0.1536524431119874, 0.16304332489057005, 0.2336198762481063, 0.24687314862343143, 0.2727520223686886, 0.28454928154864356, 0.22474735491574588, 0.15014379340053083, 0.020928344452918073, -0.04866981722713955, -0.1128089154642898, -0.07787787268894997, -0.09547377357445173, 0.0015958403439194863, 0.16701664786371734, 0.02231785796507824, 0.0489753222684834, -0.01468985125528063, 0.1290801036328539, 0.07310228092006323, 0.25046158058717494, 0.006109233467251689, -0.11443330279100228, 0.15231212856964121, -0.25301104416156645, 0.06608570494447573, -0.03975534638263868, 0.04442422193865368, 0.13961523774820533, 0.13085389700229863, 0.11247311100143431, -0.034122796106647736, 0.03762230587824195, 0.10542446760495856, -0.06763038487173337, 0.1036873091236435, -0.08298729052156925, 0.1821171242020173, 0.20670160676544938, 0.13843209134332052, 0.20877179244204483, -0.002646876823382026, 0.004548045356904698, -0.12071492555052396, 0.08708618942726648, 0.1941361134547269, 0.08732126626157526, 0.11632080815307706, 0.08379844196163616, -0.23229310793771926, -0.147985553977299, 0.016533079851948206, -0.09048485445698583, -0.1615261550299452, 0.23981670081160603, 0.012689747766835647, 0.05049542645347869, 0.07045387194986272, 0.06334619841337819, 0.11037618346103048, 0.28005400614546166, 0.10193967159166488, -0.02484471065770341, 0.13746235285461617, -0.22194148602150995, -0.02298884868563884, -0.008365820078870045, 0.03920346484934258, 0.15541680063811278, 0.2639059678667064, -0.0939202196405578, -0.029403489871308894, -0.03705611505764903, -0.026222253653608057, 0.19936002405955008, 0.08816599207982034, 0.013613968049315264, 0.082847371293706, 0.17344276526571467, 0.1533299205608105, 0.2434915722485585, 0.027136190774479364, 0.14117784362276375, 0.11860647250445472, -0.06366495018275814, -0.10099697155198371, -0.14888872911686643, -0.09726663568387513, -0.25488129235929163, -0.197002250716792, 0.09520612903145022, 0.26584517864749785, -0.019277043060112114, 0.013044058460757915, 0.03648298880396426, -0.07067944476693616, 0.1514329842040758, -0.09725215382088093, 0.03058312238201849, -0.06953970946326817, 0.030449334010543255, -0.04909464728894309, 0.07646942524783305, 0.020627553799852007, -0.0795580201853916, 0.029096704596131083, -0.030729231521593712, 0.19325418868306163, 0.10666864690675609, 0.1765105962612155, 0.09456980990505298, -0.052368092997650204, 0.06515367117674888, 0.12353660078324034, 0.11695973227720233, -0.05635474746174571, 0.0828061112595276, -0.021264485704970998, 0.048059161009510624, 0.12701743778774077, 0.15834657997286697, 0.012935717112075026, 0.08925096614632239, 0.08170947283432257, -0.07971205149383383, -0.11016966743306204, -0.1260771277362036, -0.19162177473515735, -0.32165742820274423, -0.16182026323880802, 0.07326912652784441, 0.3223577193246373, 0.06080060891911504, 0.12627405355979157, 0.0724165398381286, -0.12733952014048128, 0.20628546780329424, -0.09994959183323532, 0.011610140360541856, 0.08034532369308839, 0.07417357344103005, 0.11505035921073553, -0.14609461152631723, -0.36416828443350485, -0.15841719518295633, 0.18319110270131364, -0.016250602651026985, 0.2100968994985039, -0.02158105391027714, 0.11155890510643159, -0.032048259607201875, 0.1647359956200387, -0.07767220229244647, -0.11048420539564317, -0.012867601677948093, -0.05032722223963822, 0.11885883447063054, -0.10058468231246245, -0.14936634665631013, -0.06983728380537102, 0.018893793512510762, -0.05839235367009148, 0.11653996515000853, -0.05981837735671615, -0.06928463556415475, 0.03925520549631315, -0.04509409447242199, -0.013863040342862381, -0.14211016615157193, -0.08318357285389405, 0.05761561155430104, -0.12134382781119463, 0.26098592233410206, 0.24397193079658194, 0.14093343002253925, 0.0065020068743601896, 0.027529093706238177, -0.008430617604574465, 0.045315093089664166, 0.046329336199974605, 0.057710776806141596, 0.10288426426108596, -0.12184405434797765, -0.2607928588587531, -0.12081541272808338, 0.16876982256858378, -0.023179381395316792, 0.29634313580433647, -0.09714052927502707, 0.009861075704273544, -0.21287906956679606, 0.016579972970296345, 0.08295334066168697, 0.013512279282732471, 0.009815101893777944, -0.13480772301252042, 0.02300339990216295, 0.0034050740219368734, 0.08468810150103508, -0.0799164258186787, -0.310288282248481, -0.07401574173396251, -0.10249174175654949, -0.014160338284935164, -0.07259710456170029, -0.028384842250152833, -0.02922834039471862, -0.24317912324174584, 0.04639581644737474, -0.034479157122496766, 0.2543436596834627, -0.20605088528981796, -0.14309319968200926, -0.13536764343704247, 0.041089575859877076, 0.057497798201516945, 0.07330175784161007, -0.24382985102034335, -0.06883446979819677, -0.033847478382352976, -0.02991174161837471, 0.10801786502568332, 0.21496215430925703, 0.0345181280400862, -0.10253109011966928, -0.09553894614120116, 0.12154209720911285, 0.1330967862253495, -0.09455753871968338, -0.32459593921470575, 0.030535475557304147, -0.12990811691676363, 0.07264509641947967, -0.002938215052573695, 0.08671083065675222, -0.1536895944459473, 0.1685525450955189, 0.035094854040467444, 0.1263043303166875, 0.2335615770011503, -0.06698930158340023, 0.14815655946646183, 0.1157638341834161, 0.17912125058841694, 0.027515333046413865, 0.031014458783591263, -0.09951633102871223, -0.17317213617155416, -0.09475387848572871, 0.06900329031288888, 0.09366074494683924, -0.20056228693838374, -0.02486038145917191, 0.08510540978916799, -0.12330805746237439, 0.06993593786621459, 0.07676560553619576, 0.11234045751103015, 0.03725026882788774, -0.032628159095522524, -0.1698119341566437, 0.3490999318431106, 0.05544204610083381, 0.3213810430512986, -0.08019402284522222, 0.08932736908340229, -0.053755611277168554, 0.177415231762424, 0.13068297511062524, -0.09345605662153407, 0.05261269281412915, -0.07389530864810451, 0.033267202335175224, 0.019346655455773527, -0.032299037063551875, -0.008086804605875483, -0.03710788277208559, 0.07131076447109587, -0.1489116578774142, -0.1304024176770533, 0.0793099698232237, -0.12115816361088144, -0.14980020343058914, -0.11350031594830416, -0.20545425010131174, -0.19684816188483656, -0.15683347107759557, 0.0009497462734653428, -0.06656569584896783, -0.020824186175784326, -0.04753123270567543, 0.1512541457162298, 0.027296556153128686, 0.08597200711116562, -0.0375293415804397, -0.0575164259740855, 0.08705218701128858, -0.07358012266644907, -0.06200371219516268, -0.03479270625572368, -0.2474575452379583, 0.23806212455249134, -0.07858938610695192, 0.09844928551788536, -0.26092260681882784, -0.07438011897082811, -0.2129208666401803, 0.11689165258951256, 0.004452050469109774, 0.1329241200247032, 0.14768053575924275, 0.1852663034584991, 1.1207273216083515, 1.8274975292107254, 1.4557542386681672, 1.2640455596685851, 5.3510743041351505, 1.3853426892108813, 3.5633544313543744, 1.497934322045999, 0.8859428327114853, 1.082438606578969, 0.4247793710700938, 0.7984196107206972, 0.6829414212078465, 0.6693804241783082, 0.6799587808133589, 0.6708514508206149, 0.7778305112766691, 0.4697978335418313, 0.028839277064358643, 0.012223446492642294, 0.048177680424825786, 0.004294758376451444, 0.060980029304455195, 0.0035854860755611924, 0.08632980662856646, 0.2691454113764132, 0.18884579837029847, 0.6479830837196807, 0.3772426078257839, 0.6638246763457016, 0.6891964310654277, 0.8063503628047305, 0.5974796112681726, 1.0228700915764926, 0.6598790370948376, 1.0152643326025161, 0.5326215219795821, 1.7339821861628415, 0.8383880278245036, 0.4758490026764131, 0.7567880535822031, 1.1915771735461655, 0.8122736260704951, 0.762034680889643, 1.2845058926498791, 0.6929236922978466, 0.8537800014249455, 0.801077531729846, 0.8107673743271597, 0.5777054089103292, 0.45203505125274773, 0.7499820852889928, 0.754926695082765, 0.6371786116779825, 0.8846086950979487, 0.47271301161762536, 0.6105326265696221, 0.4054142899868186, 0.019347473155107334, 0.010258732473165102, 0.00720307564684286, 0.006990692966921607, 0.005734996208276974, 0.012247854582117764, 0.004944304408322458, 0.44445006192029557, 0.03578861808051631, 0.4390208420784934, 0.6849384240365435, 0.8105196044679795]

x_init = np.array(x_init)


print(U(x_init))
t1 = time.time()
p1 = potential_derv_fast(x_init)
t2 = time.time()
p2 = Potential_approx_deriv(x_init,0.000000001)
t3 = time.time()
print('sup|grad_exact - grad_num| =',max(p2-p1))
print('Exact time = ',t2-t1)
print('Numerical time = ',t3-t2)


#Initialise
T = 0.01
num_of_events = 3
dt = 0.001


@numba.jit(nopython=True)
def RRHMC(num_of_events,dt_max,T,x_init):
    
    #Exponential Expected Value
    rate = T
    x = x_init
    
    position_list = [x_init]
    v = tangent_space_gaussian(x)
    
    
    for i in range(num_of_events):
        
        t = time_exp(rate)
        L = np.ceil(t/dt_max)
        
        dt = t/L
        
        h = hamiltonian(x,v)
    
        xnew,vnew = RATTLE_with_Potential(x,v,t,dt,1000)
        
        h_new = hamiltonian(xnew,vnew)
        
        #metropolis hasting step
        #adding in rejection for non-positive diagonal matrices.
        d_1 = xnew[d*p:d*p+p]
        d_2 = xnew[d*p+p:]
        
        u = np.random.rand()
        if u <= np.exp(-h_new+h) and min(d_1)>0 and min(d_2)>0:
            x = xnew
            

        
        position_list.append(x)
        
        v = tangent_space_gaussian(x)
        if i%10==0:
            print(i)
        
        
    return position_list


@numba.jit(nopython=True)
def RHMC(num_of_events,dt,T,x_init):
    
#     x = list(matrix_to_vec(np.eye(d,p)))
#     x += list(np.ones(d + p))
    x = x_init
    #initialisation of x on V_{d,p} \times \mathbb{R}^p 
    #                        \times \mathbb{R}^d
    v = tangent_space_gaussian(x)
    
    position_list = [x]
    
    
    for i in range(num_of_events):
        
        h = hamiltonian(x,v)
        
        xnew,vnew = RATTLE_with_Potential(x,v,T,dt,1000)
        
        h_new = hamiltonian(xnew,vnew)
        
        #metropolis hasting step
        #with reject with positivity constraint
        d_1 = xnew[d*p:d*p+p]
        d_2 = xnew[d*p+p:]
        
        
        u = np.random.rand()
        if u <= np.exp(-h_new+h) and min(d_1)>0 and min(d_2)>0:
            x = xnew
        
        position_list.append(x)
        v = tangent_space_gaussian(x)  
        if i%10==0:
            print(i)
            
    return position_list


t = time.time()
position = RRHMC(num_of_events,dt,T,x_init)
elapsed = time.time() - t
print('RT time =',elapsed)
t = time.time()
position = RHMC(num_of_events,dt,T,x_init)
elapsed = time.time() - t
print('DT time =',elapsed)



#Initialise
#T = 0.1
num_of_events = 10000
dt = 0.001
N = 10
T = N*dt
#ensure positivity
#log normal prior

position = RRHMC(num_of_events,dt,T,x_init)

first_entry_list = []
first = 0
burn = 20000
for i in range(len(position)-burn):
    mat = f(position[i+burn])
    
    first += mat[-1,0]
    first_entry_list.append(first/(i+1))
    

last_entry_list = []
last = 0
burn = 0
for i in range(len(position)-burn):
    mat = f(position[i+burn])
    
    last += mat[0,-1]
    last_entry_list.append(last/(i+1))



burn = 0
Σ_ave = np.zeros((d,d))
for i in range(num_of_events-burn):
    Σ_ave += f(position[i+burn])/(num_of_events-burn)
    #print(len(position[i+burn]))
    
S = Σ_ave


#4 before
#plot increase in maximum element
