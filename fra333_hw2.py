#!/usr/bin/python3
from calendar import prmonth
from resource import prlimit
import numpy as np
import math
from HW2_utils import FKHW2

'''
    Name:   <SUPAWADEE> <CHANSANGDEE> <62340500056>
            <ONNALIN> <ARUN> <62340500056>
    Date:   <30-10-2022>
    Description:
'''

# Question 1
def endEffectorJacobianHW2(q):
    # R,P,R_e,P_e = FKHW2([3.8212593425706647, 5.0119152480887, 3.7921922310016782])
    R,P,R_e,P_e = FKHW2(q)  # call function FKHW2
    R0_1 = R[:,:,1]         # rotation matrix from frame 0 frame 1
    R0_2 = R[:,:,2]         # rotation matrix from frame 0 frame 2
    ray = np.array([[0.],[0.],[1.]])    # array z of frame i relative with frame i
    R01_z = R0_1.dot(ray)   # dot product of R01 & ray
    R02_z = R0_2.dot(ray)   # dot product of R02 & ray
    d0_1 = P[:,1]           # position vector frame 0 to frame 1 of rightmost column (size: 3x1) of homogeneous transformation from frame 0 to frame 1
    d0_2 = P[:,2]           # position vector frame 0 to frame 2 of rightmost column (size: 3x1) of homogeneous transformation from frame 0 to frame 2
    d0_3 = P[:,3]           # position vector frame 0 to frame 3 of rightmost column (size: 3x1) of homogeneous transformation from frame 0 to frame 3

    Jw = np.array([[],[],[]])   #initial array of Jw
    Jw  = np.column_stack((Jw,ray))     # add column ray in array Jw
    Jw  = np.column_stack((Jw,R01_z))   # add column R01_z in array Jw
    Jw  = np.column_stack((Jw,R02_z))   # add column R02_z in array Jw
    # print(Jw)

    Jv = np.array([[],[],[]])   #initial array of Jv
    Jv  = np.column_stack((Jv,(np.cross(ray.T,d0_3)).T))    # ray transpose cross product with d0_3 then transpose all, add in column array Jv 
    Jv  = np.column_stack((Jv,(np.cross(R01_z.T,(d0_3-d0_1))).T))    # R01_z transpose cross product with d0_3 minus d0_1 then transpose all, add in column array Jv 
    Jv  = np.column_stack((Jv,(np.cross(R01_z.T,(d0_3-d0_2))).T))    # R01_z transpose cross product with d0_3 minus d0_2 then transpose all, add in column array Jv 
    # print(Jv)

    J_e = np.concatenate((Jw,Jv),axis=0)    # Combine Jw & Jv
    J_e = J_e.tolist()    # Change Type array >> list
    # print(type(J_e))    # Check type of J_e
    # print(J_e)
    return J_e  
    '''
        q : format list 1x3 [[i_11, i_12, i_13]]
        q unit: rad
        type something here

        return format list 6x3
        [ [i_11, i_12, i_13],
          [i_21, i_22, i_23],
          [i_31, i_32, i_33],
          [i_41, i_42, i_43],
          [i_51, i_52, i_53],
          [i_61, i_62, i_63] ]
    '''
    
# Question 2
def checkSingularityHW2(q):
    # q = [3.1374321176951665, 5.668759349057154, 0.8985620903977458]
    R,P,R_e,P_e = FKHW2(q)  # call function FKHW2
    q1 = q[0]     # theta1
    q2 = q[1]     # theta2
    q3 = q[2]     # theta3

    l1 = 0.0892   # length of link1 (meter)
    l2 = 0.425    # length of link2 (meter)
    l3 = 0.47443  # length of link3 (meter)

    J_star = np.array([[-(l1*np.sin(q1)+l2*np.sin(q1+q2)+l3*np.sin(q1+q2+q3)), -(l2*np.sin(q1+q2)+l3*np.sin(q1+q2+q3)), -l3*np.sin(q1+q2+q3)],
                        [(l1*np.cos(q1)+l2*np.cos(q1+q2)+l3*np.cos(q1+q2+q3)), (l2*np.cos(q1+q2)+l3*np.cos(q1+q2+q3)), l3*np.cos(q1+q2+q3)],
                        [1,1,1]])       #reduced jacobian matrix
    det_J_star = np.linalg.det(J_star)  #det J_star
    
    if abs(det_J_star) < 0.001:     # if absolute of det J_star < 0.001
        flag = True                 # let flage is True (singularity)
    else:                           # else absolute of det J_star < 0.001
        flag = False                # let flage is false (non singularity)
    
    return flag
    '''
        q : format list 1x3 [[i_11, i_12, i_13]]
        q unit: rad
        type something here
        
        return format bool
    '''

# Question 3
def computeEffortHW2(q,w):
    J_e = endEffectorJacobianHW2(q) # call function endEffectorJacobianHW2
    J_e = np.array(J_e).T   # J_e transpose
    R,P,R_e,P_e = FKHW2(q)  # call function FKHW2
    n = np.dot(R_e,np.array(w[:3]))  # R_e dot product with n (moment)
    f = np.dot(R_e,np.array(w[3:]))  # R_e dot product with f (force)
    nf = np.concatenate((n,f),axis=0)  # Combine n & f 
    tau = np.dot(J_e,nf).tolist()      # # J_e transpose dot product with nf, and Change Type array >> list
    return tau

    '''
        q : format list 1x3 [[i_11, i_12, i_13]]
        q unit: rad
        type something here

        return format list 1x3
        [ [i_11, i_12, i_13] ]
    '''
