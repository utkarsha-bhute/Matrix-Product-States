# -*- coding: utf-8 -*-
"""
Created on Wed May 21 16:54:08 2025

@author: utkar
"""

import numpy as np

file1 = open("test_file.txt", "w")

'''This function normalizes the state vector v'''
def normalize(v):
    c = 0
    for i in range (0, len(v)):
        c = c + (v[i]*v[i])
    v = v / np.sqrt(c)
    return v

'''This function returns an array of left normalized A matrices for the normalized state vector vin'''
def leftCanonicalMPS(d, L, vin, count):
    global i1, i2, dc
    if (count <= L-1):

        A0, A1 = np.linalg.qr(vin.reshape((int(len(vin[0])/(d**(L-count))), d**(L-count))), mode='reduced')  

        for i1 in range (0, d):
            dc = i1
            A = []
            
            x = d**(count-1)
            if (count > (L+1)//2):
                x = d**(L+1-count)
                
            for i2 in range (0, int(len(A0)/d)):
                A.append(A0[dc])
                dc = dc + d
            matricesA.append(np.array(A).reshape(x, int(len(np.array(A).reshape((1,-1))[0])/x)))

        if (count == L-1):
            A1 = A1.T
            for i1 in range (0, d):
                matricesA.append(A1[i1].reshape((-1,1)))
                
        count = count + 1
        leftCanonicalMPS(d, L, A1.reshape((1,-1)), count)
    return matricesA
    
'''This function returns an array of right normalized B matrices for the normalized state vector vin'''
def rightCanonicalMPS(d, L, vin, count):
    global i1, i2, dc
    if (count <= L-1):
        
        B0, B1 = np.linalg.qr(vin.reshape(( d**(L-count), int(len(vin[0])/(d**(L-count))))), mode='reduced')

        for i1 in range (0, d):
            dc = i1
            B = []
            
            x = d**(count-1)
            if (count > (L+1)//2):
                x = d**(L+1-count)
  
            B1 = B1.reshape(int(len(B1.reshape((1,-1))[0])/x), x)
            for i2 in range (0, int(len(B1)/d)):
                B.append(B1[dc])
                dc = dc + d
            matricesB.append(np.array(B).reshape(int(len(np.array(B).reshape((1,-1))[0])/x), x))
            
        if (count == L-1):
            for i1 in range (0, d):
                matricesB.append(B0[i1].reshape((1,-1)))  
             
        count = count + 1
        rightCanonicalMPS(d, L, B0.reshape((1,-1)), count)
    return matricesB

'''This function is used to from the matrices of required dimensions in MCMPS ONLY'''
def appendA(A0, x, count, p):
    if (count > p):
        A0 = A0.reshape(int(len(A0.reshape((1,-1))[0])/x), x)
    for i1 in range (0, d):
        dc = i1
        A = []
        for i2 in range (0, int(len(A0)/d)):
            A.append(A0[dc])
            dc = dc + d
        y = int(len(np.array(A).reshape((1,-1))[0])/x)
        if (count <= p):
            matrices.append(np.array(A).reshape((x, y)))
        if (count > p):
            matrices.append(np.array(A).reshape((y, x)))

'''This function returns an array of left normalized A matrices, singular matrix and 
right normalized B matrices for the normalized state vector vin for partition at p
for 5 qubit state p = 2 means [12]:[345] partition'''
def mixedCanonicalMPS(d, L, vin, count, p):
    if (count <= L-1):

        if (count <= p):
            A0, A1, A2 = np.linalg.svd(vin.reshape((int(len(vin[0])/(d**(L-count))), d**(L-count))), full_matrices=False)  
            A1 = np.diag(A1)
            
            x = d**(L-count+1)
            if (count <= (L+1)//2):
                x = int(d**(count-1))
            appendA(A0, x, count, p)

            if (count < p):
                A2 = np.dot(A1, A2)

            if (count == p):
                matrices.append(A1) 
                
                if (p == L-1):
                    appendA(A2, d, count, p)
                    return matrices
                    
            count = count + 1
            mixedCanonicalMPS(d, L, A2.reshape((1,-1)), count, p)           

        if (count > p):
            x = d**(L-count+p)
            
            if (p > L//2):
                if (L%2 == 0):
                    x = d**(L+p-count-(2*(p-L//2)))
                if (L%2 == 1):
                    x = d**(L+p-count-(2*(p-L//2)))
                    
            y = int(len(vin[0])/x)

            A0, A1, A2 = np.linalg.svd(vin.reshape((x, y)), full_matrices=False)
            A1 = np.diag(A1)
            A0 = np.dot(A0, A1)
            
            appendA(A2, int(y//d), count, p)
            if (count == L-1):
                if (p >= L//2):
                    appendA(A0, y, count, p)
                if (p < L//2):
                    appendA(A0, x, count, p)
            count = count + 1
            mixedCanonicalMPS(d, L, A0.reshape((1,-1)), count, p)
    
    return matrices

'''This function gives overlap between two states'''
def overlap(d, L, z1, z2):
    d1 = np.zeros((d,d))
    for i in range (0,d):
        d1 = d1 + np.dot(np.transpose(np.conjugate(z2[i])),z1[i])

    for i0 in range (1, L-1):
        x = d**(L-(i0+1))
        if (i0 < L//2):
            x = d**(i0+1)   
        elif (L % 2 == 1 and i0 == L//2):
            x = d**(i0)  
 
        d2 = np.zeros((x,x))
        for i in range (0,d):
            d2 = d2 + np.dot(np.dot(np.transpose(np.conjugate(z2[i+i0*d])), d1), z1[i+i0*d])
        d1 = d2

    dn = np.zeros((1,1))
    for i in range (0,d):
        dn = dn + np.dot(np.dot(np.transpose(np.conjugate(z2[i+(L-1)*d])), d2),z1[i+(L-1)*d])

    return dn[0][0] 
    
'''This function gives matrix elements of operator O'''  
def matrixEle(d, L, z1, z2, O):
    print("\n\n\n\n")
    d1 = np.zeros((d,d))
    for i in range (0,d):
        for j in range (0,d):
            d1 = d1 + O[0][i][j] * np.dot(np.transpose(np.conjugate(z2[i])),z1[i])
    
    if (L == 2):
        d2 = d1
    else:
        for i0 in range (1, L-1):
            x = d**(L-(i0+1))
            if (i0 < L//2):
                x = d**(i0+1)   
            elif (L % 2 == 1 and i0 == L//2):
                x = d**(i0)  
     
            d2 = np.zeros((x,x))
            for i in range (0,d):
                for j in range (0,d):
                    d2 = d2 + O[i0][i][j] * np.dot(np.dot(np.transpose(np.conjugate(z2[i+i0*d])), d1), z1[i+i0*d])
            d1 = d2

    dn = np.zeros((1,1))
    for i in range (0,d):
        for j in range (0,d):
            dn = dn + O[L-1][i][j] * np.dot(np.dot(np.transpose(np.conjugate(z2[i+(L-1)*d])), d2), z1[j+(L-1)*d])
            
    return dn[0][0]
    
'''This function gives expectation of operator O'''  
def expValue(d, L, z1, O, i0):
    dn = 0
    for i in range (0,d):
        for j in range (0,d):
            dn = dn + O[i][j] * np.trace(np.dot(np.transpose(np.conjugate(z1[i+i0*d])), z1[j+i0*d]))
            
    return dn  
 
'''This function gives reduced density matrix for partition A'''  
def rdmA(d, L, z1, l):
    rhol = np.zeros((d,d))  
    for i in range (0,d):
        rhol = rhol + np.dot(z1[d*L-i-1], np.transpose(np.conjugate(z1[d*L-i-1])))    
    for i0 in range (1, L-l):
        x = d**(L-(i0+1))
        if (i0 < L//2):
            x = d**(i0+1)   
        elif (L % 2 == 1 and i0 == L//2):
            x = d**(i0)  
        d2 = np.zeros((x,x))
        for i in range (0,d):
            d2 = d2 + np.dot(np.dot(z1[d*L-1-i-i0*d], rhol), np.transpose(np.conjugate(z1[d*L-1-i-i0*d])))
        rhol = d2
        
    rangelist = []
    num = 0
    for i in range (1, l+1):
        templist = []  
        for j in range (0, d):
            templist.append(num)
            num = num +1
        rangelist.append(templist)
    
    rdm = []
    def for_recursive(number_of_loops, range_list, current_index=0, iter_list = []):
        if (iter_list == []):
            iter_list = [0]*number_of_loops
        if current_index == number_of_loops-1:
            for iter_list[current_index] in range_list[current_index]:
                iter_list = np.array(iter_list)
                result = np.identity(1)
                for i in range(0, len(iter_list)):
                    result = np.dot(result, z1[iter_list[i]])  
                AmatrixMul.append(result)
        else:
            for iter_list[current_index] in range_list[current_index]:
                for_recursive(number_of_loops, iter_list = iter_list, range_list = range_list,  current_index = current_index+1) 
    
    for_recursive(range_list = rangelist, number_of_loops = l)
    for i in range (0, len(AmatrixMul)):
        for j in range (0, len(AmatrixMul)):
            rdm.append(np.dot(np.dot(AmatrixMul[i], rhol), np.transpose(np.conjugate(AmatrixMul[j])))[0][0])
    rdm = np.array(rdm).reshape((d**l, d**l))
    #print(rdm,'\n')
     
    sum1 = 0
    for i in range (len(rdm)):
        sum1 = sum1 + rdm[i][i]   
    #print(sum1)    
        
'''This function gives reduced density matrix for partition B'''  
def rdmB(d, L, z1, l):
    rhol = np.zeros((d,d))  
    for i in range (0,d):
        rhol = rhol + np.dot(np.transpose(np.conjugate(z1[d*L-i-1])), z1[d*L-i-1])   
    
    for i0 in range (1, L-l):
        x = d**(L-(i0+1))
        if (i0 < L//2):
            x = d**(i0+1)   
        elif (L % 2 == 1 and i0 == L//2):
            x = d**(i0)  
        d2 = np.zeros((x,x))
        for i in range (0,d):
            d2 = d2 + np.dot(np.dot(np.transpose(np.conjugate(z1[d*L-1-i-i0*d])), rhol), z1[d*L-1-i-i0*d])
        rhol = d2

    def for_recursive(number_of_loops, range_list, current_index=0, iter_list = []):
        if (iter_list == []):
            iter_list = [0]*number_of_loops
        if current_index == number_of_loops-1:
            for iter_list[current_index] in range_list[current_index]:
                iter_list = np.array(iter_list)
                result = np.identity(1)
                for i in range(0, len(iter_list)):
                    result = np.dot(result, np.transpose(np.conjugate(z1[iter_list[i]])))  
                BmatrixMul.append(result)
        else:
            for iter_list[current_index] in range_list[current_index]:
                for_recursive(number_of_loops, iter_list = iter_list, range_list = range_list,  current_index = current_index+1) 
    
    rangelist = []
    num = 0
    for i in range (1, l+2):
        templist = []  
        for j in range (0, d):
            templist.append(num)
            num = num +1
        rangelist.append(templist)
        
    rdm = []
    for_recursive(range_list = rangelist, number_of_loops = L-l)

    for i in range (0, len(BmatrixMul)):
        for j in range (0, len(BmatrixMul)):
            rdm.append(np.dot(np.dot(BmatrixMul[i], rhol), np.transpose(np.conjugate(BmatrixMul[j])))[0][0])
    rdm = np.array(rdm).reshape((d**(L-l), d**(L-l)))
    #print(rdm,'\n')
     
    sum1 = 0
    for i in range (len(rdm)):
        sum1 = sum1 + rdm[i][i]   
    #print(sum1)    
        
''''For LCMPS'''

L = 3
d = 3

psir = [0]*d**L
psir[0] = 1
psir[13] = 1
psir[26] = 1

matricesA = []  
vin = normalize(psir)
               
matrices = leftCanonicalMPS(d, L, vin.reshape((1,-1)), 1)

file1.write("\nA Matrices:")
for i in range(len(matrices)):
    file1.write("\nA" + str(i) +" = \n" + str(matrices[i]))

'''For RCMPS'''

L = 4
d = 2

psir = [0]*d**L
psir[0] = 1
psir[15] = 1
psir[2] = 0
 
matricesB = []  
vin = normalize(psir)                  
 
matricesB = rightCanonicalMPS(d, L, vin.reshape((1,-1)), 1)

file1.write("\nB Matrices:")
for i in range(len(matricesB)):
    file1.write("\nB" + str(i) +" = \n" + str(matricesB[i]))

'''For MCMPS'''

p = 2
L = 3
d = 3

psir = [0]*d**L
psir[0] = 1
psir[13] = 1
psir[26] = 1

matrices = []  
vin = normalize(psir)                                                      

matrices = mixedCanonicalMPS(d, L, vin.reshape((1,-1)), 1, p)

file1.write("\nMatrices:")
for i in range(0,d*L+1):
    if (i < d*p):
        file1.write("\nA" + str(i) +" = \n" + str(matrices[i]))
    if (i == d*p):
        file1.write("\nS" +" = \n" + str(matrices[i]))
    if (i > d*p):
        file1.write("\nB" + str(i) +" = \n" + str(matrices[i]))
  
'''For Overlap'''

L = 3
d = 2

matricesA = [] 
psir1 = [0]*d**L
psir1[0] = 1
psir1[7] = 1

z1 = leftCanonicalMPS(d, L, normalize(psir1).reshape((1,-1)), 1) 

matricesA = [] 
psir2 = [0]*d**L
psir2[1] = 1
psir2[2] = 1
psir2[3] = 1

z2 = leftCanonicalMPS(d, L, normalize(psir2).reshape((1,-1)), 1)

# for i in range(len(z1)):
#     print(z1[i])
# print('\n')
# for i in range(len(z2)):
#     print(z2[i])
    
# #checking the inner product
# print(np.dot(np.transpose(np.conjugate(normalize(psir2))), normalize(psir1)))

c = overlap(d, L, z1, z2)
file1.write("\n\nOverlap:\n" + str(c) + "\n\n")


'''For matrix elements'''

L = 2
d = 2

matricesA = [] 
psir1 = [0]*d**L
psir1[1] = 1
psir1[2] = 1
psir1[3] = 1

z1 = leftCanonicalMPS(d, L, normalize(psir1).reshape((1,-1)), 1) 

matricesA = [] 
psir2 = [0]*d**L
psir2[1] = 1
psir2[2] = 1
psir2[0] = 1

z2 = leftCanonicalMPS(d, L, normalize(psir2).reshape((1,-1)), 1)

'''
O1 = np.array([[0,1],[1,0]])
O2 = np.array([[0,1],[1,0]])
O3 = np.array([[0,1],[1,0]])
'''
O1 = np.array([[1,0],[0,-1]])
O2 = np.array([[1,0],[0,-1]])
O3 = np.array([[1,0],[0,-1]])

O = [O1, O2]
print(O)

c = matrixEle(d, L, z1, z2, O)
file1.write("\n\nMatrix element:\n" + str(c) + "\n\n")
# print(c)


'''For expectation value'''

L = 3
d = 2
i0 = 2

matricesA = [] 
psir1 = [0]*d**L
psir1[1] = 1
psir1[5] = 1
psir1[6] = 1

z1 = leftCanonicalMPS(d, L, normalize(psir1).reshape((1,-1)), 1) 

O =np.array([[0,1],[1,0]])

c = expValue(d, L, z1, O, i0)
file1.write("\n\nExpectation value:\n" + str(c) + "\n\n")
# print(c)


'''For rdm A'''

l = 1
L = 3
d = 2
psir = [0]*d**L
psir[1] = 1
psir[2] = 1
psir[4] = 1
psir[3] = 0
    



matricesA = []
AmatrixMul = []
z1 = leftCanonicalMPS(d, L, normalize(psir).reshape((1,-1)), 1)           
rdmA(d, L, z1, l)   

'''For rdm B'''


l = 1
L = 3
d = 2
psir = [0]*d**L
psir[1] = 1
psir[2] = 1
psir[4] = 1
psir[3] = 0

matricesB = []
BmatrixMul = []
z2 = rightCanonicalMPS(d, L, normalize(psir).reshape((1,-1)), 1)           
rdmB(d, L, z2, l)          
            
   
file1.close()   

