import numpy as np
import math


def getDCM(globalAxisBasisVectors,newCoordinateSystemVectors):
    dcm=np.zeros([3,3])
    for iRow in range(3):
        for iColumn in range(3):
            dcm[iRow,iColumn]=np.dot(newCoordinateSystemVectors[iRow],globalAxisBasisVectors[iColumn])

    # check if matrix is orthogonal
    dot_product = np.dot(dcm, dcm.T)
    identity_matrix = np.identity(len(dcm))
    
    matrixIsOrthogonal=np.allclose(dot_product, identity_matrix,atol=0.01)
    
    if matrixIsOrthogonal:
        return dcm
    else:
        raise ValueError(f'Direction cosine matrix is not orthogonal {dot_product}')
    

def dcmToEulerAngles(dcm):

    '''
    Formulas from:
    Euler Angle Formulas
    David Eberly, Geometric Tools, Redmond WA 98052
    https://www.geometrictools.com/
    '''
    #RzRyRx
    if round(dcm[2,0],3)<1:   
        if round(dcm[2,0],3)>-1:
            theta2=math.asin(-dcm[2,0])
            theta3=math.atan2(dcm[1,0],dcm[0,0])
            theta1=math.atan2(dcm[2,1],dcm[2,2])
        else: #r20 =-1
            theta2=math.pi/2
            theta3=-math.atan2(-dcm[1,2],dcm[1,1])
            theta1=0
    else:#r20 =1
        theta2=-math.pi/2
        theta3=math.atan2(-dcm[1,2],dcm[1,1])
        theta1=0
    return theta1*180/math.pi,theta2*180/math.pi,theta3*180/math.pi