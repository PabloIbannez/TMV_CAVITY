import sys
import numpy as np

GPa2kcalmol_per_angstrom3 = 0.14393
kcalmol_per_angstrom3_2GPa = 1.0 / GPa2kcalmol_per_angstrom3
kcalmol_per_angstrom3_2MPa = 1000.0 / GPa2kcalmol_per_angstrom3

traj   = sys.argv[1]

minStress = 0
maxStress = 10

def interval2bgr(minVal, maxVal, val):
    #Compute the interval
    interval = maxVal - minVal
    if val < minVal:
        val = minVal
    if val > maxVal:
        val = maxVal

    #Compute the position of val in the interval
    pos = (val - minVal) / interval

    #Compute the rgb value
    r = 255.0 * pos
    g = 0.0
    b = 255.0 * (1.0 - pos)

    #Return the bgr value in hex
    bgr = '#%02x%02x%02x' % (int(b), int(g), int(r))

    #bgr to int
    return int(bgr[1:],16)

if len(sys.argv) > 2:
    stress = sys.argv[2]
else:
    stress = 'Mises'

def checkSymmetry(xx,xy,xz,yx,yy,yz,zx,zy,zz):
    #Load into numpy matrix
    stressMatrix = np.matrix([[xx,xy,xz],[yx,yy,yz],[zx,zy,zz]])
    #Check symmetry
    if np.allclose(stressMatrix,stressMatrix.T):
        return True
    else:
        return False

def computeI1(xx,xy,xz,yx,yy,yz,zx,zy,zz):
    return xx + yy + zz
def computeSigma1(xx,xy,xz,yx,yy,yz,zx,zy,zz):
    #We assume that the stress tensor is symmetric, this avoid numerical errors
    stressMatrix = np.matrix([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])

    #Compute eigenvalues
    eigenvalues = np.linalg.eigvals(stressMatrix)

    #Return the maximum eigenvalue
    return np.max(eigenvalues)


def computeMises(xx,xy,xz,yx,yy,yz,zx,zy,zz):

    I1 = computeI1(xx,xy,xz,yx,yy,yz,zx,zy,zz)
    I2 = xx*yy + xx*zz + yy*zz - xy*yx - xz*zx - yz*zy

    return np.sqrt(I1*I1 - 3*I2)

def computeTresca(xx,xy,xz,yx,yy,yz,zx,zy,zz):

    #We assume that the stress tensor is symmetric, this avoid numerical errors
    stressMatrix = np.matrix([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])

    #Compute eigenvalues
    eigenvalues = np.linalg.eigvals(stressMatrix)

    sigma1 = np.max(eigenvalues)
    sigma3 = np.min(eigenvalues)

    return 0.5*(sigma1 - sigma3)

if stress == "I1":
    computeStress = computeI1
elif stress == "sigma1":
    computeStress = computeSigma1
elif stress == "Mises":
    computeStress = computeMises
elif stress == "Tresca":
    computeStress = computeTresca
else:
    print("Stress not available. Options are: I1, sigma1, Mises, Tresca")
    sys.exit(1)


with open(traj) as f:
    frameCount    = -1
    particleCount = 0
    for line in f:
        if "#" in line:
            frameCount += 1
            particleCount = 0
            print(line, end='')
        else:
            try:
                x,y,z,r,c,vol,xx,xy,xz,yx,yy,yz,zx,zy,zz = [float(val) for val in line.split()]
            except ValueError:
                print("Error (",sys.exc_info()[0],") reading line: ", line)
                sys.exit(1)

            #if not checkSymmetry(xx,xy,xz,yx,yy,yz,zx,zy,zz):
            #    print("WARNING: Detected non-symmetric stress tensor at frame %d particle %d" % (frameCount, particleCount), file=sys.stderr)

            stressObserved = computeStress(xx,xy,xz,yx,yy,yz,zx,zy,zz)

            if vol > 0:
                s = (stressObserved/vol)*kcalmol_per_angstrom3_2MPa
                #s = interval2bgr(minStress, maxStress, s)
                print(x,y,z,r,c,s)
                #print(x,y,z,r,np.power(vol*3.0/(4.0*np.pi),1.0/3.0))
                particleCount += 1
            else:
                s = 0
                print(x,y,z,r,c,s)








