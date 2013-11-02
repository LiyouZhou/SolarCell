from matplotlib import pyplot as plt
from solarcell import SolarCell
from matplotlib.pylab import np
from scipy import optimize
import time


c1 = SolarCell()
c2 = SolarCell()

# consider by-pass diode
#d1 = Diode()
#d2 = Diode()

def fplot(f, limits, transp = 0, args = None, style = '-', label = None):
    [low, high] = limits
    X = np.arange(low, high, (high-low)/100)
    Y = []
    for x in X:
        x = [x]
        if args != None: x = x+args
        Y.append(f(*x))
    if transp: plt.plot(Y, X, style, hold = True , label = label)
    else: plt.plot(X, Y, style, hold = True, label = label )

# function to find MPP
def find_mpp( cell, T, I_L ):
    # minimise -V*I
    def obj_fun(V): return -V*cell.current(V,T, I_L)
    res = optimize. minimize_scalar(obj_fun)
    V_mpp = res.x
    I_mpp = cell.current(V_mpp, T, I_L)
    return V_mpp*I_mpp


# for cells in serise
def in_serise( V, T1, T2, I_L1, I_L2 ):
    def obj_fun(V1):
        V2 = V - V1
        I_D1 = c1.current(V1,T1,I_L1)
        I_D2 = c2.current(V2,T2,I_L2)
        return I_D1-I_D2
    V1 = optimize.newton( obj_fun, 0.8 )
    return V1

def in_parallel( V, T1, T2, I_L1, I_L2 ):
    I1 = c1.current(V,T1,I_L1)
    I2 = c2.current(V,T2,I_L2)
    return I1 + I2

# plot the resultant curves
def draw_axis():
    plt.grid(True)
    plt.axhline(0, linewidth=2, color='black')
    plt.axvline(0, linewidth=2, color='black')

[T1, T2, I_L1, I_L2] = [300,300,0.1,0.05]

plt.subplot(221)
plt.title('Serise')
draw_axis()
fplot(c1.current, [-0.8, 0.8], args=[T1,I_L1], label='Cell 1', style='b')
fplot(c2.current, [-0.8, 0.8], args=[T2,I_L2], label='Cell 2', style='g')

plt.subplot(222)
plt.title('Parallel')
draw_axis()
fplot(c1.current, [-0.8, 0.8], args=[T1,I_L1], label='Cell 1', style='b')
fplot(c2.current, [-0.8, 0.8], args=[T2,I_L2], label='Cell 2', style='g')

plt.subplot(223)
plt.title('Power')
draw_axis()
mpp1 = find_mpp(c1, T1, I_L1)
mpp2 = find_mpp(c2, T2, I_L2)
plt.axhline(mpp1, linewidth=2, color='b', label = 'Cell 1 Max' )
plt.axhline(mpp2, linewidth=2, color='g', label = 'Cell 2 Max')
plt.axhline(mpp1 + mpp2, linewidth=2, color='r', label = 'Cell 1+2 Max')

VI_serise   = [[],[],[]]
VI_parallel = [[],[],[]]
for V in np.arange(0,1.6,0.05):
    # for cells in serise
    V1 = in_serise( V, T1, T2, I_L1, I_L2 )
    V2 = V - V1
    I = c1.current(V1,T1,I_L1)
    VI_serise[0].append(V1)
    VI_serise[1].append(V2)
    VI_serise[2].append(I)

    # for cells in parallel
    I = in_parallel( V, T1, T2, I_L1, I_L2 )
    I1 = c1.current( V, T1, I_L1)
    I2 = I - I1
    VI_parallel[0].append(I1)
    VI_parallel[1].append(I2)
    VI_parallel[2].append(V)

# output of serise arrangement
V1 = VI_serise[0]
V2 = VI_serise[1]
V  = [a + b for a, b in zip(V1, V2)]
I  = VI_serise[2]
power = [a * b for a, b in zip(V, I)]

plt.subplot(221)
plt.plot(V1,I,'c.')
plt.plot(V2,I,'m.')
plt.plot(V,I,'r', label='combined')

plt.subplot(223)
plt.plot(V, power, 'b')

# output power of parallel arrangement
I1 = VI_parallel[0]
I2 = VI_parallel[1]
I  = [a + b for a, b in zip(I1, I2)]
V  = VI_parallel[2]
power = [a * b for a, b in zip(V, I)]

plt.subplot(222)
plt.plot(V,I1,'c.')
plt.plot(V,I2,'m.')
plt.plot(V,I,'r', label='combined')

plt.subplot(223)
plt.plot(V, power, 'yellow')

plt.subplot(221)
plt.axis([0,1.6,0,(I_L1+I_L2)*1.1])
plt.legend(prop={'size':9})
plt.subplot(222)
plt.axis([0,1.6,0,(I_L1+I_L2)*1.1])
plt.legend(prop={'size':9})
plt.subplot(223)
plt.axis([0,1.6,0,(mpp1+mpp2)*1.1])
plt.legend(prop={'size':9})
plt.show()