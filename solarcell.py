from matplotlib.pylab import np
from scipy import optimize
from math import exp, log
import time
import sys


k = 1.3806488 * 10 ** -23 # boltsmann constant
q = 1.602176565 * 10 ** -19 # electron charge


class solarcell(object):

    def __init__(self, Rs=0.1, Rsh=1e6, n2=2, I01=10**-14, I02=10**-12):
        self.n1 = 1
        self.n2 = n2
        self.I01 = I01
        self.I02 = I02
        self.D1 = diode( self.n1, self.I01 )
        self.D2 = diode( self.n2, self.I02 )
        self.Rs = Rs
        self.Rsh = Rsh
        self.T = 300;
        self.I_L = 0;

    def f(self, I):
        V_D = V - I * self.Rs
        y = self.I_L - self.D1.current(V_D,self.T) - \
                  self.D2.current(V_D,self.T) - \
                  V_D/self.Rsh - I
        return y

    def df(self, I):
        V_D = V - I * self.Rs
        y = -((q*self.Rs)/(k*self.T)) * ( self.D1.current(V_D,self.T) +
                                          self.D2.current(V_D,self.T) ) \
            - self.Rs/self.Rsh - 1
        return y

    def current(self, V, T, I_L):
        self.I_L = I_L
        global k,q
        try: I = optimize.newton(self.f, -0.1, fprime=self.df)
        except RuntimeError: I = 0
        return I

    def calibrate(self, data):
        optimize.curve_fit(f, xdata, ydata, p0=None, sigma=None, **kw)


class diode(object):
    def __init__(self, n, I_s ):
        # implementation of ideal diode equation
        self.n = n # ideality factor
        self.I_s = I_s # saturation current Amps

    def current( self, V_D, T ):
        global k
        global q
        V_T = k * T / q # thermal voltage
        I = self.I_s * ( exp( V_D/( self.n*V_T ) ) - 1 ) # ideal diode equation
        return I


if __name__ == '__main__':
    from matplotlib import pyplot as plt

    d = solarcell()
    plt.figure()
    plt.ion()
    plt.axis([0,0.8,10**-12,10])
    plt.grid('on',which='both')
    plt.xlabel('voltage/V')
    plt.ylabel('current/A')
    plt.show()

    for T in range(300,310,10):
        I_li = []
        V_li = []
        for V in np.arange(0.01,0.7,0.01):
            I = -d.current(V,T,0)
            sys.stdout.flush()
            I_li.append(I)
            V_li.append(V)
    plt.plot(V_li, I_li, hold = True )
    plt.gca().set_yscale('log')

    while(1):
        plt.draw()
        time.sleep(1)