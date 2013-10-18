from matplotlib.pylab import np
from scipy import optimize
from math import exp, log
import time
import sys


k = 1.3806488 * 10 ** -23 # boltsmann constant
q = 1.602176565 * 10 ** -19 # electron charge


class Diode(object):
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


class SolarCell(object):

    def __init__(self, Rs=0.1, Rsh=1e6, n2=2, I01=10**-14, I02=10**-12):
        self.n1 = 1
        self.n2 = n2
        self.I01 = I01
        self.I02 = I02
        self.D1 = Diode( self.n1, self.I01 )
        self.D2 = Diode( self.n2, self.I02 )
        self.Rs = Rs
        self.Rsh = Rsh
        self.T = 300;
        self.I_L = 0;

    def set_param(self, name, value):
        if "I0" or "n" in name:
            self.D1 = Diode( self.n1, self.I01 )
            self.D2 = Diode( self.n2, self.I02 )
        setattr(self, name, value)
        attr = getattr( self, name )

        print "setting %s to %e: %e"%(name, value, attr)

    def f(self, I):
        V_D = self.__V - I * self.Rs
        y = self.I_L - self.D1.current(V_D,self.T) - \
                  self.D2.current(V_D,self.T) - \
                  V_D/self.Rsh - I
        return y

    def df(self, I):
        V_D = self.__V - I * self.Rs
        y = -((q*self.Rs)/(k*self.T)) * ( self.D1.current(V_D,self.T) +
                                          self.D2.current(V_D,self.T) ) \
            - self.Rs/self.Rsh - 1
        return y

    def current(self, V, T, I_L):
        self.__V = V
        self.T = T
        self.I_L = I_L
        global k,q
        try: I = optimize.newton(self.f, -0.1, fprime=self.df)
        except RuntimeError: I = float(1)
        return float(I)

    def calibrate(self, xdata, ydata):
        def objective_function( V, I01): #Rs, Rsh, n2, I01, I02):
            #if I01 < 10**-14: I01=1
            #print type(Rs)
            #self.set_param( "Rs" , Rs  )
            #self.set_param( "Rsh", Rsh )
            #self.set_param( "n2" , n2  )
            self.set_param( "I01", I01 )
            #self.set_param( "I02", I02 )
            I=[]
            for v in V:
                I.append(-self.current(v, self.T, 0))
            return np.array(I)

        def residule(I01):
            print "x=",I01
            I01=10**I01
            #I01 = I01/(10**13)
            retval = sum((np.array(ydata) - objective_function(xdata,I01))**2)
            #if retval < 1: retval=log(retval)
            print "retval", log(retval)
            return log(retval)