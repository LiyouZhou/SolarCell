from matplotlib.pylab import np
from scipy import optimize
from math import exp, log


k = 1.3806488 * 10 ** -23 # boltsmann constant
q = 1.602176565 * 10 ** -19 # electron charge


class Diode(object):
    def __init__(self, n=1, I_s=10**-14 ):
        # implementation of ideal diode equation
        self.n = n # ideality factor
        self.I_s = I_s # reverse saturation current Amps

    def current( self, V_D, T ):
        global k, q
        V_T = k * T / q # thermal voltage
        I = self.I_s * ( exp( V_D/( self.n*V_T ) ) - 1 ) # ideal diode equation
        return I


class SolarCell(object):
    def __init__(self,
                 Rs     = 0.1,
                 Rsh    = 1e6,
                 n2     = 2,
                 I01    = 10**-14,
                 I02    = 10**-12,
                 bypass = False,
                 I_0bp  = 10**-12):
        self.n1     = 1
        self.n2     = n2
        self.I01    = I01
        self.I02    = I02
        self.D1     = Diode( self.n1, self.I01 )
        self.D2     = Diode( self.n2, self.I02 )
        self.Rs     = Rs
        self.Rsh    = Rsh
        self.T      = 300
        self.I_L    = 0
        self.bypass = bypass
        if self.bypass:
            self.D_bp = Diode( 1, I_0bp )

    def set_param(self, name, value):
        setattr(self, name, value)
        attr = getattr( self, name )
        if "I0" or "n" in name:
            self.D1 = Diode( self.n1, self.I01 )
            self.D2 = Diode( self.n2, self.I02 )
        print "setting %s to %e: %e"%(name, value, attr)

    def f(self, V, I): # equation for double diode model
        V_D = V + I * self.Rs
        y = self.I_L - self.D1.current(V_D,self.T) - \
                       self.D2.current(V_D,self.T) - \
            V_D/self.Rsh - I
        return y

    def df_dI(self, V, I):
        V_D = V + I * self.Rs
        y = -q*self.Rs/k/self.T*\
            ((self.D1.current(V_D,self.T)+self.D1.I_s)/self.n1
            +(self.D2.current(V_D,self.T)+self.D2.I_s)/self.n2)\
            - self.Rs/self.Rsh - 1
        return y

    def df_dV(self, V, I):
        V_D = V + I * self.Rs
        print 'V_D',V_D, 'I_dash', I, 'V',V
        print 'partial result',(self.D1.current(V_D,self.T))
        y = -q/k/self.T*\
            ((self.D1.current(V_D,self.T)+self.D1.I_s)/self.n1
            +(self.D2.current(V_D,self.T)+self.D2.I_s)/self.n2)\
            -1/self.Rsh
        print 'y=',y
        return y

    def current(self, V, T, I_L):
        def  obj_fun( I ): return  self.f(V,I)
        def dobj_fun( I ): return  self.df_dI(V,I)
        self.T = T
        self.I_L = I_L
        try: I = optimize.newton( obj_fun, -0.1, fprime=dobj_fun )
        except RuntimeError as e:
            I = float(-1)
            print e
        if self.bypass:
            return I + self.D_bp.current( -V, self.T )
        else: return float(I)

    def voltage(self, I, T, I_L):
        def  obj_fun( V ):
            if self.bypass:
                I_dash = I - self.D_bp.current( -V, self.T )
            else:
                I_dash = I
            return self.f(V,I_dash)
        def dobj_fun( V ):
            return self.df_dV(V,I)
        global k,q
        self.T = T
        self.I_L = I_L

        try:
            if self.bypass:
                V = optimize.newton( obj_fun, 0.8 )
            else:
                V = optimize.newton( obj_fun, 0.8, fprime=dobj_fun )
        except RuntimeError as e:
            V = float(-1000)
            print e
        return float(V)

    def calibrate(self, xdata, ydata):
        pass

    def objective_function(self, V, I01, I02 ):  # Rs, Rsh, n2, I01, I02):
        #if I01 < 10**-14: I01=1
        #print type(Rs)
        #self.set_param( "Rs" , Rs  )
        #self.set_param( "Rsh", Rsh )
        #self.set_param( "n2" , n2  )
        self.set_param( "I01", I01 )
        self.set_param( "I02", I02 )
        I=[]
        for v in V:
            I.append(-self.current(v, self.T, 0))
        return np.array(I)

    def residule(self, x,xdata,ydata):
        print "-"*10,"\n","x=",x

        I01=10**-x[0]
        I02=10**-x[1]
        retval = sum((np.array(ydata)-self.objective_function(xdata,I01,I02))**2)
        return retval

if __name__ == '__main__':
    c = SolarCell(bypass=1)
    c.voltage(0.1,300,0)
