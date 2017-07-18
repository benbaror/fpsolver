import matplotlib.pyplot as plt
from numpy import *
from constants import *
class FP(object):
    """
    """

    def __init__(self, path,label,nM,FP_NR='',ac=''):
        """
        """
        self.path = path
        self.label = label
        self.nM = nM
        Q = loadtxt(path + '/Q_last.out')
        self.Q_x = Q[:,1]
        self.Q = Q[:,5:5+4*(nM-1)+1:4]
        g = loadtxt(path + '/FP_last.out')
        self.g_x = g[:,1]
        self.g = g[:,2:2+3*(nM-1)+1:3]
        n = loadtxt(path + '/density.out')
        self.n_r = n[:,0]
        self.n = n[:,1:nM+1]
        self.Trlx = n[:,-1]
        self.rh = self.n_r[-1]
        self.Th = self.Trlx[-1]
        self.nh = sum(self.n[-1,:])
        self.nstar = self.n[-1,:].max()
        self.I0 = 8*pi**2/3./sqrt(2)*self.rh**3*(G*M0)**2\
                  /(G*M0*4e6/self.rh/pc)**(3/2.)\
                  *log(4e6)*self.nstar**2*yr*pc**-3 
        n = loadtxt(path + '/density_ub.out')
        self.I = self.Q*self.I0 # 1/Yr
        self.nub_r = n[:,0]
        self.nub = n[:,1:nM+1]
        n = loadtxt(path + '/density_b.out')
        self.nb_r = n[:,0]
        self.nb = n[:,1:nM+1]
        N = loadtxt(path + '/number.out')            
        self.N_r = N[:,0]
        self.N = N[:,1:nM+1]
        self.get_constants()
        self.convergence = self.get_convergence()
        self.a_crit()
        if ac != '':
            self.ac = ac
        if FP_NR:
            self.FP_NR = FP_NR
        if 'FP_NR' in self.__dict__.keys():
            self.Qgw = array([interp(1/self.ac,self.Q_x,self.Q[:,i]) for i in xrange(len(self.Q[0,:]))])
            self.Igw = self.Qgw*self.I0
            self.Qgw_NR = array([interp(1/self.FP_NR.ac_NR,self.FP_NR.Q_x,
                                        self.FP_NR.Q[:,i]) \
                                     for i in xrange(len(self.FP_NR.Q[0,:]))])
            self.Igw_NR = self.Qgw_NR*self.FP_NR.I0
            b = where(self.Q_x > 0.5)
            if (nM == 4):
                self.xd = zeros(4)
                for i in range(4):
                    self.xd[i] = self.Q_x[b][argmin(abs(self.Q[b,i]-self.FP_NR.Q[b,i])/self.Q[b,i])]
            elif (nM == 1):
                self.xd = self.Q_x[b][argmin(abs(self.Q[b]-self.FP_NR.Q[b])/self.Q[b])]
        self.get_times()
        if (nM == 4):
            self.labels = ['WD','MS','NS','BH']
            self.colors = ['r','g','b','m']
            self.markers = ['','','','']
            self.linestyles = ['-','-','-','-']
        elif (nM == 1):
            self.labels = label
            self.colors = 'r'
            self.markers = ''
            self.linestyles = '-'

    def get_constants(self):
        f = open(self.path+'/FP.inp')
        lines = f.readlines()
        f.close()
        self.AGR = 1.0
        for line in lines:
            if line[1:4] in ['ARR','AGR','ANR','AEP']:
                exec('self.'+line[1:-2])


    def get_convergence(self):
        f = open(self.path+'/FP.out')
        lines = f.readlines()
        f.close()
        for line in lines:
            if 'Global' in line:
                return line
        pass
    def a_crit(self):
        Mbh = 4e6
        Ms = 1.0
        Q = Mbh/Ms
        Nh = 4e6
        ANR = self.ANR/log(Q)
        AGR = self.AGR
        ARR = self.ARR
        acusp = 1.75
        p = acusp - 3/2.
        s0 = ANR*85*pi/3./2.**9*Q/Nh
        n_factor = 8*sqrt(2)/31.366
        rh = 2.0
        if rh == 0:
            rh = 2*sqrt(Mbh/3e6) ### Debug this
        k = 2/(3.-2.*p)
        j02 = 1/5.*7.*(1+sqrt(161))/16.
        self.ac_NR = (s0*n_factor)**k*rh
        self.ac_e = (1+8/3.*AGR*ANR/ARR*(j02-1)/log(j02))**(-k)*self.ac_NR #Without ecc dep
        self.ac = (1+8.8*pi**2/12.*AGR*ANR/ARR*(j02-1)/log(j02))**(-k)*self.ac_NR
        pass
        


    def Na(self):
        x = self.g_x
        g = self.g
        n = self.nb
        y = log(x)
        g = g.flatten()
        Ne = g*x**(-5/2.)*x
        d = [(Ne[i]+Ne[i+1])/2.*(-y[i-1]+y[i]) for i in xrange(len(x)-2)]
        I = array([sum(d[i:]) for i in xrange(len(x))])
        return pi**(3/2.)*n[-1]*2**(3)/2.*I

    def Plot_q(self,color,marker,linestyle,label):
        plt.loglog(self.Q_x,self.Q/self.FP_NR.Q,
                   color=color,
                   marker=marker,
                   linestyle=linestyle,
                   label=label)

        leg = plt.legend(loc = 3)
        for t in leg.get_texts():
            t.set_fontsize('small')
        plt.xlim([2,1e3])
        #plt.ylim([1e-3,10])
        plt.xlabel('$X$')
        plt.ylabel('$Q$')

    def Plot_Q(self,index=''):
        if (self.nM == 1):
            plt.loglog(self.Q_x,self.Q,
                       label=self.labels,
                       color=self.colors,
                       marker=self.markers,
                       linestyle=self.linestyles)
            if 'FP_NR' in self.__dict__.keys():
                plt.loglog(self.Q_x,self.Q/self.FP_NR.Q,
                           color=self.colors,
                           marker=self.markers,
                           linestyle=self.linestyles)
                plt.loglog(self.Q_x,self.FP_NR.Q,
                           color=self.FP_NR.colors,
                           marker=self.FP_NR.markers,
                           linestyle=self.linestyles)

        else:
            if index == '':
                List = xrange(self.nM)
            else:
                List = [index]
            for i in List:
                plt.loglog(self.Q_x,self.Q[:,i],
                           label=self.labels[i],
                           color=self.colors[i],
                           marker=self.markers[i],
                           linestyle=self.linestyles[i])
                if 'FP_NR' in self.__dict__.keys():
                    plt.loglog(self.Q_x,self.Q[:,i]/self.FP_NR.Q[:,i],
                               color=self.colors[i],
                               marker=self.markers[i],
                               linestyle=self.linestyles[i])
                    plt.loglog(self.Q_x,self.FP_NR.Q[:,i],
                               color=self.colors[i],
                               marker=self.markers[i],
                               linestyle=self.linestyles[i])


        leg = plt.legend(loc = 3)
        for t in leg.get_texts():
            t.set_fontsize('small')
        plt.xlim([2,1e3])
        #plt.ylim([1e-3,10])
        plt.xlabel('$X$')
        plt.ylabel('$Q$')

    def Plot_g(self):
        if (self.nM == 1):
            plt.loglog(self.g_x,self.g,
                       label=self.labels,
                       color=self.colors,
                       marker=self.markers,
                       linestyle=self.linestyles)
        else:
            for i in xrange(self.nM):
                plt.loglog(self.g_x,self.g[:,i],
                           label=self.labels[i],
                           color=self.colors[i],
                           marker=self.markers[i],
                           linestyle=self.linestyles[i])
            
        leg = plt.legend(loc = 2)
        for t in leg.get_texts():
            t.set_fontsize('small')
        plt.xlim([1e-1,2e4])
        #plt.ylim([1,30])
        plt.xlabel('$X$')
        plt.ylabel('$g$')
################################################
    def Plot_n(self):
        if (self.nM == 1):
            plt.loglog(self.n_r,self.n,
                       label=self.labels,
                       color=self.colors,
                       marker=self.markers,
                       linestyle=self.linestyles)
        else:
            for i in xrange(self.nM):
                plt.loglog(self.n_r,self.n[:,i],
                           label=self.labels[i],
                           color=self.colors[i],
                           marker=self.markers[i],
                           linestyle=self.linestyles[i])

        leg = plt.legend(loc = 3)
        for t in leg.get_texts():
            t.set_fontsize('small')
        plt.xlim([1e-4,2])
        #plt.ylim([2e4,2e12])
        plt.xlabel('$r[pc]$')
        plt.ylabel('$n[pc^{-3}]$')

################################################
    def Plot_N(self):
        if (self.nM == 1):
            plt.loglog(self.N_r,self.N,
                       label=self.labels,
                       color=self.colors,
                       marker=self.markers,
                       linestyle=self.linestyles)
        else:
            for i in xrange(self.nM):
                plt.loglog(self.N_r,self.N[:,i],
                           label=self.labels[i],
                           color=self.colors[i],
            marker=self.markers[i],
            linestyle=self.linestyles[i])

        leg = plt.legend(loc = 2)
        for t in leg.get_texts():
            t.set_fontsize('small')
        plt.xlim([1e-4,2])
        #plt.ylim([1,4e6])
        plt.xlabel('$r[pc]$')
        plt.ylabel('$N$')

    def get_times(self):
        """
        
        Arguments:
        - `g`:
        """
        rh2rg = 1.0e7
        x = self.g_x
        g = self.g
        M = array([6.000E-01,1.000E+00,1.400E+00,1.000E+01])
        Cm = array([1.000E-01,1.000E+00,1.000E-02,1.000E-03])
        if self.nM == 1:
            M = array([1.000E+00])
            Cm = array([1.000E+00])
        trlx = 1/sum(g*M**2,1)
        tm = sqrt(32)*sum(Cm*M)/sum(g*M,1)*x**(3.0/2.0)
        Jc2 = rh2rg/32./x
        tGR_inv = 3./8.*log(Jc2)/Jc2
        trr = (1/tm+tGR_inv)*trlx
        self.tm = tm
        self.tGR_inv = tGR_inv
        self.trlx = trlx
        self.trr  = trr 



