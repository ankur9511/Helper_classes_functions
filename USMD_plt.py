import pymbar
import numpy as np
import glob
import matplotlib.pyplot as plt
import re
def axplt(ax,x=np.array([]),val=np.array([]),llabel="",
          xlabel="",ylabel="",text="",showleg="none",
          color='C4',yticks=np.array([]),xticks=np.array([]),
          xtickparams="",ytickparams="",
          xlim=np.array([]),ylim=np.array([]),
          lw=1,marker='o',markersize=0,linestyle='-'):
    if x.size and val.size:
        ax.plot(x,val,label=llabel,c=color,linewidth=lw,marker=marker,markersize=markersize,linestyle=linestyle)
        if yticks.size:
            ax.set_yticks(yticks)
        if ylim.size:
            ax.set_ylim(ylim)
        if xlim.size:
            ax.set_xlim(xlim)
        if ylabel: 
            ax.set_ylabel(ylabel)
        if xticks.size:
            ax.set_xticks(xticks)
        ax.tick_params(direction='in',
                       bottom=True,
                       top=True,
                       right=True,
                       left=True)
        if xlabel:
            ax.set_xlabel(xlabel)
        if showleg != "none":
            leg = ax.legend(fancybox=False,frameon=False,loc='best')
            for t,c in zip(leg.get_texts(),leg.get_lines()):
                t.set_ha('left')
                t.set_color(c.get_color())
                t.set_fontweight('normal')
                for item in leg.legendHandles:
                    item.set_visible(False)


def readarrtxt(fname):
    with open(fname,'r') as f:
        x = np.array(f.read().split()).astype(float)
    return x

class US_1D:
    
    kB = 1.381e-23 * 6.022e23 / 1000.0
    beta = 1.0/kB
    
    def __init__(self,temperature=300,dname="*",
                 histarr=np.arange(-0.8,0.2,0.05),
                 eq=500,nprod=1000):
        self.temperature = temperature
        self.dirnames = []
        tdirnames = []
        tdirnames = glob.glob(dname)
        for name in tdirnames:
            if len(glob.glob(name + '/*CV*'))>=1 and len(glob.glob(name + '/*plumed*.dat'))>=1:
                self.dirnames+=[name]
        self.histedgeall = histarr
        self.chi_min = self.histedgeall[0]
        self.chi_max = self.histedgeall[-1]
        self.nbins = self.histedgeall[0:].size - 1
        self.eq = eq
        self.nprod = nprod
        self.chi_kn = []
        self.results = []
        self.bin_center_i = []

    def extraction(self,CVcol=2,depcol="none",CVname='/CV',plumedname='/plumed.dat'):
        K = len(self.dirnames)
        T_k = np.ones(K,float)*self.temperature
        N_max = self.nprod
        N_k = np.zeros([K], dtype = int) 
        K_k = np.zeros([K]) 
        chi0_k = np.zeros([K]) 
        self.chi_kn = np.zeros([K,N_max]) 
        self.u_kn = np.zeros([K,N_max]) 
        g_k = np.zeros([K])
        beta_k = 1.0/(self.kB*T_k)
        kn = -1
        for names in self.dirnames:
                CV = names + CVname
                plumedcontrol = names + plumedname
                with open(plumedcontrol) as f:
                    fread = f.read()
                    try:
                        ATmatch,kappamatch = re.findall("AT=.*\s",
                                           re.findall(".*RESTRAINT.*\n",
                                           fread)[-1])[0].split()[0:2]
                        kn = kn+1
                        chi0_k[kn] = float(ATmatch[3:])
                        K_k[kn] = float(kappamatch[6:])
                        flines=fread.split('\n')
                    except IndexError:
                        ATmatch = []
                        kappamatch = []
                    with open(CV) as f:
                        fread = f.read().split('\n')
                        lenfread = len(fread)
                        fread = iter(fread)
                        j = 0
                        jc = 0
                        flag = 1
                        n = -1
                        while j < lenfread:
                            i = next(fread)
                            j = j+1
                            if ("FIELDS" in i):
                                flag = -1
                                continue
                            if (flag < 0):
                                flag = flag+1
                                continue
                            flag = 1
                            jc = jc + 1
                            if jc > self.eq and jc < self.eq + self.nprod:
                                flagtdep = 0.0
                                try:
                                    for tdep in [depcol]:
                                        if tdep != "none":
                                            if float(i.split()[int(tdep)]) != 0.0:
                                                flagtdep = 1.0
                                    if flagtdep == 0.0 :
                                        for col in [CVcol]: 
                                            if col != "none":
                                                n = n+1
                                                self.chi_kn[kn,n] = float(i.split()[col])
                                                N_k[kn] = n+1
                                except IndexError:
                                    pass
        N_max = np.max(N_k)
        u_kln = np.zeros([K,K,N_max])
        delta = (self.chi_max - self.chi_min) / float(self.nbins)
        self.bin_center_i = np.zeros([self.nbins], np.float64)
        for i in range(self.nbins):
            self.bin_center_i[i] = self.chi_min + delta/2 + delta * i
        bin_kn = np.zeros([K,N_max], np.int32)
        for kn in range(K):
            for n in range(N_k[kn]):
                bin_kn[kn,n] = int((self.chi_kn[kn,n] - self.chi_min) / delta)
        for kn in range(K):
            for n in range(N_k[kn]):
                dchi = self.chi_kn[kn,n] - chi0_k
                u_kln[kn,:,n] = self.u_kn[kn,n] + beta_k[kn] * (K_k/2.0) * dchi**2
        mbar = pymbar.MBAR(u_kln, N_k)
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_axes([0,0,1,1])
        histogram_bin_k_nbin = np.zeros((K,self.nbins))
        for kn in range(K):
            for n in range(N_k[kn]):
                if bin_kn[kn,n] >=0 and bin_kn[kn,n] < self.nbins :
                    histogram_bin_k_nbin[kn,bin_kn[kn,n]] += 1.0
            y = histogram_bin_k_nbin[kn,:]
            y[y==0.0]=np.nan
            axplt(ax,self.bin_center_i,y,color='C1',lw=1,linestyle='-',marker='*',markersize=0,xlabel="q",ylabel="sample count")
        ax.legend(bbox_to_anchor=(1.1,1.))
        self.results = mbar.computePMF(self.u_kn, bin_kn, self.nbins, return_dict=True)
