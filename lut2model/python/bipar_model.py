from numpy import *
import scipy
import astropy.io.fits as pyfits
import time
from scipy.special import erf
import sys
import lut2model
from scipy.interpolate import RectBivariateSpline,interp1d

import os,glob

import matplotlib
matplotlib.use("Agg")
from matplotlib import pylab as p

import plot as plot_module

plot=plot_module.plot

version_attributes=[]

version_attributes.append(dict(name=None,format="l2m1"))

class auto_get_attr(object):
    def __getattr__(self,pn):
        if pn[:2]=="__" or pn[:4]=="get_":
            raise  AttributeError("no attribute: "+pn)
    
        cn='get_'+pn
        if hasattr(self,cn):
            return getattr(self,cn)()
        else:
            raise  AttributeError("no attribute: "+pn)

rtoffset=5

version_attributes.append(dict(name=None,format="dm_bal4_rt1"))
version_attributes.append(dict(name="rtbump",format="rtbump%r",default=True))
version_attributes.append(dict(name="rtoffset",format="rtoffset%.5lg",default=5))

rtbump=True

class detector(auto_get_attr):
    depth=0.2 # cm

    mu_e =870.0 # -5 C    ; cm2 s-1 V-1
    mu_t = 58.0 #-5 C      ; cm2 s-1 V-1

    tau_e = 2.75E-6 # s 0
    tau_t = 6.0E-6 # s 0

    V=120.

    def export_detector(self):
        return self.mu_e,self.mu_t,self.tau_e,self.tau_t,self.V,self.offset
    
    def import_detector(self,imp):
        self.mu_e,self.mu_t,self.tau_e,self.tau_t,self.V,self.offset=imp

    def get_d(self):
        return self.depth

    def get_l_e(self):
        return self.mu_e*self.V/self.d*self.tau_e
    
    def get_l_t(self):
        return self.mu_t*self.V/self.d*self.tau_t

    def __repr__(self):
        return "(detector: V=%(V).5lg mu_e=%(mu_e).5lg mu_t=%(mu_t).5lg tau_e=%(tau_e).5lg tau_t=%(tau_t).5lg offset=%(offset).5lg rt_offset=%(rt_offset).5lg"%dict([[a,getattr(self,a)] for a in dir(self)])

# interation physics

    def get_depth_distribution(self,energy,x):
        Absorbtion=absorption()
        peneration=Absorbtion.get(energy) #!!!!
        
        return 1./peneration*exp(-x/peneration)/(1.0-exp(-self.depth/peneration))

# electronics

    gain=1.0/0.1238  # 1.0/0.175554
    offset=0.0

    rt_offset=rtoffset

    def get_rt_bip(self,rt,energy):
        """
        model electronics measurement of the RT
        """
        if rtbump==True:
            return rt*3.0e7+15+1.0*exp(-0.5*((energy-65)/22.0)**2.0)+self.rt_offset
        elif rtbump==False:
            return rt*3.0e7+15+self.rt_offset #!!!
        elif rtbump=="step":
            b=lambda en:1.0*exp(-0.5*((en-65)/22.0)**2.0)+self.rt_offset
        
            if energy<65:
                energy=65
            return rt*3.0e7+15+self.rt_offset + b(energy)
        #return rt*3.0e7+15+6.7*exp(-0.5*((energy-65)/22.0)**2.0)+self.rt_offset

    def get_ballistic_losses_old(self,rt_bip):
        """
        model electronics measurement of the PH corrected by RT

        this is to check with GC and evolution!
        """
        r=(1.46-rt_bip/178.0) * (1.0-0.11*erf((rt_bip-121.0)/34.0)) #* (q/2048)
        r[r<0]=0
        return r
    
    def get_ballistic_losses(self,rt_bip,energy): 
        """
        model electronics measurement of the PH corrected by RT

        this is to check with GC and evolution!
        """
        r=(1.46-(rt_bip/152.0) ) * (1.0-0.11*erf((rt_bip-121.0)/34.0)) 
        #r=(1.46-(rt_bip/155.0) ) * (1.0-0.11*erf((rt_bip-121.0)/34.0)) 
        #r=(1.46-(rt_bip/178.0) ) * (1.0-0.11*erf((rt_bip-121.0)/34.0)) 
        #r=(1.46-(rt_bip/178.0)* (1-exp(-energy/200)) ) * (1.0-0.11*erf((rt_bip-121.0)/34.0)) 
        r[r<0]=0

    #    m=rt_bip<32
    #    print zip(r[m],rt_bip[m],erf((rt_bip-121.0)/34.0)[m])

        return r

    def get_q_bip(self,q,rt_bip,energy):
        """
        channel in PH is defined by gain and offset and corrected for ballistic losses
        """
        f=self.get_ballistic_losses(rt_bip,energy) 

        # bending is not strong at 60 kev. suppress it


        
        return (q*self.gain+self.offset)*f

    def model_risetime(self,td_e,td_t):
        """
        make a single rise time from the hole and electron transition times, used later to convert to the electronic channel
        """
        rt=empty_like(td_t)

        rt[:]=td_e
        rt[td_t>td_e]=td_t[td_t>td_e]

        return rt

    x0f=1.
    Ef=1.

    def get_rt_q(self,nbins=3000,extra=False):
        """
            get rise time and collected charge parametrised by the interaction depth
        """

        d=self.depth

        V=self.V
        
        mu_e=self.mu_e
        mu_t=self.mu_t
        l_e=self.l_e
        l_t=self.l_t
        
        x=linspace(0,d,nbins)
        
    # Ramo's theorem consequences
        E0=V/d
        E=ones_like(x)*E0
        x0=self.x0f*d
        m=x>x0
        E[m]=E0*self.Ef
        #E[m]=E0*(1-self.Ef*exp(-x[m]/x0))
        #E[m]=E0+(x0-x[m])*(-E0+E0*self.Ef)/(x0-d)
        
        l_t=mu_t*E*self.tau_t
        l_e=mu_e*E*self.tau_e

        td_t = x/(mu_t*E)
        td_e = (d-x)/(mu_e*E)
        rt=self.model_risetime(td_e,td_t)

    # Hecht equation
        q = l_e*(1-exp(-(d-x)/l_e))+l_t*(1-exp(-x/l_t))

        if extra:
            return x,rt,q,td_e,td_t,E
        else:
            return x,rt,q

    def get_mono_energy_bipar(self,energy): 
        """
        get detector response to the energy: peneration, transport and electronics
        """

        # transport
        x,rt,q=self.get_rt_q()

        q=q*energy
        #print "get_rt_q returns max q*energy",q.max()
        
        # apply electronics
        rt_bip=self.get_rt_bip(rt,energy)
        #print "get_rt_bin returns max rt",rt_bip.max()

        q_bip=self.get_q_bip(q,rt_bip,energy)
        #print "get_q_bip returns max q2",q_bip.max()
        
        # peneration
        intensity=self.get_depth_distribution(energy,x)

        savetxt("bip_line.txt",column_stack((x,rt_bip,q_bip,intensity)))

        return x,rt_bip,q_bip,intensity

class detector_nu(detector):
    def get_mono_energy_bipar(self,energy): 
        """
        get detector response to the energy: peneration, transport and electronics
        """

        # transport
        x,rt,q=self.get_rt_q()

        q=q*energy
        
        # apply electronics
        rt_bip=self.get_rt_bip(rt,energy)
        q_bip=self.get_q_bip(q,rt_bip,energy)
        
        # peneration
        intensity=self.get_depth_distribution(energy,x)

        savetxt("bip_line.txt",column_stack((x,rt_bip,q_bip,intensity)))

        return x,rt_bip,q_bip,intensity

class detector_nobal(detector):
    def get_q_bip(self,q,rt_bip,energy):
        f=1
        return (q*self.gain+self.offset)*f

class detector_tabulated(detector):
    mu_e_ground =1020. # cm2 s-1 V-1
    mu_t_ground = 57.0 # cm2 s-1 V-1

    def __init__(self,i_rev):
        self.i_rev=i_rev

    def get_mu_t():
        return self.mu_t_ground

    def get_mu_e(self):
        mu_e = self.mu_e_ground*(1.0-0.1*i_rev/64.0)

    def load_tables(self):
        # tab_tau = mrdfits('D:\data\S2DATA_rev041_1051\tau_e_tau_t_history_2.fits',1)
        pass

    def get_tau_e(self):
        tau_e = self.tab_tau[i_rev].tau_e
    
    def get_tau_t(self):
        tau_t = self.tab_tau[i_rev].tau_t

version_attributes.append(dict(name=None,format="murnoedge5"))

class absorption:
    def __init__(self):
        self.read()

    def read(self):
    # which one??
        self.e,self.mur,self.muer=map(array,zip(*[map(float,(lambda x:x[-3:] if len(x)>3 else x)(l.split()))  for l in open(os.environ['EDDOSA_TOOLS_ROOT']+"/lut2model/python/resources/CdTe_abs.txt").readlines()[10:]]))
        ktr=50/1000.
        self.muer[self.e<ktr]=(self.e/ktr)[self.e<ktr]**(-2)*self.muer[self.e>=ktr][0]
        self.mur[self.e<ktr]=(self.e/ktr)[self.e<ktr]**(-2)*self.mur[self.e>=ktr][0]
        savetxt("absorption.txt",column_stack((self.e,self.mur,self.muer)))


    def get(self,e):
        r=1/(6.*10**interp(log10(e/1000.),log10(self.e),log10(self.mur))) 
        return r# or mur?

version_attributes.append(dict(name=None,format="sigmasv8_rtprogr"))
#bipar_model_version+="sigmasv9"

Rcc=True
#Rcc=False

sigmalemod=""

version_attributes.append(dict(name="sigmalemod",format="sigmale%s",default=""))

def get_sigmas(energy,rt_bip,q_bip):
    sig_e=3.5*(energy/60.0)**(0.5)*(1.55+0.6*erf(((rt_bip-49))/18.0)) #; 1.8 keV # 2.5? 3.5??
    sig_t_s=1/sqrt(2)*5.5*exp(-energy/60.0/5.0)*(1.0+(rt_bip/25.0)**0.95) # 
    sig_t_s=max(sig_t_s,sig_t_s*0+1.) #> 1.0           ; short rise-time error
    #sig_t_l=1/sqrt(2)*sig_t_s*(1.0+2.0*23/energy) #*(1.+rt_bip/20.)   # long rise-time error
    sig_t_l=1/sqrt(2)*sig_t_s*(1.0+2.0*23/energy)*(1.+rt_bip/20.)   # long rise-time error # MS!
    if Rcc:
        sig_cc=q_bip*0+pi/4*(1-exp(-q_bip/60.)) * exp(-(q_bip/200.)**0.5)

        #print "Py",q_bip,rt_bip,sig_cc
    else:
        sig_cc=q_bip*0

    if sigmalemod=="ex1":
        sig_e=sig_e*(1+(20/q_bip))
    
    if sigmalemod=="ex2e":
        sig_e=sig_e/(1+(20/q_bip))
    
    if sigmalemod=="ex3":
        sig_t_l=sig_t_l*(1+(20/q_bip))
    
    if sigmalemod=="ex4":
        sig_t_l=sig_t_l*(1+(40/q_bip))
        sig_t_s=sig_t_s*(1+(20/q_bip))
    
    if sigmalemod=="ex5":
        sig_t_l=sig_t_l*(1+(40/q_bip)**2)

    if sigmalemod=="ex6":
        sig_t_l=sig_t_l*(1+5*exp(-q_bip/20))

    return sig_e,sig_t_s,sig_t_l,sig_cc

def get_resolutions(energy,rt_bip,q_bip,intensity):
    bip=outer(zeros(2048),zeros(256))
    dim_bip=bip.shape

    for rt,q,intens in zip(rt_bip,q_bip,intensity):
        resol1=get_resol(energy,rt) # energy?

        dim=resol1.shape

        m=round((dim[0]-1)/2)
        n=round((dim[1]-1)/2)

     #   print "bipar size:",dim_bip
     #   print "resolution size:",dim
  
        e1K=0
        e1=round(q)-m
  
        if e1 < 0:
            e1K=-e1
            e1=0
    
        e2K=2*m
        e2=round(q)+m

        if e2 > (dim_bip[0]-1):
            e2K=dim[0]-1-(e2-(dim_bip[0]))
            e2=dim_bip[0]-1

        if e1 > e2:
            e1=e2

      #  print "e1,e2:",e1,e2
      #  print "e1K,e2K:",e1K,e2K

        rt1K=0
        rt1=round(rt)-n
        if rt1 < 0:
            rt1K=-rt1
            rt1=0

        rt2K=2*n
        rt2=round(rt)+n
        
      #  print "rt1,rt2:",rt1,rt2
      #  print "rt1K,rt2K:",rt1,rt2
  
        if rt2 > dim_bip[1]-1:
            rt2K=2*n-(rt2-dim_bip[1])
            rt2=dim_bip[1]-1

        e1=int(e1)
        e2=int(e2)
        e1K=int(e1K)
        e2K=int(e2K)
        rt1=int(rt1)
        rt2=int(rt2)
        rt1K=int(rt1K)
        rt2K=int(rt2K)

        try:
            bip[e1:e2,rt1:rt2]=bip[e1:e2,rt1:rt2]+intens*resol1[e1K:e2K,rt1K:rt2K]/(resol1[e1K:e2K,rt1K:rt2K]).sum()
        except:
            pass

    return bip


# make more better this

K_dim=2921
e_psf=arange(K_dim)
t_psf=arange(K_dim)

Unity_t=t_psf-t_psf+1.
Unity_e=e_psf-e_psf+1.

t_psf=outer(Unity_e,t_psf)
e_psf=outer(e_psf,Unity_t)

def get_resol(energy,rt):
    sig_e,sig_t_s,sig_t_l=get_sigmas(energy,rt)

    #print ":",rt,sig_e,sig_t_s,sig_t_l

    m=round(3.5*sig_e)  # de -3.5 sigma to + 3.5 sigma
    n=round(3.5*sig_t_l)

    #print "coordinates:",m,n

#    t_psf,e_psf=get_grid()

   # resol0_t=exp(-0.5*((t_psf[0:2*n]-n)/sig_t_s)**2)
   # resol0_e=exp(-0.5*((t_psf[0:2*m]-m)/sig_e)**2)
   # resol1_t=exp(-0.5*((t_psf[0:2*n]-n)/sig_t_l)**2)
   # resol1_e=resol0_e

   # resol0=outer(resol0_t,resol0_e)
   # resol1=outer(resol1_t,resol0_e)


    resol0=exp(-0.5*(((e_psf[0:2*m,0:2*n]-m)/sig_e)**2+((t_psf[0:2*m,0:2*n]-n)/sig_t_s)**2))
    resol1=exp(-0.5*(((e_psf[0:2*m,0:2*n]-m)/sig_e)**2+((t_psf[0:2*m,0:2*n]-n)/sig_t_l)**2))
    resol0[:,n:2*n]=resol1[:,n:2*n]
    resol0 = resol0/resol0.sum()

    return resol0


def produce_perenergy(i,j):
    for energy in logspace(log10(20),3,500)[i*j:i*j+j]:
        print "energy:",energy
        print "make monoenergetic bipar"
        x,rt_bip,q_bip,intensity=detector().get_mono_energy_bipar(energy)
        print "apply resolution"
        #bip=get_resolutions(energy,rt_bip,q_bip,intensity)
        print "saving"

        #save("run/model_bipars/bipar_energy_%.5lg.npy"%energy,bip)

def produce_one_energy(energy):
    x,rt_bip,q_bip,intensity=detector().get_mono_energy_bipar(energy)  

    if False:
        p.clf()
        p.scatter(q_bip,x,c=intensity,lw=0)
        p.colorbar()
        plot("x_q.png")
        
        p.clf()
        p.scatter(q_bip,rt_bip,c=intensity,lw=0)
        p.colorbar()
        plot("rt_q.png")
            
        start= time.clock()
        bip=get_resolutions(energy,rt_bip,q_bip,intensity)
        print "done in",time.clock()-start,"seconds"

        p.clf()
        p.contourf(log10(bip[:200,:100]).transpose(),levels=linspace(-10,1,100))
        p.colorbar()
        plot("bip.png")
    
    sigmas=zip(*[get_sigmas(energy,rt) for rt in rt_bip])
    p.clf()
    p.scatter(q_bip,sigmas[0],lw=0)
    p.scatter(q_bip,sigmas[1],lw=0,color="red")
    p.scatter(q_bip,sigmas[2],lw=0,color="green")
    plot("sig_q.png")

    input=array([rt_bip*0+energy,rt_bip,q_bip,intensity,sigmas[0],sigmas[0],sigmas[1],sigmas[2],sigmas[3]])
    
    bip=zeros((2048,256))
    lut2model.render_bipar(input,bip)

    p.clf()
    p.contourf(log10(bip[:200,:100]).transpose(),levels=linspace(-10,1,100))
    p.colorbar()
    plot("bip.png")

def for_spectrum(energy,flux):
    for _energy in zip(energy,flux):
        print "energy:",_energy
        print "make monoenergetic bipar"
        x,rt_bip,q_bip,intensity=detector().get_mono_energy_bipar(_energy)
        print "apply resolution"
        #bip=get_resolutions(energy,rt_bip,q_bip,intensity)
        print "saving"

        #save("run/model_bipars/bipar_energy_%.5lg.npy"%energy,bip)

def make_lut2_simple_stack():
    rt2d,ph2d=mgrid[0:2048:1,0:256:1]
    lut2=ph2d

    savez("lut2",lut2)

    pyfits.ImageHDU(lut2).writeto("lut2_1d.fits",clobber=True)

def get_peak(det,energy):
    x,rt_bip,q_bip,intensity=det.get_mono_energy_bipar(energy)  

    #p.clf()
    #p.scatter(q_bip,x,c=intensity,lw=0)
    #p.colorbar()
    #plot("x_q.png")
    
#        p.scatter(q_bip,rt_bip,c=intensity,lw=0)

    print "energy",energy,'q',q_bip[rt_bip.argmin()]
    return q_bip[rt_bip.argmin()]

def get_go(det):
    energies=[59,511]

    #p.clf()
    positions=[]
    for energy in energies:
        x,rt_bip,q_bip,intensity=det.get_mono_energy_bipar(energy)  

        #p.clf()
        #p.scatter(q_bip,x,c=intensity,lw=0)
        #p.colorbar()
        #plot("x_q.png")
        
#        p.scatter(q_bip,rt_bip,c=intensity,lw=0)
        print "energy",energy,'q',q_bip[rt_bip.argmin()]
        
        positions.append(q_bip[rt_bip.argmin()])

    g=(positions[1]-positions[0])/(energies[1]-energies[0])
    o=positions[0]-energies[0]*g

    return o,g

def get_chan(det,energies=[57.981,150,511],render_model="default"):
    rtlim=40

    #p.clf()
    positions=[]
    for energy in energies:
        #x,rt_bip,q_bip,intensity=det.get_mono_energy_bipar(energy)  
        bip,rt_bip,q_bip,intensity=generate_bipar(det,energy,render_model="default")

        #p.clf()
        #p.scatter(q_bip,x,c=intensity,lw=0)
        #p.colorbar()
        #plot("x_q.png")

        bip_projection=bip[:,:rtlim].sum(axis=1)
        av=(bip_projection*arange(2048)).sum()/bip_projection.sum()
        
#        p.scatter(q_bip,rt_bip,c=intensity,lw=0)
        print "energy",energy,'q',q_bip[rt_bip.argmin()],av
            

 #       print "energy",energy,'q',q_bip[rt_bip.argmin()]
        
        positions.append(av/2.)

    #  E=(Chan-go_law[0]-ref_E)/go_law[1]+ref_E + go_law[2]*pow((Chan-go_law[0]-ref_E)/(ref_Eh-ref_E),2);


 #   E_func=lambda Chan,offset,gain,gain2: ((Chan-offset-energies[0])/gain+energies[0] + gain2*((Chan-offset-energies[0])/(energies[-1]-energies[0]))**2)


    return positions

def plot_parameter_sensitivity():
        
    
    def plot_dets(detectors,savein,pf):
        gos=[]
        ps=[]


        for i,det in enumerate(detectors):

            print i,len(detectors)

            gos.append(get_go(det))

       #     p.colorbar()

            mu_e,mu_t,tau_e,tau_t=det.mu_e,det.mu_t,det.tau_e,det.tau_t
            #p.title("e: %.5lg mu_e: %.5lg mu_t: %.5lg tau_e:%.5lg tau_t=%.5lg"%(energy,mu_e,mu_t,tau_e,tau_t))
            #p.xlim(10,2000)
            #p.ylim(1,200)

            try:
                os.makedirs("parameter_grid_%s"%savein)
            except:
                pass
        #    plot("parameter_grid_%s/rt_q_%s.png"%(savein,"e%.5lg_mu_e%.5lg_mu_t%.5lg_tau_e%.5lg_tau_t%.5lg"%(energy,mu_e,mu_t,tau_e,tau_t)))
                
            #start= time.clock()
            #bip=get_resolutions(energy,rt_bip,q_bip,intensity)
            #print "done in",time.clock()-start,"seconds"

            #p.clf()
            #p.contourf(log10(bip).transpose(),levels=linspace(-10,1,100))
            #p.contourf(log10(bip[:200,:100]).transpose(),levels=linspace(-10,1,100))
            #p.colorbar()
            #plot("bip.png")
            ps.append(pf(det))

            print pf(det),o,g

        p.clf()
        p.plot(ps,zip(*gos)[0]/average(zip(*gos)[0]))
        p.xlabel(savein)
        p.ylabel("offset, %")
        p.ylim([0.5,1.5])
        p.semilogx()
        plot("offset_%s.png"%savein)
        save("offset_%s.npy"%savein,(ps,zip(*gos)[0]/average(zip(*gos)[0])))

        
        p.clf()
        p.plot(ps,zip(*gos)[1]/average(zip(*gos)[1]))
        p.xlabel(savein)
        p.ylabel("gain, %")
        p.ylim([0.5,1.5])
        p.semilogx()
        plot("gain_%s.png"%savein)
        save("gain_%s.npy"%savein,(ps,zip(*gos)[0]/average(zip(*gos)[0])))

        return gos

    e_rat=logspace(-1,1,2)*detector().mu_e/detector().tau_e
    e_prod=logspace(-1,1,2)*detector().mu_e*detector().tau_e


    a,g_map_e=meshgrid(e_rat,e_prod)
    g_map_e*=0

    for ia,a in enumerate(e_rat):
        for ib,b in enumerate(e_prod):
            dr=detector()
            dr.mu_e=sqrt(a*b)
            dr.tau_e=sqrt(b/a)

            g,o=get_go(dr)

            print ia,ib,a,b,g,o

            g_map_e[ia,ib]=g

    save("g_map_e.npy",g_map_e)
    p.clf()
    p.contourf(e_rat,e_prod,g_map_e/average(g_map_e),levels=linspace(0.7,1.3))
    p.ylabel("e_rat")
    p.xlabel("e_prod")
    p.colorbar()
    plot("g_map_e.png")


    detectors=[]
    ps=[]
    for a in range(-40,20):
        dr=detector()
        dr.mu_e*=exp(a/20.)
        detectors.append(dr)
    gos=plot_dets(detectors,"mu_e",lambda x:x.mu_e)
    
    detectors=[]
    ps=[]
    for a in range(-40,20):
        dr=detector()
        dr.mu_e*=exp(a/20.)
        dr.tau_e*=exp(a/20.)
        detectors.append(dr)
    gos=plot_dets(detectors,"e_prod",lambda x:x.mu_e*x.tau_e*dr.V/dr.d)
    
    detectors=[]
    ps=[]
    for a in range(-20,20):
        dr=detector()
        dr.mu_e*=exp(a/20.)
        dr.tau_e/=exp(a/20.)
        detectors.append(dr)
    gos=plot_dets(detectors,"e_rat",lambda x:x.mu_e/x.tau_e)

    
    detectors=[]
    for a in range(-20,20):
        dr=detector()
        dr.mu_t*=exp(a/20.)
        detectors.append(dr)
    plot_dets(detectors,"mu_t",lambda x:x.mu_t)
    
    detectors=[]
    for a in range(-20,20):
        dr=detector()
        dr.tau_e*=exp(a/20.)
        detectors.append(dr)
    plot_dets(detectors,"tau_e",lambda x:x.tau_e)
    
    detectors=[]
    for a in range(-20,20):
        dr=detector()
        dr.tau_t*=exp(a/20.)
        detectors.append(dr)
    plot_dets(detectors,"tau_t",lambda x:x.tau_t)

def estimate_parameters(offset,gain):
    from scipy.optimize import minimize

    print "requested gain:",gain
    print "default gain:",get_go(detector())[1]

    def f(d_tau_e):
        d=detector()
        d.offset=offset
        d.tau_e*=d_tau_e
        o,g=get_go(d)
        return (g-gain)**2

    r=minimize(f,1,method='nelder-mead',options={'xtol': 1e-8, 'disp': True})

    print "tau_e change",r.x

    d=detector()
    d.offset=offset
    d.tau_e*=r.x
    d.tau_t*=1. # controls the angle
    
    o,g=get_go(d)
    print "fitted gain, offset",g,o
        
    d.offset=o
    o,g=get_go(d)
    print "final gain, offset",g,o
    return d

def estimate_parameters_from_peaks(peak1,peak2):
    from scipy.optimize import minimize

    eline1=59 # !!!
    eline2=511 # !!!

    d_default=detector()
    d_default.tau_t*=1
    d_default.mu_t*=0.5
    print "requested peaks:",peak1,peak2
    print "default peaks:",get_peak(d_default,eline1),get_peak(d_default,eline2)

    def f((d_tau_e,offset)):
        d=detector()
        d.import_detector(d_default.export_detector())
        d.offset=offset
        d.tau_e=d_default.tau_e*d_tau_e
        return (peak1-get_peak(d,eline1))**2+(peak2-get_peak(d,eline2))**2

    r=minimize(f,[1,0],method='nelder-mead',options={'xtol': 1e-8, 'disp': True})

    print "tau_e change",r.x

    d=detector()
    d.import_detector(d_default.export_detector())
    d.tau_e*=r.x[0]
    d.offset=r.x[1]
    
    o,g=get_go(d)
    print "fitted gain, offset",g,o
        
    return d

def make_bipar_monoenergetic(detector_model,energy,resolutionfactor=1.,resolution_step=0,render_model="default"):
    x,rt_bip,q_bip,intensity=detector_model.get_mono_energy_bipar(energy)  

    if resolutionfactor==1.:
        sigmas=zip(*[get_sigmas(energy,rt,q) for rt,q in zip(rt_bip,q_bip)])
    else:
        sigmas=[array(x)/resolutionfactor for x in zip(*[get_sigmas(energy,rt,q) for rt,q in zip(rt_bip,q_bip)])]
        
    config=zeros_like(rt_bip)
    config[0]=resolution_step # step
    input=array([rt_bip*0+energy,rt_bip,q_bip,intensity,sigmas[0],sigmas[0],sigmas[1],sigmas[2],sigmas[3],config])
    bip=zeros((2048,256))

    if render_model=="m0":
        lut2model.render_bipar_m0(input,bip)
    elif render_model=="default":
        lut2model.render_bipar(input,bip)

    return bip

version_attributes.append(dict(name=None,format="gen2"))

class response3d:
    def __init__(self,energies,npha=2048,nrt=256,grouping=0.01):
        self.npha=npha
        self.nrt=nrt
        self.energies=energies
        self.grouping=grouping

    def writeto(self,fn):
        tbhdu = pyfits.new_table(pyfits.ColDefs([
                                              pyfits.Column(name='ENERGY',
                                                            format='E',
                                                            array=self.energies)]
                                            ))
        tbhdure = pyfits.new_table(pyfits.ColDefs([
                                              pyfits.Column(name='ENERGY',
                                                            format='E',
                                                            array=self.response_energies)]
                                            ))

        imghdu=pyfits.PrimaryHDU(self.response)
        tbhdu.header['grouping']=self.grouping
        pyfits.HDUList([imghdu, tbhdu, tbhdure]).writeto(fn,clobber=True)
    
    def loadfrom(self,fn):
        f=pyfits.open(fn)
        self.response=f[0].data
        print "loading 3d response from",fn,"sum",self.response.sum()
        self.energies=f[1].data['ENERGY']
        self.response_energies=f[2].data['ENERGY']
        self.grouping=f[1].header['grouping']
        self.group_energies(init=False)

    def init(self):
        self.response=zeros((len(self.energy_groups),self.npha,self.nrt))
        self.response_energies=zeros(len(self.energy_groups))

    def group_energies(self,init=True):
        energies=self.energies
        energy_group=[energies[1]]
        energy_groups=[[energies[0]]]

        print "total energies:",len(energies)
        ie=0
        ieg=0
        for energy in energies[2:]:
            print "energy",ie,energy,
            if energy < energy_group[0]*(self.grouping+1):
                energy_group.append(energy)
                ie+=1
            else:
                energy_groups.append(energy_group)
                print "group filled",energy_group,len(energy_group),ieg
                energy_group=[energy]
                ie+=1
                ieg+=1
            print "group",ieg
            print "total in groups:",sum([len(g) for g in energy_groups])+len(energy_group)
        energy_groups.append(energy_group)

        print "total in groups:",sum([len(g) for g in energy_groups])

        self.energy_groups=energy_groups

        if init:
            self.init()

    def get_group_for_energy(self,energy):
        ie_start=0
        ie_stop=0
        ieg=0
        if energy<self.energy_groups[1][0]:
            print "energy is below minimal",self.energy_groups[1][0]
            return 1
            
        for energy_group in self.energy_groups:
         #   print "group",energy_group[0],energy_group[-1]
            if energy<energy_group[0]:
                return ieg
            if energy_group[0]<=energy<=energy_group[-1]:
                break

            ie_start=ie_stop
            ie_stop+=len(energy_group)
            ieg+=1
        if ieg>=len(self.energy_groups):
            print "energy not found in groups",energy
            print "energy is above maximal",self.energy_groups[-1][-1]
            return 
        return ieg

    def get_energy(self,energy,silent=False,disable_ai=True):
        ieg=self.get_group_for_energy(energy)

        if ieg is None:
            return

        if not silent:
            print "get energy",energy,"in group",ieg,self.energy_groups[ieg][0],self.energy_groups[ieg][-1]
            print "response:",self.response[ieg-1].sum(),self.response[ieg].sum(),self.response.sum(),self.response.sum()
        

        energy1=self.energy_groups[ieg-1][-1]
        energy2=self.energy_groups[ieg][-1]
        if energy>30 and energy2-energy1>energy2*0.02 and not disable_ai: 
            print "advanced interpolation"
            pha=arange(2048)
            rt=arange(256)

            response1=self.response[ieg-1]

            response2=self.response[ieg]

            r1=RectBivariateSpline(pha/energy1*energy,rt,response1)(pha,rt)
            r2=RectBivariateSpline(pha/energy2*energy,rt,response2)(pha,rt)
            
            r=r1  + (r2 - r1) * (energy-energy1)/(energy2-energy1)

            if False:
                pyfits.PrimaryHDU(response1).writeto("r1_raw.fits",clobber=True)
                pyfits.PrimaryHDU(r1).writeto("r1.fits",clobber=True)
                pyfits.PrimaryHDU(response2).writeto("r2_raw.fits",clobber=True)
                pyfits.PrimaryHDU(r2).writeto("r2.fits",clobber=True)
                pyfits.PrimaryHDU(r2).writeto("r.fits",clobber=True)
        else:
            print "simple interpolation"
            return self.response[ieg-1]  + (self.response[ieg] - self.response[ieg-1]) * (energy-self.energy_groups[ieg-1][-1])/(self.energy_groups[ieg][-1]-self.energy_groups[ieg-1][-1])
        return r
    
    def set_energy(self,energy,bip):
        ieg=self.get_group_for_energy(energy)
        self.response[ieg]=bip
        self.response_energies[ieg]=energy

    def get_energy_distribution_in_channel(self,pha,rt,energies):
        distribution=[]
        for energy in energies:
            r=self.get_energy(energy,silent=True)
            if r is None:
                distribution.append(0)
            else:
                distribution.append(r[pha,rt])
        return distribution # no fast!

    def get_bipar_for_model(self,model):
        bip=None
        for energy in self.energies:
            r=self.get_energy(energy,silent=True)

            if r is None:
                continue

            if bip is None:
                bip=r
            else:
                bip+=r
        
        return bip

def generate_bipar(detector_model,energy,render_model="m0",resolution_step=1):
    x,rt_bip,q_bip,intensity=detector_model.get_mono_energy_bipar(energy)  
    sigmas=zip(*[get_sigmas(energy,rt,q) for rt,q in zip(rt_bip,q_bip)])

    config=zeros_like(rt_bip)
    config[0]=resolution_step # step
    input=array([rt_bip*0+energy,rt_bip,q_bip,intensity,sigmas[0],sigmas[0],sigmas[1],sigmas[2],sigmas[3],config])

    bip=zeros((2048,256))

    if render_model=="m0":
        lut2model.render_bipar_m0(input,bip)
    elif render_model=="default":
        lut2model.render_bipar(input,bip)
    elif render_model=="m0+default":
        if energy<50:   
            lut2model.render_bipar_m0(input,bip)
        else:
            lut2model.render_bipar(input,bip)
    return bip,rt_bip,q_bip,intensity

def make_lut2_v3(detector_model=None,resolution_step=1,nenergies=401,logenergies=True,energies=None,use_resolution=True,render_model="m0",nenergies_r2=500,randomize_energies=True,source_model=lambda e:e**(-2.),group_energies=0.01,**aa):
    bip_n=zeros((256,256))
    bip_energy=zeros((256,256))
    bip_gain=zeros((256,256))

    emin,emax=(lambda x:(x['E_MIN'],x['E_MAX']))(pyfits.open(os.environ["INTEGRAL_DATA"]+"/resources/rmf_256bins.fits")['EBOUNDS'].data)
    energies=(emin+emax)/2.
    denergies=(emax-emin)/2.

    ph=arange(2048)
    
    bips=[]
    
    for energy,denergy in zip(energies,denergies)[::resolution_step]:
        print "energy",energy
        de=0.1*energy/60.
        weight=source_model(energy)

        bip,rt_bip,q_bip,intensity=generate_bipar(detector_model,energy,render_model)
        bip_pd,rt_bip_pd,q_bip_pd,intensity_pd=generate_bipar(detector_model,energy+de,render_model)

        bips.append([energy,bip])

        rt_profile=[]
        for rt in range(256)[::-1]:
            av_ph=nansum(bip[:,rt]*ph)/nansum(bip[:,rt])
            gain=(nansum(bip_pd[:,rt]*ph)/nansum(bip_pd[:,rt])-av_ph)/de

            if isnan(av_ph) or isinf(av_ph) or av_ph<=0:
                if len(rt_profile)>0:
                    av_ph=rt_profile[-1][0]
                    gain=rt_profile[-1][1]
                else:
                    continue
            
            if av_ph>ph.shape[0] or av_ph<0:
                continue
            
            rt_profile.append([av_ph,gain])

            if isnan(source_model(energy)):
                print "skipping",energy,source_model(energy)
                continue
            

            si_ph=where((av_ph>emin*2) & (av_ph<emax*2))[0] # pick single
            if len(si_ph)==0: continue
            si_ph=si_ph[0]

            bin_energy=energy+(energies[si_ph]*2-av_ph)/gain

            
            print "profile",energy,bin_energy,weight,av_ph,gain,si_ph
            
            bip_energy[si_ph,rt]+=bin_energy*weight
            bip_gain[si_ph,rt]+=gain*weight
            bip_n[si_ph,rt]+=weight

    lut21d=bip_energy/bip_n
    lut21d_gain=bip_gain/bip_n
       
    pyfits.HDUList([
            pyfits.PrimaryHDU(array([b for en,b in bips]).astype(float64)),
            pyfits.ImageHDU(array([en for en,b in bips]).astype(float64))
                    ]).writeto("response_3d.fits",clobber=True)

    pyfits.HDUList([
            pyfits.PrimaryHDU(lut21d.astype(float64)),
            pyfits.ImageHDU(lut21d_gain.astype(float64))
                    ]).writeto("lut2_1d.fits",clobber=True)
    
    return lut21d,bip_energy,bip_gain,bip_n

def make_lut2_simple_faster(writelut23d=True,detector_model=None,write3dresp=True,resolution_step=1,nenergies=401,logenergies=True,energies=None,use_resolution=True,render_model="m0",nenergies_r2=500,randomize_energies=True,source_model=lambda e:e**(-2.),group_energies=0.01):
    bip_n=zeros((2048,256))
    bip_tot=zeros((2048,256))
    bip_tot2=zeros((2048,256))
    bip_model=zeros((2048,256))

    #if detector_model is None:
    #    detector_model=detector()

    # todo: compress responses by removing empty part

    # define energies
    if energies is None:
        if logenergies:
            energies=logspace(log10(10),log10(1000),nenergies) # !!!
        else:
            energies=linspace(10,1000,nenergies) # !!
    
    savetxt("lut2_energies_orig.txt",energies)
            
    #if randomize_energies:
    #    de=energies[1:]-energies[:-1]
    #    energies[:-1]+=de*random.rand(de.shape[0])

    #print "randomized energies",energies
    #savetxt("lut2_energies.txt",energies)
    

    space3d=500
    lut23d=zeros((2048,256,space3d))

 #   response=zeros((energies.shape[0],2048,256),dtype=float32)

    r3d=response3d(energies,grouping=group_energies)
    r3d.group_energies()

    response_peak_ph=zeros((energies.shape[0],256,2),dtype=float32)
    response_nores=zeros((energies.shape[0],256,2),dtype=float32)

    pyfits.PrimaryHDU()

    
    pyfits.PrimaryHDU(generate_bipar(detector_model,45,render_model="m0",resolution_step=resolution_step)[0]).writeto("bipar_45kev.fits",clobber=True)
    pyfits.PrimaryHDU(generate_bipar(detector_model,35,render_model="m0",resolution_step=resolution_step)[0]).writeto("bipar_35kev.fits",clobber=True)
    pyfits.PrimaryHDU(generate_bipar(detector_model,60,render_model="m0",resolution_step=resolution_step)[0]).writeto("bipar_60kev.fits",clobber=True)
    pyfits.PrimaryHDU(generate_bipar(detector_model,511,render_model="m0",resolution_step=resolution_step)[0]).writeto("bipar_511kev.fits",clobber=True)


    ie=0
    ieg=0
    for energy_group in r3d.energy_groups:
    #for ie,energy in enumerate(zip(energies[:-1],energies[1:])):

        print "energy group",energy_group[0],"-",energy_group[-1],len(energy_group)

        print "generating for",energy_group[-1]

        bip,rt_bip,q_bip,intensity=generate_bipar(detector_model,energy_group[-1],render_model)

        r3d.set_energy(energy_group[-1],bip.astype(float32))
    
        for energy in energy_group:
            print ie,"/",len(energies),ieg,"/",len(r3d.energy_groups),"energy",energy,"interpolating"
            
            bip=r3d.get_energy(energy)

            #response_peak_ph[ie][:,0]=argmax(response[ie],axis=0)
            #response_peak_ph[ie][:,1]=amax(response[ie],axis=0)
            response_peak_ph[ie][:,0]=sum(outer(arange(2048),ones(256))*bip,axis=0)/sum(bip,axis=0)
            response_peak_ph[ie][:,1]=sum(bip,axis=0)
            
            response_nores[ie][:,0]=interp1d(rt_bip,q_bip,bounds_error=False)(arange(256))
            response_nores[ie][:,1]=interp1d(rt_bip,intensity,bounds_error=False)(arange(256))

            # generate lut2 2d

            #bip=response[ie,:,:]
            bip_n+=bip


            if isnan(source_model(energy)):
                print "skipping",energy,source_model(energy)
                continue
            
            print "filling",energy,source_model(energy)
            bip_tot+=bip*energy*source_model(energy)
            bip_model+=bip*source_model(energy)
            bip_tot2+=(bip*energy)**2
            print "average bipar:",average(bip)
            ie+=1
        ieg+=1

    pyfits.PrimaryHDU(response_peak_ph).writeto("response_peak_ph.fits",clobber=True)
    #pyfits.PrimaryHDU(response.astype(float16)).writeto("response_3d.fits")

    # plot lines

    def plot_lines(what=None,maxen=1000,step=0.1,tag=""):
        p.clf()
        ep=0
        for ie,(e1,e2) in enumerate(zip(energies[:-1],energies[1:])):
            if (e1-ep)/e1<step: continue
            if e1>maxen: break
            ep=e1
            plot_module.plot_line_colored(what[ie][:,0],arange(256),what[ie][:,1]/(what[:,0,1]).max(),what[ie][:,1]**0.01,cmap='autumn')
            #plot_module.plot_line_colored(what[ie][:,0],arange(256),ie*ones_like(what[ie][:,1]),1+3*what[ie][:,1]/average(what[:,0,1]))
            p.xlim([0,maxen*2.5])
            p.ylim([0,256])
            #p.plot(what[ie][:,0],arange(256))
            #p.scatter(what[ie][:,0],arange(256),c=what[ie][:,1])
        plot("curved%s.png"%tag)

    plot_lines(response_nores,1000,0.1,"_nores")
    plot_lines(response_nores,200,0.1,"_nores_till200")
    plot_lines(response_nores,50,0.02,"_nores_till50")
    plot_lines(response_peak_ph,1000,0.1,"")
    plot_lines(response_peak_ph,200,0.1,"_till200")
    plot_lines(response_peak_ph,50,0.02,"_till50")
    
    # /plot lines

    lut21d=(bip_tot/bip_model)
    savez("lut2",lut21d)
    imghdu=pyfits.PrimaryHDU(lut21d.astype(float64))
    imghdu.writeto("lut2_1d.fits",clobber=True)
    #pyfits.HDUList([imghdu, tbhdu]).writeto("lut2_1d.fits",clobber=True)
    
    if writelut23d:
        print "will generate LUT2 3D"
            
        model_bipar=r3d.get_bipar_for_model(source_model)

        pyfits.PrimaryHDU(model_bipar.astype(float64)).writeto("model_bipar.fits",clobber=True)

        spent_in_sampling=0

        lut23d=zeros((2048,256,500))
        lut23d_n=zeros((2048,256))

        for energy in r3d.energies:
            print "energy:",energy
            bip=r3d.get_energy(energy)
            bip_n=bip/model_bipar*500 # number of 


            for (pha,rt),n in ndenumerate(bip_n):
                if isnan(n) or isinf(n) or int(n)==0: continue
                n=int(n)
                n1,n2=lut23d_n[pha,rt],lut23d_n[pha,rt]+n
                #print "was",lut23d_n[pha,rt],"adding",n,"in",n1,n2
                lut23d[pha,rt][n1:n2]=energy*ones(int(n))
                lut23d_n[pha,rt]+=n
                

        if False:
            for (pha,rt),dummy in ndenumerate(zeros((2048,256))):
            #for ie,energy in enumerate(zip(energies[:-1],energies[1:])):

                #print pha,rt
                energy=0.5*(energies[:-1]+energies[1:])
                energy_distribution_flat=r3d.get_energy_distribution_in_channel(pha,rt,energy)

                energy_distribution_source=source_model(energy)*energy_distribution_flat
                energy_distribution_source/=energy_distribution_source.sum()

                from  scipy.stats import rv_discrete
                xk=energy
                pk=energy_distribution_source

                mask=~isnan(pk)

                #print xk[mask],pk[mask]

                if xk[mask].shape[0]==0:
                    print pha,rt,"only nan"
                    continue
    #
                start= time.clock()
                a=rv_discrete(name="a",values=(range(len(xk[mask])),pk[mask]))
                lut23d[pha,rt]=xk[mask][a.rvs(size=space3d)]
                spent_in_sampling+=time.clock()-start

                if rt==30:
                    print "sampling spent",spent_in_sampling,"seconds"
                    spent_in_sampling=0
                    print pha,rt,"energies",average(lut23d[pha,rt]),median(lut23d[pha,rt]),var(lut23d[pha,rt])**0.5
            
        fd=pyfits.open(os.environ['CURRENT_IC']+"/ic/ibis/mod/isgr_3dl2_mod_0001.fits")
        for i in range(500):
            fd[1].data[i]=lut23d[:,:,i].transpose()[:,::2]*30.
            #fd[1].data[i]=lut2.transpose()[:,::2]*30.
        fd.writeto("lut2_3d.fits",clobber=True)


    #print "lut2 1d summary",average(lut2)

        
    if write3dresp:
        r3d.writeto("response_3d.fits")

    save_line(r3d,59.)
    save_line(r3d,511.)
    
    #response_response(lut21d,r3d,energies,savefits=True,limitrt=None)
    #response_response(lut21d,r3d,energies,savefits=True,limitrt=40,nenergies=nenergies_r2)
    #response_response(lut21d,r3d,energies,savefits=True,limitrt=100,nenergies=nenergies_r2)
    #response_response(lut21d,r3d,energies,savefits=True,limitrt=116,nenergies=nenergies_r2)


    return lut21d,bip_tot,bip_tot2,bip_n

def save_line(r3d,energy):
    pass
    

version_attributes.append(dict(name=None,format="r2v2"))

    
def response_response(lut2=None,r3d=None,energies=None,savefits=False,limitrt=None,nenergies=500):
    lut2=pyfits.open("lut2_1d.fits")[0].data if lut2 is None else lut2
    #response=pyfits.open("response_3d.fits")[0].data if r3d is None else r3d
    lut2=lut2.astype(float64)
    #response=response.astype(float64)

    print lut2.shape

    tag=""

    lut2_s=[]


    if not isinstance(limitrt,list):
        limitrt=[limitrt]
    
    r1t=[]
    tags=[]
    for lrt in limitrt:
        if lrt is None:
            lrt=255
        lut2_s.append(copy(lut2))
        r1t.append([])
        if isinstance(lrt,tuple):
            tags.append("_%.3lg-%.3lg"%lrt)
            minlrt,lrt=lrt
            lut2_s[-1][:,lrt:]=0
            lut2_s[-1][:,:minlrt]=0
        else:
            lut2_s[-1][:,lrt:]=0
            tags.append("_%.3lg"%lrt)


    print "all tags:",tags

    i=0
    bin_edges=arange(2049)*0.5 #linspace(0,2048,nenergies)
    bin_centers=(bin_edges[:-1]+bin_edges[1:])/2.
    for i,energy in enumerate(energies):
        e_r=r3d.get_energy(energy)

        if e_r is None: continue
            
        print i,e_r.shape
        i+=2

        for i,lut2 in enumerate(lut2_s):
            print "lut2...",i
            mask=logical_and(~isnan(lut2),~isnan(e_r))
            mask=logical_and(mask,logical_and(~isinf(lut2),~isinf(e_r)))
            mask=logical_and(mask,e_r>0)
            r1=histogram(lut2[mask],weights=e_r[mask],bins=bin_edges)
            ntot=where(mask)[0].flatten().shape[0]
            #print e_r[mask][e_r[mask]>0]
            #print lut2[mask]
            #print lut2[mask]*e_r[mask]
        
            r1t[i].append(r1[0])

 #       if i==100 and False:
  #          pyfits.PrimaryHDU(array(lut2).astype(float64)).writeto("lut2.fits",clobber=True)
   #         pyfits.PrimaryHDU(array(e_r).astype(float64)).writeto("er.fits",clobber=True)
    #        pyfits.PrimaryHDU((array(e_r)*array(lut2)).astype(float64)).writeto("lut2er.fits",clobber=True)

     #       savetxt("lut2line.txt",column_stack((e_r[mask],lut2[mask])))

         #   print lut2[mask][e_r[mask].argmax()]
          #  print r1[0]

        #print "total non-nan",ntot,"average energy",average(e_r[mask])
        #print r1[0]
 #       if where(isnan(r1[0]))[0]!=[]:
#            print r1[0]
            #raise

    files=[]
    if savefits:

        for i,tag in enumerate(tags):
            imghdu=pyfits.PrimaryHDU(array(r1t[i]).astype(float64))
            tbhdu = pyfits.new_table(pyfits.ColDefs([
                                                  pyfits.Column(name='ENERGY',
                                                                format='E',
                                                                array=energies),
                                                ]
                                                ))
            tbhdu_bins = pyfits.new_table(pyfits.ColDefs([
                                                  pyfits.Column(name='ENERGY',
                                                                format='E',
                                                                array=bin_centers),
                                                ]
                                                ))

            files.append("response_2d%s.fits"%tag)
            pyfits.HDUList([imghdu, tbhdu]).writeto(files[-1],clobber=True)

    return r1t,files
    

def make_lut2_simple():
    perenergyroot="/Integral/throng/savchenk/projects/integral/LUT2/run/model_bipars"
    bip_tot=reduce(lambda x,y:x+y,[load(fn) for fn in glob.glob(perenergyroot+"/bipar_energy_*")])
    bip_n=bip_tot*0
    bip_tot2=bip_tot*0
    for fn in glob.glob("bipar_energy_*"):
        en=float(fn.replace("bipar_energy_","").replace(".npy",""))
        print en
        bip=load(fn)
        bip[bip<average(bip)]=0
        bip[bip>0]=1
        bip_tot+=bip*en
        bip_n+=bip
        bip_tot2+=(bip*en)**2

def make_lut2_simple():
    perenergyroot="/Integral/throng/savchenk/projects/integral/LUT2/run/model_bipars"
    bip_tot=reduce(lambda x,y:x+y,[load(fn) for fn in glob.glob(perenergyroot+"/bipar_energy_*")])
    for fn in glob.glob("bipar_energy_*"):
        en=float(fn.replace("bipar_energy_","").replace(".npy",""))
        print en
        bip=load(fn)
        bip[bip<average(bip)]=0
        bip[bip>0]=1
        bip_tot+=bip*en
        bip_n+=bip
        bip_tot2+=(bip*en)**2

def make_lut2_1d():
    perenergyroot="run/model_bipars"
    bip_tot=load(glob.glob(perenergyroot+"/bipar_energy_*")[0])*0
    #bip_tot=reduce(lambda x,y:x+y,[load(fn) for fn in glob.glob(perenergyroot+"/bipar_energy_*")])
    bip_n=bip_tot*0
    bip_tot2=bip_tot*0
    for fn in glob.glob(perenergyroot+"/bipar_energy_*"):
        print "fn",fn
        en=float(fn.split("/")[-1].replace("bipar_energy_","").replace(".npy",""))
        print en
        bip=load(fn)
        bip_tot+=bip*en
        bip_n+=bip
        bip_tot2+=(bip*en)**2
    
    lut2=bip_tot/bip_n

    save("lut2",lut2)

    fd=pyfits.open(os.environ['CURRENT_IC']+"/ic/ibis/mod/isgr_3dl2_mod_0001.fits")
    for i in range(500):
        fd[1].data[i]=lut2.transpose()[:,::2]*30.
    fd.writeto("lut2test.fits")

    return lut2,bip_tot,bip_tot2,bip_n

def get_version():
    version=""

    for va in version_attributes:
        if va['name'] is None:
            version+=va['format']
            continue
        
        gvalue=globals()[va['name']]
        
        if 'default' in va and gvalue==va['default']:
            continue

        version+=va['format']%gvalue

    if Rcc:
        version+=""
        #v+=".rcc"
    else:
        version+=".nrcc"
    return version
    
#produce_perenergy(int(sys.argv[1]),int(sys.argv[2]))
#make_lut2_1d()

#def inspect_lut2():
