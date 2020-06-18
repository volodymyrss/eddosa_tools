#!/usr/bin/python

try:
    from xspec import *
    import xspec
except:
    pass

import os
import glob
import pickle
import dataanalysis as da

meroot=os.path.dirname(os.path.realpath(__file__))

def list_ht():
    ht=glob.glob(meroot+"/xspec_model/humantouch/*")
    print(meroot+"/xspec_model/humantouch/*",ht)
    return ht

#print "isgri_background",meroot+"/xspec_model"
#xspec.AllModels.lmod("isgri_background",meroot+"/xspec_model") # move it

class SpectrumRevBkg(da.DataAnalysis):
    #input_revolution=Revolution

    def main(self):
        self.spectrum=meroot+"/xspec_model/background_spec.fits"

class BestGuess(da.DataAnalysis):
    input_name="named"

    def main(self):
        self.xcm_filename=meroot+"/xspec_model/humantouch/complete_0044.xcm"

class Fit1D(da.DataAnalysis):
    input_spectrum=SpectrumRevBkg
    input_bestguess=BestGuess

    override_parameters=None

    enable_fit=True

    plot=True

    def main(self):
        self.load()
        self.select_range()
        self.define_model()

        if self.enable_fit:
            self.fit()

        self.save()
        self.save_xcm()

        if self.plot:
            self.save_plots()
        


        self.estimate_parameters()

    def estimate_parameters(self):
        pass

    def load(self):
        AllData.clear()

        fn=self.input_spectrum.get_spectrum().get_path()
        self.S=Spectrum(fn)

        print("loading",fn)
        #self.S.response=

    def select_range(self):
        pass

    def fit(self):
        AllData.show()
        Fit.show()
        Fit.query="yes"

        Fit.method=["migrad",100]
        Fit.perform()
        
        Fit.method=["leven",100,0.1,0.1,"delay"]
        Fit.perform()
    
    def save(self):
        M=self.M
        
        model_dict={}

        model_dict['chi2_reduced']=Fit.statistic/Fit.dof
        model_dict['dof']=Fit.dof
        model_dict['components']={}

        for comp_name in M.componentNames:
            print("component:",comp_name)
            model_dict['components'][comp_name]={}
            comp=getattr(M,comp_name)
            for par_name in comp.parameterNames:
                par=getattr(comp,par_name)
                print("parameter:",par_name)
                print("--- ",par.values)
                model_dict['components'][comp_name][par_name]=list(par.values)+list(par.error)

        import json
        json.dump(model_dict,open("fitresults.json","w"),sort_keys=True,indent=4, separators=(',', ': '))

        self.fitresults=model_dict
    
        self.explicit_output=['fitresults']

    def save_plots(self):
        print("plotting")

        Plot.xAxis = "keV"
        Plot.device = "/cps"

        Plot("ldata")
        os.system("ps2pdf pgplot.ps; cp pgplot.ps data.ps")
        os.system("convert pgplot.ps pgplot.pdf")
        os.system("convert pgplot.pdf data.png")
        
        Plot.device = "/cps"
        Plot("eeuf")
        #os.system("ps2pdf pgplot.ps")
        os.system("convert pgplot.ps pgplot.pdf; cp pgplot.ps unfolded.ps")
        os.system("convert pgplot.pdf unfolded.png")
        
        Plot.device = "/cps"
        Plot("emodel")
        #os.system("ps2pdf pgplot.ps")
        os.system("convert pgplot.ps pgplot.pdf; cp pgplot.ps model.ps")
        os.system("convert pgplot.pdf model.png")

    def chan_to_energy(self,chan):
        bkgfit=self.fitresults['components']['isgribackground']

        # this must be compatible with the xspec model!

        energies=[57.9817,511]

        offset=bkgfit['offset'][0]
        gain=bkgfit['gain'][0]
        gain2=bkgfit['gain2'][0]

        print("fit parameters:",offset,gain,gain2)
        offset=offset+gain*energies[0]-energies[0]

        E_func=lambda Chan: ((Chan-offset-energies[0])/gain+energies[0] + gain2*((Chan-offset-energies[0])/(energies[-1]-energies[0]))**2)
        
        print("for channel",chan,"energy",E_func(chan))
        
        return E_func(chan)


class Fit1Dnear60keV(Fit1D):
    def save(self):
        M=self.M

        allres=""

        Fit.error("maximum 100 3 4 6 9")

        print(Fit.statistic/Fit.dof)

        sig=float(M.gaussian.Sigma),M.gaussian.Sigma.error

        le=float(M.gaussian.LineE),M.gaussian.LineE.error
        le1=float(M.gaussian_3.LineE),M.gaussian_3.LineE.error
        le2=float(M.gaussian_4.LineE),M.gaussian_4.LineE.error
        le3=float(M.gaussian_5.LineE),M.gaussian_5.LineE.error

        print(sig,le,le1,le2,le3)
            
        allres+=str(dict(le=le,sig=sig,le1=le1,le2=le2,le3=le3,chi2=Fit.statistic/Fit.dof))
        open("fitresults.txt","w").write(allres.replace(",",",\n"))
        
        analysis.save(self)
    
    def select_range(self):
        #S.response="/Integral/data/ic_collection/ic_tree-20130108/ic/ibis/rsp/isgr_rmf_grp_0025.fits"
        S.response=self.input_spectrum.rmf.get_cached_path()
        self.S.ignore("1-35,161-2048")
    
    def define_model(self):
        self.M=Model("pow+gaus+gaus+gaus+gaus+gaus+gaus")
        M=self.M

        M.gaussian.LineE="58.3 0.01 57 57 60 61"

        M.gaussian_3.LineE="67. 0.01 65 65 70 70"
        M.gaussian_3.Sigma="=4"

        M.gaussian_4.LineE="75 0.01 70 70 80 81"
        M.gaussian_4.Sigma="=4"
        
        M.gaussian_5.LineE="84 0.01 80 80 100 100"
        M.gaussian_5.Sigma="=4"

        M.gaussian_6.LineE="84 0.01 80 80 100 100"
        M.gaussian_6.Sigma="7 0.01 0.1 0.1 20 20"

        M.gaussian_7.LineE="50 0.01 40 40 55 55"
        M.gaussian_7.Sigma="10 0.01 7 7 30 30"

class Fit1Dglobal(Fit1D):
    version="v3"

    def save(self):
       # Fit.error("maximum 100 3")

        Fit1D.save(self)

    gain=None
    offset=None
    
    def select_range(self):
        #S.response="/Integral/data/ic_collection/ic_tree-20130108/ic/ibis/rsp/isgr_rmf_grp_0025.fits"

        if hasattr(self.input_bestguess,'parameters'):
            pp=self.input_bestguess.parameters
            gain,offset=pp['gain'],pp['offset']
        
            if isinstance(gain,tuple):
                gain=gain[0]
            
            if isinstance(offset,tuple):
                offset=offset[0]

            if self.gain is not None: gain=self.gain
            if self.offset is not None: offset=self.offset

            e1,e2=40*gain-offset,800*gain-offset

            print("e1,e2",e1,e2)
            self.S.ignore("**-%.5lg,%.5lg-**"%(e1,e2))
        else:
            self.S.ignore("**-35.,800.-**")
    
    def define_model(self):
        print("isgri_background",meroot+"/xspec_model")

        cwd=os.getcwd()
        os.chdir(meroot+"/xspec_model")
        print("cwd:",os.getcwd())
        xspec.AllModels.lmod("isgri_background",meroot+"/xspec_model") # move it
        os.chdir(cwd)

        self.M=Model("expdec + gaussian + spexpcut*expdec + gaussian + isgribackground")


        if hasattr(self.input_bestguess,'parameters'):
            if self.override_parameters is None:
                self.override_parameters={}
            self.override_parameters=dict(list(self.override_parameters.items())+list(self.input_bestguess.parameters.items()))
        
        self.load_xcm(self.input_bestguess.xcm_filename)
                
    
    def load_xcm(self,xcm_filename):
        pars=None
        self.xcm_header=""
        for l in open(xcm_filename).readlines():
            if pars is None:
                self.xcm_header+=l
            else:
                pars.append(list(map(float,l.split())))
    
            if l.startswith("model"):
                pars=[]
        
        i=0

        for component_Name in self.M.componentNames:
            print(component_Name)
            component=getattr(self.M,component_Name)
            print(dir(component))
            for parameter_Name in component.parameterNames:
                parameter=getattr(component,parameter_Name)

                if hasattr(self,'override_parameters') and self.override_parameters is not None and parameter_Name in self.override_parameters:
                    par_o=self.override_parameters[parameter_Name]
                
                    print("overriding",par_o,"to",pars[i])
                    if isinstance(par_o,float):
                        pars[i]=[par/pars[i][0]*par_o for par in pars[i]]
                    if isinstance(par_o,tuple): # should be good size
                        pars[i]=par_o

                parameter.values=pars[i]     
                print("setting   ",component_Name,parameter_Name)
                print("setting to",pars[i])
                i+=1
        
    def save_xcm(self):
        xcm_filename="model.xcm"

        print("saving",xcm_filename)

        f=open(xcm_filename,"w")
        f.write(self.xcm_header)

        for component_Name in self.M.componentNames:
            print(component_Name)
            component=getattr(self.M,component_Name)
            print(dir(component))
            for parameter_Name in component.parameterNames:
                parameter=getattr(component,parameter_Name)

                f.write((" ".join(["%.5lg"%v for v in parameter.values]))+"\n")
        f.close()
        self.xcm=da.DataFile(xcm_filename)

class Fit1DglobalLines(Fit1Dglobal):
    version="v3"

    plot=True

    lines_table=[ 57.9817,  59.3182,  67.2443,  72.8042,   74.9694,  84.9360, 150.824, 245.395, 511.]
    
    def freeze_all(self):
        for component_Name in self.M.componentNames:
            print(component_Name)
            component=getattr(self.M,component_Name)
            print(dir(component))
            for parameter_Name in component.parameterNames:
                parameter=getattr(component,parameter_Name)
                parameter.frozen=True

    def estimate_parameters(self):

        #self.M.isgribackground.offset="0.0 -1"
        self.M.isgribackground.offset.frozen=True
        self.M.isgribackground.gain.frozen=True
        self.gain=self.M.isgribackground.gain.values[0]
        
        #self.M.isgribackground.shiftc1="0 -1 -20 -20 20 20"
        self.M.isgribackground.gain2="0.0 -1"
        #self.M.isgribackground.fractionc2="0 -1"
        #self.M.isgribackground.gain="1.0 -1"

        self.fit()

        self.freeze_all()
        

        self.M.isgribackground.shift1="0 0.001 -5 -5 5 5"
        self.M.isgribackground.shift2="=27"
        self.M.isgribackground.shift3="0 0.001 -5 -5 5 5"
        self.M.isgribackground.shift4="0 0.001 -5 -5 5 5"
        self.M.isgribackground.shift5="0 0.001 -5 -5 5 5"
        self.M.isgribackground.shift6="0 0.001 -5 -5 5 5"
        self.M.isgribackground.shiftc1="0 0.001 -20 -20 20 20"
        self.M.isgribackground.shiftc2="0 -1 -20 -20 20 20"
        self.M.isgribackground.shift511="0 0.001 -50 -50 50 50"
        self.M.isgribackground.fraction1="1 0.001"
        self.M.isgribackground.fraction2="1 0.001"
        self.M.isgribackground.fraction3="1 0.001"
        self.M.isgribackground.fraction4="1 0.001"
        self.M.isgribackground.fraction5="1 0.001"
        self.M.isgribackground.fractionc1="1 0.001"
        self.M.isgribackground.fractionc2="0 -1"
        self.M.isgribackground.fraction511="1 0.001"

        self.fit()

        Fit.error("maximum 100 27-35")

        line_tags=["1","2","3","4","5","6","c1","c2","511"]
        
        gain=self.M.isgribackground.gain.values[0]
        offset=self.M.isgribackground.offset.values[0]

        f=open("lines.txt","w")
        for lt,lref in zip(line_tags,self.lines_table):
            lsh=getattr(self.M.isgribackground,'shift'+lt).values[0]
            er=getattr(self.M.isgribackground,'shift'+lt).error
            f.write("%.5lg %.5lg %.5lg %.5lg %.5lg %.5lg\n"%(lref,lsh,er[0],er[1],gain,offset))

class Fit1Dnear511keV(Fit1D):
    def save(self):
        Fit.error("maximum 100 3")

        analysis.save(self)
    
    def select_range(self):
        #S.response="/Integral/data/ic_collection/ic_tree-20130108/ic/ibis/rsp/isgr_rmf_grp_0025.fits"
        self.S.ignore("**-400.,600.-**")
    
    def define_model(self):
        self.M=Model("pow+gaus")
        M=self.M

        M.gaussian.LineE="511 0.01 480 480 560 560"
        M.gaussian.Sigma="10 0.01 1 1 100 100"

class analysis_20keV(Fit1D):
    def save(self):
        Fit.error("maximum 100 3")

        analysis.save(self)
    
    def select_range(self):
        #S.response="/Integral/data/ic_collection/ic_tree-20130108/ic/ibis/rsp/isgr_rmf_grp_0025.fits"
        self.S.ignore("50.-**")
    
    def define_model(self):
        self.M=Model("expabs*pow+gaus")
        M=self.M
        
        M.expabs.LowECut="20. 0.01"

        M.gaussian.LineE="24 0.01 20 20 30 30"
        M.gaussian.Sigma="4. 0.01"


        
def do_dir(rootdir):
    rootdir=os.path.abspath(rootdir)
    
    A=analysis_20keV(rootdir+"/fit/20keV/allrt")
    A.specname=rootdir+"/background_spec.fits"
    A.run_procedure("root",catch=True)

    A=analysis_60keV(rootdir+"/fit/60keV/allrt")
    A.specname=rootdir+"/background_spec.fits"
    A.run_procedure("root",catch=True)
    
    A=analysis_511keV(rootdir+"/fit/511keV/allrt")
    A.specname=rootdir+"/background_spec.fits"
    A.run_procedure("root",catch=True)
    
    A40=analysis_60keV(rootdir+"/fit/60keV/lrt40")
    A40.specname=rootdir+"/background_spec_lrt40.fits"
    A40.run_procedure("root",catch=True)
    
    A40=analysis_511keV(rootdir+"/fit/511keV/lrt40")
    A40.specname=rootdir+"/background_spec_lrt40.fits"
    A40.run_procedure("root",catch=True)

    return A

#if __name__=="__main__":
#    import sys
#
#    rootdir=sys.argv[1]

#    do_dir(rootdir)

