#!/usr/bin/python

from xspec import *
import os
import pickle
import aobj

xspec.AllModels.lmod("isgri_background","/Integral/throng/savchenk/projects/integral/integral-analysis/osa10_group_13278/xspec_model/") # move it


class analysis(aobj.Persistent,aobj.Procedural):
    shortname="fit"
    response="/Integral/throng/common/energy_calibration/unitary_response.fits"

    def setup_procedures(self):
        self.add_procedure("load")
        self.add_procedure("select_range")
        self.add_procedure("define_model",deps=["load","select_range"])
        self.add_procedure("fit",deps=["define_model"])
        self.add_procedure("save",deps=["fit"])
        self.add_procedure("save_plots",deps=["fit"])
        self.add_procedure("root",deps=["save_plots","save"])

    @aobj.athostdir
    def load(self):
        AllData.clear()

        self.S=Spectrum(self.specname)
        self.S.response=self.response

    def select_range(self):
        pass


    def fit(self):
        Fit.query="yes"
        Fit.perform()
    
    @aobj.athostdir
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

        pickle.dump(model_dict,open("fitresults.pickle","w"))

    @aobj.athostdir
    def save_plots(self):
        Plot.xAxis = "keV"
        Plot.device = "/cps"

        Plot("data")
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
        
class analysis_60keV(analysis):
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

class analysis_511keV(analysis):
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

class analysis_20keV(analysis):
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

if __name__=="__main__":
    import sys

    rootdir=sys.argv[1]

    do_dir(rootdir)

