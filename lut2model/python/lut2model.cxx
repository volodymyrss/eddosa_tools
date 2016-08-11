#include <stdio.h>
#include <math.h>
#include <sys/timeb.h> 
#include <stdlib.h>

void lut2model(int m, double *inarr, int n, double *outarr) {
    int i;

    printf("m: %i\n",m);
    printf("n: %i\n",n);

    for (i=0; i< m; i++)
        outarr[i] = inarr[i]*0+n;
}

void render_bipar(double *inarr, int in1, int in2, double *inplarr, int inpl1, int inpl2 ) {
    int i,j;

    printf("render_bipar default\n");

    struct timeb start, end;
    int diff;
    ftime(&start);

    // diagnose input
    printf("m: %i\n",in1);
    printf("n: %i\n",in2);

    if (in1!=10) {
        printf("wrong number of coloumns in the input: %i expected %i\n",in1,10);
        exit(1);
    };
    
    printf("config 0: %.5lg\n",inarr[in2*9]);
    printf("first row:");
    for (int it=0;it<in1;it++) {
        printf(" %.5lg",inarr[it*in2]);
    };
    printf("\n");

    printf("config 0: %.5lg\n",inarr[in2]);

    int min_step_q=inarr[in2*9];
    int min_step_rt=inarr[in2*9];

    double *bip=inplarr;
    int bip_size_ph=inpl1;
    int bip_size_rt=inpl2;

    int old_q_i=-1;
    int old_rt_i=-1;
    double intensity_accumulated=0;

    double intensity_threshold=1e-5;

    for (int trajectory_index=0;trajectory_index<in2;trajectory_index++) {
        double energy=inarr[trajectory_index+in2*0];
        double rt=inarr[trajectory_index+in2*1];
        double q=inarr[trajectory_index+in2*2];
        double intensity=inarr[trajectory_index+in2*3];
        double sig_e_s=inarr[trajectory_index+in2*4];
        double sig_e_l=inarr[trajectory_index+in2*5];
        double sig_t_s=inarr[trajectory_index+in2*6];
        double sig_t_l=inarr[trajectory_index+in2*7];
        double sig_cc=inarr[trajectory_index+in2*8];

    
        //printf("energy: %.5lg; rt: %.5lg; q: %.5lg\n",energy,rt,q);

        double fraction,distance;
        double sig_e,sig_t;

        int q_i=int(q);
        double frac_q=q-q_i;
        int rt_i=int(rt);
        double frac_rt=rt-rt_i;
        int q_i_d;
        int rt_i_d;

        if (old_q_i<0) old_q_i=q_i;
        if (old_rt_i<0) old_rt_i=rt_i;

        if (fabs(q_i-old_q_i)<min_step_q && fabs(rt_i-old_rt_i)<min_step_rt) { // approximation
        //if (q_i==old_q_i && rt_i==old_rt_i) {
            //printf("the same\n");
            intensity_accumulated+=intensity;
            
            continue;  
        };


        intensity=intensity+intensity_accumulated;
        
        //printf("intensity %.5lg\n",intensity);

        intensity_accumulated=0;
        q_i=old_q_i;
        rt_i=old_rt_i;
        old_q_i=int(q);
        old_rt_i=int(rt);

        double sigma_factor=5.;

        int ph_min=q_i-sig_e_s*sigma_factor;
        int ph_max=q_i+sig_e_l*sigma_factor;
        if (ph_min<0) ph_min=0;
        if (ph_max>=bip_size_ph) ph_max=bip_size_ph-1;
        int q_step=1;
        
        int rt_min=rt_i-sig_t_s*sigma_factor;
        int rt_max=rt_i+sig_t_l*sigma_factor;
        if (rt_min<0) rt_min=0;
        if (rt_max>=bip_size_rt) rt_max=bip_size_rt-1;
        int rt_step=1;

        //printf("range: %i - %i, %i - %i\n",ph_min,ph_max,rt_min,rt_max);

        for (q_i_d=ph_min;q_i_d<ph_max;q_i_d++) {
            double deviation_ph=q_i_d-q_i;
            sig_e=(deviation_ph>0)?sig_e_l:sig_e_s;

            for (rt_i_d=rt_min;rt_i_d<rt_max;rt_i_d++) {
                double deviation_rt=rt_i_d-rt_i; 
                sig_t=(deviation_rt>0)?sig_t_l:sig_t_s;
                
                
                distance=(deviation_ph/sig_e)*(deviation_ph/sig_e)/2.+(deviation_rt/sig_t)*(deviation_rt/sig_t)/2.+(deviation_ph/sig_e)*(deviation_rt/sig_t)*sin(sig_cc);
                fraction=exp(-distance)/(2*M_PI*sig_e*sqrt(sig_t_l*sig_t_s)); // normalize if want
                //fraction=exp(-distance); // / (2*M_PI*sig_e*sig_t); // normalize if want

                bip[(q_i_d)*bip_size_rt+rt_i_d]+=intensity*fraction;
            };
        };

        if (intensity<intensity_threshold) break;
    };

    ftime(&end);
    diff = (int) (1000.0 * (end.time - start.time)
        + (end.millitm - start.millitm));

    printf("Finnished in %i milliseconds. \n", diff);
};

void render_bipar_m0(double *inarr, int in1, int in2, double *inplarr, int inpl1, int inpl2 ) {
    int i,j;
    
    printf("render_bipar m0\n");

    struct timeb start, end;
    int diff;
    ftime(&start);

    printf("m: %i\n",in1);
    printf("n: %i\n",in2);

    if (in1!=10) {
        printf("wrong number of coloumns in the input: %i expected %i\n",in1,10);
        exit(1);
    };
    
    printf("config 0: %.5lg\n",inarr[in2*8]);
    printf("first row:");
    for (int it=0;it<in1;it++) {
        printf(" %.5lg",inarr[it*in2]);
    };
    printf("\n");

    printf("config 0: %.5lg\n",inarr[in2]);

    int min_step_q=inarr[in2*9];
    int min_step_rt=inarr[in2*9];

    double *bip=inplarr;
    int bip_size_ph=inpl1;
    int bip_size_rt=inpl2;

    int old_q_i=-1;
    int old_rt_i=-1;
    
    double intensity_threshold=1e-7;

    for (int trajectory_index=0;trajectory_index<in2;trajectory_index++) {
        double energy=inarr[trajectory_index+in2*0];
        double rt=inarr[trajectory_index+in2*1];
        double q=inarr[trajectory_index+in2*2];
        double intensity=inarr[trajectory_index+in2*3];
        double sig_e_s=inarr[trajectory_index+in2*4];
        double sig_e_l=inarr[trajectory_index+in2*5];
        double sig_t_s=inarr[trajectory_index+in2*6];
        double sig_t_l=inarr[trajectory_index+in2*7];
        double sig_cc=inarr[trajectory_index+in2*8];

    
        //printf("energy: %.5lg; rt: %.5lg; q: %.5lg\n",energy,rt,q);

        double fraction,distance;
        double sig_e,sig_t;

        int q_i=int(q);
        int rt_i=int(rt);
        int q_i_d;
        int rt_i_d;

        //printf("intensity %.5lg\n",intensity);

        double sigma_factor=5.;

        int ph_min=q_i-sig_e_s*sigma_factor;
        int ph_max=q_i+sig_e_l*sigma_factor;
        if (ph_min<0) ph_min=0;
        if (ph_max>=bip_size_ph) ph_max=bip_size_ph-1;
        int q_step=1;
        
        int rt_min=rt_i-sig_t_s*sigma_factor;
        int rt_max=rt_i+sig_t_l*sigma_factor;
        if (rt_min<0) rt_min=0;
        if (rt_max>=bip_size_rt) rt_max=bip_size_rt-1;
        int rt_step=1;

        //printf("range: %i - %i, %i - %i\n",ph_min,ph_max,rt_min,rt_max);

        for (q_i_d=ph_min;q_i_d<ph_max;q_i_d++) {
            double deviation_ph=q_i_d-q;
            sig_e=(deviation_ph>0)?sig_e_l:sig_e_s;

            for (rt_i_d=rt_min;rt_i_d<rt_max;rt_i_d++) {
                double deviation_rt=rt_i_d-rt;
                sig_t=(deviation_rt>0)?sig_t_l:sig_t_s;
                
                distance=(deviation_ph/sig_e)*(deviation_ph/sig_e)/2.+(deviation_rt/sig_t)*(deviation_rt/sig_t)/2.+(deviation_ph/sig_e)*(deviation_rt/sig_t)*sin(sig_cc);
                fraction=exp(-distance)/(2*M_PI*sig_e*sqrt(sig_t_l*sig_t_s)); // normalize if want

                bip[(q_i_d)*bip_size_rt+rt_i_d]+=intensity*fraction;
            };
        };

        if (intensity<intensity_threshold) break;
    };

    ftime(&end);
    diff = (int) (1000.0 * (end.time - start.time)
        + (end.millitm - start.millitm));

    printf("Finnished in %i milliseconds. \n", diff);
};

/*def get_resolutions(energy,rt_bip,q_bip,intensity):
    bip=outer(zeros(2048),zeros(256))
    dim_bip=bip.shape

    for rt,q,intens in zip(rt_bip,q_bip,intensity):
        resol1=get_resol(energy,rt) # energy?

        dim=resol1.shape

        m=round((dim[0]-1)/2)
        n=round((dim[1]-1)/2)

        print "bipar size:",dim_bip
        print "resolution size:",dim
  
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

        print "e1,e2:",e1,e2
        print "e1K,e2K:",e1K,e2K

        rt1K=0
        rt1=round(rt)-n
        if rt1 < 0:
            rt1K=-rt1
            rt1=0

        rt2K=2*n
        rt2=round(rt)+n
        
        print "rt1,rt2:",rt1,rt2
        print "rt1K,rt2K:",rt1,rt2
  
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

    print ":",rt,sig_e,sig_t_s,sig_t_l

    m=round(3.5*sig_e)  # de -3.5 sigma to + 3.5 sigma
    n=round(3.5*sig_t_l)

    print "coordinates:",m,n

#    t_psf,e_psf=get_grid()

    resol0=exp(-0.5*(((e_psf[0:2*m,0:2*n]-m)/sig_e)**2+((t_psf[0:2*m,0:2*n]-n)/sig_t_s)**2))
    resol1=exp(-0.5*(((e_psf[0:2*m,0:2*n]-m)/sig_e)**2+((t_psf[0:2*m,0:2*n]-n)/sig_t_l)**2))
    resol0[:,n:2*n]=resol1[:,n:2*n]
    resol0 = resol0/resol0.sum()

    return resol0

*/
