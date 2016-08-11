#include "xsTypes.h"

#include "time.h"

#include "iostream"

void K_Lines_model ( double Chan,const RealArray a,double * y,int na );

using namespace std;

static double start_time=-1;
double current_time;
int current_time_tick;

extern "C" void isgribackgroundfunc (const RealArray& energy, const RealArray& parameter, int spectrum, RealArray& flux, RealArray& fluxError, const string& init)

{

// Model code:  Should resize flux RealArray to energy.size() - 1.  Do the same for

// fluxError array if calculating errors, otherwise leave it at size 0. 
//
//

  /*  cout << "ISGRI background model requested with parameters " <<  parameter.size() << endl;

    for (int ipar=0;ipar<parameter.size();ipar+=1) {
        cout << parameter[ipar] << " ";
    };
    cout << endl;*/

    current_time=clock()/CLOCKS_PER_SEC;
    if (start_time<0) {
        start_time=current_time ;
    } else {
        if ((current_time-start_time)>current_time_tick+1) {
                current_time_tick++;
                cout << "fitting tick " << current_time_tick << endl;
                for (int i=0;i<parameter.size();i++) cout << " " << parameter[i];
                cout << endl;
        };
        if (current_time_tick>1000) {
            cout << "taking too long to fit! assuming broken" << endl;
            exit(1);
        };
    };


    flux.resize(energy.size() - 1);

    for (int i=0;i<energy.size()-1;i++) {
        Real emin=energy[i];
        Real emax=energy[i+1];

        //cout << emin << " " << emax << endl;

        double f1,f2;

        K_Lines_model(emin,parameter,&f1,parameter.size()-2);
        K_Lines_model(emax,parameter,&f1,parameter.size()-2);
        
       // cout << f1 << " " << f2 << endl;

        flux[i]=f1;
    };

    //K_Lines_model_3

}



/*===============================================================*/
/*   Model the background fluorescence lines                     */
/*                                                               */
/*===============================================================*/

#include <stdio.h>
#include <math.h>

void K_Lines_model ( double Chan,const RealArray a,double * y,int na ) // not only K
{

#define  N_GAUSS 9

    /* |  W-Ka2 |   W-Ka1 |   W-Kb  |   Pb-Ka2  | Pb-Ka1  |   Pb-Kb  |*/
    double K_En[N_GAUSS]=   { 57.9817,  59.3182,  67.2443,  72.8042,   74.9694,  84.9360, 150.824, 245.395, 511. };  /* (keV) */
    static const double W_Ka_fraction[2]= {0.365,0.635};  /* W Ka2 Wka1 composition  */

    int i;
    float x[N_GAUSS];
    float Gauss[N_GAUSS];

    float N[N_GAUSS];

    float dGauss[N_GAUSS];

    float Sigma[N_GAUSS];

    //float F;
    float E;


    float Gauss_bck;
    float x_bck;



   /* if ( na!=15 )
    {
        printf ( "\n K_Lines_model_3: ERROR! The number of parameters must be 15 not %d!\n",na );
        exit ( 1 );
    }*/


    RealArray go_law=a[slice(0,3,1)];
    RealArray resolution_law=a[slice(3,3,1)];
    RealArray line_fractions=a[slice(6,8,1)];
    RealArray line_shifts=a[slice(14,9,1)];


    Real ref_E=K_En[0];
    Real ref_Eh=K_En[8];
    
    go_law[0]=go_law[0]+go_law[1]*ref_E-ref_E;

/// GO law
// Chan= a0 + a1*(E-ref_E)  // + a2*(E-ref)
// Chan= a0 + a1*(E-ref_E) + a2*a1^2*(E-ref)^2
   // E=go_law[0] + 1./go_law[1]*(Chan-(ref_E-go_law[0]))+ go_law[2]/go_law[1]*Chan*Chan; /* convert channel to energy (keV) */
   //
    
   
    if (fabs(go_law[2])<1e-7) {
        E=(Chan-go_law[0]-ref_E)/go_law[1]+ref_E;
    } else {
        E=(Chan-go_law[0]-ref_E)/go_law[1]+ref_E + go_law[2]*pow((Chan-go_law[0]-ref_E)/(ref_Eh-ref_E),2);
    };
    
    //if (fabs(E-511)<1) printf("for channel %.5lg energy is %.5lg; %.5lg %.5lg %.5lg %.5lg\n",Chan,E,go_law[0],go_law[1],go_law[2],ref_E);
    

    /*---------------------------------
     *  compute the gaussian terms
     *---------------------------------*/

    for ( i=0;i<N_GAUSS;i++ )
    {
/// Sigma law
        Sigma[i]=resolution_law[0] + resolution_law[1]*pow(Chan/120.,resolution_law[2]); //  + resolution_law[2]*resolution_law[1]*(Chan-120.)*(Chan-120.);
        //Sigma[i]=resolution_law[0] + resolution_law[1]*(Chan-120.)  + resolution_law[2]*resolution_law[1]*(Chan-120.)*(Chan-120.);
        //Sigma[i]=resolution_law[0] + resolution_law[1]*Chan + resolution_law[2]*resolution_law[1]*Chan*Chan;
        
      //  cout << "line" << K_En[i] << " shift " << line_shifts[i] << endl;

        K_En[i]+=line_shifts[i];
        

        x[i]= ( K_En[i]-E ) /Sigma[i];


        if ( fabs ( x[i] ) <12. ) /* avoid underflow when computing the exponential */
            Gauss[i]=exp ( -x[i]*x[i]/2. );
        else
            Gauss[i]=0;
    }


    N[0]=line_fractions[0]*W_Ka_fraction[0];   /* W - Ka1 */
    N[1]=line_fractions[0]*W_Ka_fraction[1];   /* W - Ka2 */

    for ( i=2;i<N_GAUSS;i++ )
        N[i]=line_fractions[i-1];

 /*   for (i=0;i<N_GAUSS;i++)
     printf("Sigma=%lg N=%lg; En=%lg\n",Sigma[i],N[i],K_En[i]);*/

    //waitkey();
    /*----------------------------------------
     *  compute the gaussan lines function
     *----------------------------------------*/

    *y=0;
    for ( i=0;i<N_GAUSS;i++ )
        *y+=N[i]*Gauss[i];

    return;


} /* end of K_Lines_model_3() */


