#include <TCanvas.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TPad.h>
#include <TMath.h>

void calculus() {
    double p_min = 2.5; // GeV/c
    double p_max = 2.5;  // GeV/c

    const double K = 0.307075; // MeV·cm^2/mol
    const double z = 1.0;   // muon charge
    const double Z_over_A = 0.54141; // BC-404 plastic
    const double I = 6.47e-5; // MeV
    const double me = 0.510998918;  // MeV
    const double mmu = 105.6583745; // MeV
    const double rho = 1.032; // g/cm^3


     const double C = 3.1997;
    const double a = 0.16101;
    const double k = 3.2393;
    const double X_1 = 2.4855;

    // ======== BETA AND GAMMA ==============
    double E = sqrt(p_min * p_min  + mmu*0.001 * mmu*0.001 ); // total energy in MeV
        double beta = p_min / E;
        double gamma = E / (mmu*0.001);
    cout<<"beta: " << beta << endl;
    cout<<"gamma: " << gamma << endl;

    //for correction
        double bg = gamma*beta; //bg stands for beta-gamma
        double X = log10(bg);
        cout<<"beta-gamma: "<< bg <<endl;
        cout<<"X: "<< X <<endl;


        double delta = (2 * log(10.0) * X ) - C + (a * pow((X_1 - X), k));
        cout<<"delta: "<< delta <<endl;

    // ======== MASSES RATIO AND Tmax ==============
    double me_over_mmu = me / mmu;
    cout<<"masses: " <<  me_over_mmu << endl;

        double Tmax_num = 2 * me * beta * beta * gamma * gamma;

        double Tmax_den = 1 + 2 * gamma * me_over_mmu + me_over_mmu * me_over_mmu;
     double Tmax = Tmax_num / Tmax_den;
    cout<<"T-max: " << Tmax << endl;

    //======== LOGARTIM ARGUMENT ========
     double arg =(2 * me * beta * beta * gamma * gamma * Tmax) / (I * I);
       // if (arg <= 0) arg = 1e-10; // avoid log(0)
    cout<<"argument log: " << arg << endl;
    //======== LOGARITM CALCULUS ========
    double log_term = log(arg);

    //======== BRACKET =====================
     double bracket = 0.5 * log_term - beta * beta - 0.5*delta;
    cout<<"bracket: " << bracket  << endl;
    // ====== ENERGY LOSS ==================
     double dEdx = K * z * z * Z_over_A *(1/ (beta * beta)) * bracket;
     cout<<"first term: " <<  K * z * z * Z_over_A *(1/ (beta * beta))  << endl;
     double result = dEdx*rho ;

     double two_bars = result*2;
cout << "" << endl;
cout << "dE/dx = " << result << endl;
cout << "for two bars = " << two_bars << endl;
}