#include <fstream>
#include <string>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TH2.h"
#include "TStyle.h"

#include <sstream>

using namespace std;

void lastOnes_analysis()
{
    //------------- add your ROOT file -------------

    TFile *file = TFile::Open("mu_0.1_2.5_20gap_Abs_70thick.root", "READ"); 


//------------- Get the TTree of the file -------------

    TTree *tree = (TTree *)file->Get("Edep");


//------------- Define variables in order to associate them with branchess from the TTree -------------

    Double_t particleMomentum = 0.0; 

    int Hit = 0;

    Double_t bar_copy_A = 0.0;
    Double_t bar_copy_B = 0.0;

    Double_t coor_x_a = 0.0;
    Double_t coor_z_a = 0.0;
    Double_t coor_x_b = 0.0;
    Double_t coor_z_b = 0.0;

    Double_t angle_a = 0.0;
    Double_t angle_b = 0.0;


    Double_t total_energy_deposition = 0.0;
    Double_t total_dEdx = 0.0;
    int total_detected_gammas = 0.0;
    Double_t total_generated_gammas = 0.0;
   

    Double_t total_edep_layer_A = 0.0;
    Double_t total_edep_layer_B = 0.0;

    Double_t total_dEdx_layer_A = 0.0;
    Double_t total_dEdx_layer_B = 0.0;

    Double_t total_detected_photons_layer_A = 0.0;
    Double_t total_detected_photons_layer_B = 0.0;

    Double_t total_generated_photons_layer_A = 0.0;
    Double_t total_generated_photons_layer_B = 0.0;


//------------- Link the variables with branchess  -------------

    tree->SetBranchAddress("Particle_Momentum_GeV", &particleMomentum); 
    tree->SetBranchAddress("HIT_particle_passed_both_layers", &Hit);

    tree->SetBranchAddress("Copy_num_Bar_Traversed_A", &bar_copy_A);
    tree->SetBranchAddress("Copy_num_Bar_Traversed_B", &bar_copy_B);


    tree->SetBranchAddress("Position_x_Detected_On_Layer_A", &coor_x_a);
    tree->SetBranchAddress("Position_z_Detected_On_Layer_A", &coor_z_a);
    tree->SetBranchAddress("Position_x_Detected_On_Layer_B", &coor_x_b);
    tree->SetBranchAddress("Position_z_Detected_On_Layer_B", &coor_z_b);

    tree->SetBranchAddress("Angle_respect_y_A", &angle_a);                 // angle respect to Y axis   
    tree->SetBranchAddress("Angle_respect_y_B", &angle_b);


//================ TOTAL PER LAYER ==========================
    tree->SetBranchAddress("total_edep_on_A", &total_edep_layer_A);
    tree->SetBranchAddress("total_edep_on_B", &total_edep_layer_B);

    tree->SetBranchAddress("Total_dEdx_A", &total_dEdx_layer_A);
    tree->SetBranchAddress("Total_dEdx_B", &total_dEdx_layer_B);

    tree->SetBranchAddress("Total_Detected_Photons_A", &total_detected_photons_layer_A);
    tree->SetBranchAddress("Total_Detected_Photons_B", &total_detected_photons_layer_B);

    tree->SetBranchAddress("Total_Generated_Photons_A", &total_generated_photons_layer_A);
    tree->SetBranchAddress("Total_Generated_Photons_B", &total_generated_photons_layer_B);


//================ TOTAL  ==========================
    tree->SetBranchAddress("Total_Energy_Deposition", &total_energy_deposition);
    tree->SetBranchAddress("Total_dEdx", &total_dEdx);
    tree->SetBranchAddress("Total_Photons_Detected", &total_detected_gammas);
    tree->SetBranchAddress("Total_Photons_Generated", &total_generated_gammas);


//----------------  ----------------------------------------------------------------------



//--------------- TOTAL detected photons vs momentun DOTS ----------------------------------------------------------------
TCanvas* canva_17 = new TCanvas("canvas_17",                               
    "Total detected photons  vs Momentum [GeV/c] with Absorber for Muons", 
    1600, 1200);

TGraph* graph_17 = new TGraph();


graph_17->SetMarkerStyle(20);
graph_17->SetMarkerColor(kBlue);
graph_17->SetMarkerSize(0.4);


// Counter for points in the graph
Int_t nPoints = 0;


Long64_t nEntries = tree->GetEntries();

for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);


if ( total_detected_gammas > 0) {
graph_17->SetPoint(nPoints, particleMomentum, total_detected_gammas);
nPoints++;
}
}

graph_17->Draw("AP"); 

// Set titles                                                                         
graph_17->SetTitle("Total detected photons vs Momentum [GeV/c] with absorber Muons");
graph_17->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_17->GetYaxis()->SetTitle("Total Detected Photons");


canva_17->SetGrid();
canva_17->Update();
canva_17->Draw();


//                                              TOTAL detected photons vs momentun (TH2D) 


TCanvas* canvas_18 = new TCanvas("canva_18", "Total detected photons  vs Momentum [GeV/c] with Absorber for Muons", 1600, 1200);
TH2D *graph_18 = new TH2D("graph_18", "Total detected photons  vs Momentum [GeV/c] with Absorber for Muons",250, 0.1, 2.5, 250, 0, 5000);
tree->Draw("Total_Photons_Detected:Particle_Momentum_GeV>>graph_18", "Total_Photons_Detected > 0","COLZ");

graph_18->SetTitle("Total detected photons vs Momentum [GeV/c] with Absorber for Muons");
graph_18->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_18->GetYaxis()->SetTitle("Total detected photons");

canvas_18->SetGrid();
graph_18->Draw("COLZ");
canvas_18->Update();
canvas_18->Update();
canvas_18->Draw();



//------------------------------------------------- TOTAL generated photons vs momentun DOTS ----------------------------------------------------------------
TCanvas* canva_19 = new TCanvas("canvas_19",                               
    "Total genrated photons vs Momentum [GeV/c] with Absorber for Muons", 
    1600, 1200);

TGraph* graph_19 = new TGraph();


graph_19->SetMarkerStyle(20);
graph_19->SetMarkerColor(kAzure+7);
graph_19->SetMarkerSize(0.4);



for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);


if ( total_generated_gammas > 0) {
graph_19->SetPoint(nPoints, particleMomentum, total_generated_gammas);
nPoints++;
}
}

graph_19->Draw("AP"); 

// Set titles                                                                         
graph_19->SetTitle("Total generated photons vs Momentum [GeV/c] with absorber Muons");
graph_19->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_19->GetYaxis()->SetTitle("Total Generated Photons");


canva_19->SetGrid();
canva_19->Update();
canva_19->Draw();


//                                              TOTAL generated photons vs momentun (TH2D) 


TCanvas* canvas_20 = new TCanvas("canvas_20", "Total generated photons vs Momentum [GeV/c] with Absorber for Muons", 1600, 1200);
TH2D *graph_20 = new TH2D("graph_20", "Total generated photons vs Momentum [GeV/c] with Absorber for Muons",250, 0.1, 2.5, 250, 0, 220000);
tree->Draw("Total_Photons_Generated:Particle_Momentum_GeV>>graph_20", "Total_Photons_Generated > 0","COLZ");

graph_20->SetTitle("TotalGenerated photons vs Momentum [GeV/c] with Absorber for Muons");
graph_20->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_20->GetYaxis()->SetTitle("Total generatedted photons");

canvas_20->SetGrid();
graph_20->Draw("COLZ");
canvas_20->Update();
canvas_20->Update();
canvas_20->Draw();


//-------------------------------------------------|  [A] TOTAL detected photons vs momentun DOTS  |----------------------------------------------------------------
TCanvas* canva_21 = new TCanvas("canvas_21",                               
    "Total detected photons vs Momentum [GeV/c] with Absorber for Muons in Layer A", 
    1600, 1200);

TGraph* graph_21 = new TGraph();


graph_21->SetMarkerStyle(33);
graph_21->SetMarkerColor(kRed-7);
graph_21->SetMarkerSize(0.4);



for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);


if ( total_detected_photons_layer_A > 0) {
graph_21->SetPoint(nPoints, particleMomentum, total_detected_photons_layer_A);
nPoints++;
}
}

graph_21->Draw("AP"); 

// Set titles                                                                         
graph_21->SetTitle("Total detected photons vs Momentum [GeV/c] with absorber Muons in layer A");
graph_21->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_21->GetYaxis()->SetTitle("Total Detected Photons");


canva_21->SetGrid();
canva_21->Update();
canva_21->Draw();


//                                              |  [A] TOTAL detected photons vs momentun (TH2D)  |


TCanvas* canvas_22 = new TCanvas("canvas_22", "Total detected photons vs Momentum [GeV/c] with Absorber for Muons in Layer A", 1600, 1200);
TH2D *graph_22 = new TH2D("graph_22", "Total detected photons vs Momentum [GeV/c] with Absorber for Muons in layer A",250, 0.1, 2.5, 250, 0, 3500);
tree->Draw("Total_Detected_Photons_A:Particle_Momentum_GeV>>graph_22", "Total_Detected_Photons_A > 0","COLZ");

graph_22->SetTitle("Total detected photons vs Momentum [GeV/c] with Absorber for Muons in layer A");
graph_22->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_22->GetYaxis()->SetTitle("Total detected photons in layer A");

canvas_22->SetGrid();
graph_22->Draw("COLZ");
canvas_22->Update();
canvas_22->Update();
canvas_22->Draw();


//-------------------------------------------------|  [A] TOTAL generated photons vs momentun DOTS  |----------------------------------------------------------------
TCanvas* canva_23 = new TCanvas("canvas_23",                               
    "Total generated photons vs Momentum [GeV/c] with Absorber for Muons in Layer A", 
    1600, 1200);

TGraph* graph_23 = new TGraph();


graph_23->SetMarkerStyle(33);
graph_23->SetMarkerColor(kPink+8);
graph_23->SetMarkerSize(0.4);



for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);


if ( total_generated_photons_layer_A > 0) {
graph_23->SetPoint(nPoints, particleMomentum, total_generated_photons_layer_A);
nPoints++;
}
}

graph_23->Draw("AP"); 

// Set titles                                                                         
graph_23->SetTitle("Total generated photons vs Momentum [GeV/c] with absorber Muons in layer A");
graph_23->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_23->GetYaxis()->SetTitle("Total generated Photons");


canva_23->SetGrid();
canva_23->Update();
canva_23->Draw();


//                                              |  [A] TOTAL generated photons vs momentun (TH2D)  |


TCanvas* canvas_24 = new TCanvas("canvas_24", "Total generated photons vs Momentum [GeV/c] with Absorber for Muons in Layer A", 1600, 1200);
TH2D *graph_24 = new TH2D("graph_24", "Total generated photons vs Momentum [GeV/c] with Absorber for Muons in layer A",250, 0.1, 2.5, 250, 0, 140000);
tree->Draw("Total_Generated_Photons_A:Particle_Momentum_GeV>>graph_24", "Total_Generated_Photons_A > 0","COLZ");

graph_24->SetTitle("Total generated photons vs Momentum [GeV/c] with Absorber for Muons in layer A");
graph_24->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_24->GetYaxis()->SetTitle("Total generated photons in layer A");

canvas_24->SetGrid();
graph_24->Draw("COLZ");
canvas_24->Update();
canvas_24->Update();
canvas_24->Draw();



//-------------------------------------------------|  [B] TOTAL detected photons vs momentun DOTS  |----------------------------------------------------------------
TCanvas* canva_25 = new TCanvas("canvas_25",                               
    "Total detected photons vs Momentum [GeV/c] with Absorber for Muons in Layer B", 
    1600, 1200);

TGraph* graph_25 = new TGraph();


graph_25->SetMarkerStyle(33);
graph_25->SetMarkerColor(kRed-7);
graph_25->SetMarkerSize(0.4);



for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);


if ( total_detected_photons_layer_B > 0) {
graph_25->SetPoint(nPoints, particleMomentum, total_detected_photons_layer_B);
nPoints++;
}
}

graph_25->Draw("AP"); 

// Set titles                                                                         
graph_25->SetTitle("Total detected photons vs Momentum [GeV/c] with absorber Muons in layer B");
graph_25->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_25->GetYaxis()->SetTitle("Total Detected Photons");


canva_25->SetGrid();
canva_25->Update();
canva_25->Draw();


//                                              |  [B] TOTAL detected photons vs momentun (TH2D)  |


TCanvas* canvas_26 = new TCanvas("canvas_26", "Total detected photons vs Momentum [GeV/c] with Absorber for Muons in Layer B", 1600, 1200);
TH2D *graph_26 = new TH2D("graph_26", "Total detected photons vs Momentum [GeV/c] with Absorber for Muons in layer B",250, 0.1, 2.5, 250, 0, 4600);
tree->Draw("Total_Detected_Photons_B:Particle_Momentum_GeV>>graph_26", "Total_Detected_Photons_B > 0","COLZ");

graph_26->SetTitle("Total detected photons vs Momentum [GeV/c] with Absorber for Muons in layer B");
graph_26->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_26->GetYaxis()->SetTitle("Total detected photons in layer B");

canvas_26->SetGrid();
graph_26->Draw("COLZ");
canvas_26->Update();
canvas_26->Update();
canvas_26->Draw();


//-------------------------------------------------|  [B] TOTAL generated photons vs momentun DOTS  |----------------------------------------------------------------
TCanvas* canva_27 = new TCanvas("canvas_27",                               
    "Total generated photons vs Momentum [GeV/c] with Absorber for Muons in Layer B", 
    1600, 1200);

TGraph* graph_27 = new TGraph();


graph_27->SetMarkerStyle(33);
graph_27->SetMarkerColor(kBlue-4);
graph_27->SetMarkerSize(0.4);



for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);


if ( total_generated_photons_layer_B > 0) {
graph_27->SetPoint(nPoints, particleMomentum, total_generated_photons_layer_B);
nPoints++;
}
}

graph_27->Draw("AP"); 

// Set titles                                                                         
graph_27->SetTitle("Total generated photons vs Momentum [GeV/c] with absorber Muons in layer B");
graph_27->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_27->GetYaxis()->SetTitle("Total generated Photons");


canva_27->SetGrid();
canva_27->Update();
canva_27->Draw();


//                                              |  [B] TOTAL generated photons vs momentun (TH2D)  |


TCanvas* canvas_28 = new TCanvas("canvas_28", "Total generated photons vs Momentum [GeV/c] with Absorber for Muons in Layer B", 1600, 1200);
TH2D *graph_28 = new TH2D("graph_28", "Total generated photons vs Momentum [GeV/c] with Absorber for Muons in layer B",1000, 0.1, 2.5, 1000, 0, 150000);
tree->Draw("Total_Generated_Photons_B:Particle_Momentum_GeV>>graph_28", "Total_Generated_Photons_B > 0","COLZ");

graph_28->SetTitle("Total generated photons vs Momentum [GeV/c] with Absorber for Muons in layer B");
graph_28->GetXaxis()->SetTitle("Particle Momentum [GeV/c]");
graph_28->GetYaxis()->SetTitle("Total generated photons in layer B");

canvas_28->SetGrid();
graph_28->Draw("COLZ");
canvas_28->Update();
canvas_28->Update();
canvas_28->Draw();






//================================================  AVERAGE DETECTED & GENERATED PHOTONS  ========================================================
vector<Double_t> total_detected(20, 0); 
vector<Double_t> total_detected_A(20, 0);
vector<Double_t> total_detected_B(20, 0);

vector<Double_t> total_generated(20, 0); 
vector<Double_t> total_generated_A(20, 0);
vector<Double_t> total_generated_B(20, 0);

vector<int> HitCounts(20, 0);

std::vector<double> momentumEdges = {0.10, 0.22, 0.34, 0.46, 0.58, 0.7, 0.82, 0.94, 1.06, 1.18, 1.30, 1.42, 1.54, 1.66, 1.78, 1.9, 2.02, 2.14, 2.26, 2.38, 2.50};


for (Long64_t i = 0; i < nEntries; ++i) {

   tree->GetEntry(i);
   int momentumRange = -1;

   for (int j = 0; j < 20; ++j) {
       if (particleMomentum >= momentumEdges[j] && particleMomentum < momentumEdges[j+1]) {
           momentumRange = j;
           break;
           }

       }
       if(Hit == 1){
        HitCounts[momentumRange]++;
       total_detected[momentumRange]+= total_detected_gammas;
       total_generated[momentumRange]+= total_generated_gammas;
       }  

       total_detected_A[momentumRange] += total_detected_photons_layer_A;
       total_generated_A[momentumRange] += total_generated_photons_layer_A;

       total_detected_B[momentumRange] += total_detected_photons_layer_B;
       total_generated_B[momentumRange] += total_generated_photons_layer_B;
   }



vector<Double_t> yValues_avg_detected;
vector<Double_t> yValues_avg_generated;
//===========================================
vector<Double_t> xValues_II;           //===
//===========================================



cout << "=========|  PROMEDIOS FOTONES DETECTADOS Y GENERADOS TOTALES  |=========" << endl;
cout << "\n"<< endl;

for (int i = 0; i < 20; ++i) {

    Double_t momento = (i * 0.12) + 0.1 ;
        xValues_II.push_back(momento);

       double avgDetec = (HitCounts[i] > 0) ? total_detected[i] / HitCounts[i] : 0.0;
       yValues_avg_detected.push_back(avgDetec);
        
        cout << "Rango " << i + 1 << " (" << (i * 0.12) + 0.1 << " GeV - " << ((i + 1) * 0.12) + 0.1
            << " GeV):     " << HitCounts[i] << " hits,        " << total_detected[i] << "  suma de fotones detectados. " 
            << "       Promedio:     " << avgDetec << " fotones detectados " << endl; 
   }
cout << "\n"<< endl; 


for (int i = 0; i < 20; ++i) {

       double avgGener = (HitCounts[i] > 0) ? total_generated[i] / HitCounts[i] : 0.0;
       yValues_avg_generated.push_back(avgGener);
        
        cout << "Rango " << i + 1 << " (" << (i * 0.12) + 0.1 << " GeV - " << ((i + 1) * 0.12) + 0.1
            << " GeV):     " << HitCounts[i] << " hits,        " << total_generated[i] << "  suma de fotones generados. " 
            << "       Promedio:     " << avgGener << " fotones generados" << endl; 
   }
cout << "\n"<< endl; 




//------------------------------------------------------------------------------------------------------------------
vector<Double_t> yValues_avg_detected_A;
vector<Double_t> yValues_avg_generated_A;
cout << "=========|  [A] PROMEDIOS FOTONES DETECTADOS Y GENERADOS TOTALES  |=========" << endl;
cout << "\n"<< endl;

for (int i = 0; i < 20; ++i) {

       double avgDetec_A = (HitCounts[i] > 0) ? total_detected_A[i] / HitCounts[i] : 0.0;
       yValues_avg_detected_A.push_back(avgDetec_A);
        
        cout << "Rango " << i + 1 << " (" << (i * 0.12) + 0.1 << " GeV - " << ((i + 1) * 0.12) + 0.1
            << " GeV):     " << HitCounts[i] << " hits,        " << total_detected_A[i] << "  suma de fotones detectados en A. " 
            << "       Promedio:     " << avgDetec_A << " fotones detectados en A. " << endl; 
   }
cout << "\n"<< endl; 


for (int i = 0; i < 20; ++i) {

       double avgGener_A = (HitCounts[i] > 0) ? total_generated_A[i] / HitCounts[i] : 0.0;
       yValues_avg_generated_A.push_back(avgGener_A);
        
        cout << "Rango " << i + 1 << " (" << (i * 0.12) + 0.1 << " GeV - " << ((i + 1) * 0.12) + 0.1
            << " GeV):     " << HitCounts[i] << " hits,        " << total_generated_A[i] << "  suma de fotones generados en A. " 
            << "       Promedio:     " << avgGener_A << " fotones generados en A" << endl; 
   }
cout << "\n"<< endl;




//------------------------------------------------------------------------------------------------------------------
vector<Double_t> yValues_avg_detected_B;
vector<Double_t> yValues_avg_generated_B;
cout << "=========|  [B] PROMEDIOS FOTONES DETECTADOS Y GENERADOS TOTALES  |=========" << endl;
cout << "\n"<< endl;

for (int i = 0; i < 20; ++i) {

       double avgDetec_B = (HitCounts[i] > 0) ? total_detected_B[i] / HitCounts[i] : 0.0;
       yValues_avg_detected_B.push_back(avgDetec_B);
        
        cout << "Rango " << i + 1 << " (" << (i * 0.12) + 0.1 << " GeV - " << ((i + 1) * 0.12) + 0.1
            << " GeV):     " << HitCounts[i] << " hits,        " << total_detected_B[i] << "  suma de fotones detectados en B. " 
            << "       Promedio:     " << avgDetec_B << " fotones detectados en B. " << endl; 
   }
cout << "\n"<< endl; 


for (int i = 0; i < 20; ++i) {

       double avgGener_B = (HitCounts[i] > 0) ? total_generated_B[i] / HitCounts[i] : 0.0;
       yValues_avg_generated_B.push_back(avgGener_B);
        
        cout << "Rango " << i + 1 << " (" << (i * 0.12) + 0.1 << " GeV - " << ((i + 1) * 0.12) + 0.1
            << " GeV):     " << HitCounts[i] << " hits,        " << total_generated_B[i] << "  suma de fotones generados en B. " 
            << "       Promedio:     " << avgGener_B << " fotones generados en B" << endl; 
   }
cout << "\n"<< endl;



//...................................|  AVERAGE DETECTED / GENERATED PHOTONS PLOTS (TOTAL, for A and for B) |...................................


//Average total detected

TCanvas* canvas_29 = new TCanvas("canvas_29", "Average total detected photons vs Transverse Momentum for muons with absorber", 1600, 1200);
   TGraph* graph_29 = new TGraph(xValues_II.size(), &xValues_II[0], &yValues_avg_detected[0]);
   //TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
   graph_29->SetTitle("Average total detected photons vs Transverse Momentum for Muons with absorber");
   graph_29->SetMarkerStyle(22);
   graph_29->SetMarkerColor(kTeal+9); 
   graph_29->SetMarkerSize(1.5);
   graph_29->GetXaxis()->SetTitle("Momentum [GeV/c]");
 
   graph_29->GetYaxis()->SetTitle("Average detected photons");

   graph_29->Draw("AP");

   //NoAbs->SetMarkerStyle(21);
   //NoAbs->SetMarkerColor(kBlue);
   
   //NoAbs->Draw("P SAME");

   // Añadir leyenda
   TLegend* legend_29 = new TLegend(0.7, 0.7, 0.9, 0.9);
   legend_29->AddEntry(graph_29, "WITH ABSORBER", "lp");
   //legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
   legend_29->Draw();

   canvas_29->SetGrid();
   canvas_29->Update();



//Average total generated

TCanvas* canvas_30 = new TCanvas("canvas_30", "Average total generated photons vs Transverse Momentum for muons with absorber", 1600, 1200);
   TGraph* graph_30 = new TGraph(xValues_II.size(), &xValues_II[0], &yValues_avg_generated[0]);
   //TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
   graph_30->SetTitle("Average total generated photons vs Transverse Momentum for Muons with absorber");
   graph_30->SetMarkerStyle(23);
   graph_30->SetMarkerColor(kGreen-7); 
   graph_30->SetMarkerSize(1.5);
   graph_30->GetXaxis()->SetTitle("Momentum [GeV/c]");
 
   graph_30->GetYaxis()->SetTitle("Average generated photons");

   graph_30->Draw("AP");

   //NoAbs->SetMarkerStyle(21);
   //NoAbs->SetMarkerColor(kBlue);
   
   //NoAbs->Draw("P SAME");

   // Añadir leyenda
   TLegend* legend_30 = new TLegend(0.7, 0.7, 0.9, 0.9);
   legend_30->AddEntry(graph_30, "WITH ABSORBER", "lp");
   //legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
   legend_30->Draw();

   canvas_30->SetGrid();
   canvas_30->Update();



//Average detected | A |

TCanvas* canvas_31 = new TCanvas("canvas_31", "Average detected photons vs Transverse Momentum for muons with absorber in Layer A", 1600, 1200);
TGraph* graph_31 = new TGraph(xValues_II.size(), &xValues_II[0], &yValues_avg_detected_A[0]);
//TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
graph_31->SetTitle("Average detected photons vs Transverse Momentum for Muons with absorber in A");
graph_31->SetMarkerStyle(34);
graph_31->SetMarkerColor(kRed-4); 
graph_31->SetMarkerSize(1.5);
graph_31->GetXaxis()->SetTitle("Momentum [GeV/c]");

graph_31->GetYaxis()->SetTitle("Average detected photons in A");

graph_31->Draw("AP");

//NoAbs->SetMarkerStyle(21);
//NoAbs->SetMarkerColor(kBlue);

//NoAbs->Draw("P SAME");

TLegend* legend_31 = new TLegend(0.7, 0.7, 0.9, 0.9);
legend_31->AddEntry(graph_31, "WITH ABSORBER", "lp");
//legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
legend_31->Draw();

canvas_31->SetGrid();
canvas_31->Update();


//Average generated | A |

TCanvas* canvas_32 = new TCanvas("canvas_32", "Average generated photons vs Transverse Momentum for muons with absorber in Layer A", 1600, 1200);
TGraph* graph_32 = new TGraph(xValues_II.size(), &xValues_II[0], &yValues_avg_generated_A[0]);
//TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
graph_32->SetTitle("Average generated photons vs Transverse Momentum for Muons with absorber in A");
graph_32->SetMarkerStyle(43);
graph_32->SetMarkerColor(kOrange+8); 
graph_32->SetMarkerSize(1.5);
graph_32->GetXaxis()->SetTitle("Momentum [GeV/c]");

graph_32->GetYaxis()->SetTitle("Average generated photons in A");

graph_32->Draw("AP");

//NoAbs->SetMarkerStyle(21);
//NoAbs->SetMarkerColor(kBlue);

//NoAbs->Draw("P SAME");

TLegend* legend_32 = new TLegend(0.7, 0.7, 0.9, 0.9);
legend_32->AddEntry(graph_32, "WITH ABSORBER", "lp");
//legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
legend_32->Draw();

canvas_32->SetGrid();
canvas_32->Update();



//Average detected | B |

TCanvas* canvas_33 = new TCanvas("canvas_33", "Average detected photons vs Transverse Momentum for muons with absorber in Layer B", 1600, 1200);
TGraph* graph_33 = new TGraph(xValues_II.size(), &xValues_II[0], &yValues_avg_detected_B[0]);
//TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
graph_33->SetTitle("Average detected photons vs Transverse Momentum for Muons with absorber in B");
graph_33->SetMarkerStyle(33);
graph_33->SetMarkerColor(kMagenta+2); 
graph_33->SetMarkerSize(1.5);
graph_33->GetXaxis()->SetTitle("Momentum [GeV/c]");

graph_33->GetYaxis()->SetTitle("Average detected photons in B");

graph_33->Draw("AP");

//NoAbs->SetMarkerStyle(21);
//NoAbs->SetMarkerColor(kBlue);

//NoAbs->Draw("P SAME");

TLegend* legend_33 = new TLegend(0.7, 0.7, 0.9, 0.9);
legend_33->AddEntry(graph_33, "WITH ABSORBER", "lp");
//legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
legend_33->Draw();

canvas_33->SetGrid();
canvas_33->Update();


//Average generated | B |

TCanvas* canvas_34 = new TCanvas("canvas_34", "Average generated photons vs Transverse Momentum for muons with absorber in Layer B", 1600, 1200);
TGraph* graph_34 = new TGraph(xValues_II.size(), &xValues_II[0], &yValues_avg_generated_B[0]);
//TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
graph_34->SetTitle("Average generated photons vs Transverse Momentum for Muons with absorber in B");
graph_34->SetMarkerStyle(47);
graph_34->SetMarkerColor(kViolet+8); 
graph_34->SetMarkerSize(1.5);
graph_34->GetXaxis()->SetTitle("Momentum [GeV/c]");

graph_34->GetYaxis()->SetTitle("Average generated photons in B");

graph_34->Draw("AP");

//NoAbs->SetMarkerStyle(21);
//NoAbs->SetMarkerColor(kBlue);

//NoAbs->Draw("P SAME");

TLegend* legend_34 = new TLegend(0.7, 0.7, 0.9, 0.9);
legend_34->AddEntry(graph_34, "WITH ABSORBER", "lp");
//legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
legend_34->Draw();

canvas_34->SetGrid();
canvas_34->Update();


// EFFICIENCY

}



