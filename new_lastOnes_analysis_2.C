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

void new_lastOnes_analysis_2()
{
    //------------- add your ROOT file -------------

    TFile *file = TFile::Open("mu_0.1_2.5_20gap_Abs_70thick.root", "READ"); 
   // TFile *file = TFile::Open("Output0_t0.root", "READ"); 


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



//=============================================================== Bethe-Bloch Curvature ===============================================================

/* Here I plot the total energy loss vs the particle momentum in terms of MeV because it was multiplied by 2 cm
 and density of BC-404. Here I defined two plots, one with dots and the other with density.*/



//                                        TOTAL energy depostion vs p_t dots


TCanvas* canva_1 = new TCanvas("canva_1",                               
    "Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons", 
    1600, 1200);

TGraph* graph_1 = new TGraph();


graph_1->SetMarkerStyle(20);
graph_1->SetMarkerColor(kBlue);
graph_1->SetMarkerSize(0.4);


// Counter for points in the graph
Int_t nPoints = 0;

//===========================================
Long64_t nEntries = tree->GetEntries(); //===
//===========================================

for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);


if ( total_energy_deposition > 0) {
graph_1->SetPoint(nPoints, particleMomentum, total_energy_deposition);
nPoints++;
}
}

graph_1->Draw("AP"); 

// Set titles                                                                         
graph_1->SetTitle("Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with absorber Muons");
graph_1->GetXaxis()->SetTitle("Particle Momentum [GeV]");
graph_1->GetYaxis()->SetTitle("Total Energy Deposition [MeV]");


canva_1->SetGrid();
canva_1->Update();
canva_1->Draw();



//                                              p_t VS Total Energy Loss (TH2D) 


TCanvas* canva_bethe_mu_si = new TCanvas("canva_2", "Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons", 1600, 1200);
TH2D *edepmuon_abs = new TH2D("edepmuon_abs", "Muons Transverse Momentun [GeV] vs Total Deposited Energy [MeV] for Muons with absorber ",250, 0.1, 2.5, 250, 0, 20);
tree->Draw("Total_Energy_Deposition:Particle_Momentum_GeV>>edepmuon_abs", "Total_Energy_Deposition > 0","COLZ");

edepmuon_abs->SetTitle("Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons");
edepmuon_abs->GetXaxis()->SetTitle("Particle Momentum [GeV]");
edepmuon_abs->GetYaxis()->SetTitle("Total Energy Deposition [MeV]");

canva_bethe_mu_si->SetGrid();
edepmuon_abs->Draw("COLZ");
canva_bethe_mu_si->Update();
canva_bethe_mu_si->Update();
canva_bethe_mu_si->Draw();



//------------------------------------------------ p_t vs Energy Loss layer A dots ------------------------------------------------

TCanvas* canva_3 = new TCanvas("canva_3",                               
    "Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons in Layer A", 
    1600, 1200);

TGraph* graph_3 = new TGraph();


graph_3->SetMarkerStyle(20);
graph_3->SetMarkerColor(kRed);
graph_3->SetMarkerSize(0.4);


for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);

if ( total_edep_layer_A >= 0) {
graph_3->SetPoint(nPoints, particleMomentum, total_edep_layer_A);
nPoints++;
}
}

graph_3->Draw("AP"); 

// Set titles                                                                         
graph_3->SetTitle("Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with absorber Muons for layer A");
graph_3->GetXaxis()->SetTitle("Particle Momentum [GeV]");
graph_3->GetYaxis()->SetTitle("Total Energy Deposition [MeV]");


canva_3->SetGrid();
canva_3->Update();
canva_3->Draw();



//                                           p_t VS  Energy Loss in layer A (TH2D) 


TCanvas* canva_bethe_mu_si_A = new TCanvas("canva_4", "Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons for layer A", 1600, 1200);
TH2D *edepmuon_abs_A = new TH2D("edepmuon_abs_A", "Muons Transverse Momentun [GeV] vs Total Deposited Energy [MeV] for Muons with absorber for layer A",250, 0.1, 2.5, 250, 0, 13);
tree->Draw("total_edep_on_A:Particle_Momentum_GeV>>edepmuon_abs_A", "total_edep_on_A > 0","COLZ");

edepmuon_abs_A->SetTitle("Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons for layer A");
edepmuon_abs_A->GetXaxis()->SetTitle("Particle Momentum [GeV]");
edepmuon_abs_A->GetYaxis()->SetTitle("Total Energy Deposition [MeV]");

canva_bethe_mu_si_A->SetGrid();
edepmuon_abs_A->Draw("COLZ");
canva_bethe_mu_si_A->Update();
canva_bethe_mu_si_A->Update();
canva_bethe_mu_si_A->Draw();




//------------------------------------------------ p_t vs Energy Loss layer B dots ------------------------------------------------

TCanvas* canva_5 = new TCanvas("canva_5",                               
    "Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons in Layer B", 
    1600, 1200);

TGraph* graph_5 = new TGraph();


graph_5->SetMarkerStyle(20);
graph_5->SetMarkerColor(kViolet);
graph_5->SetMarkerSize(0.4);


for (Long64_t i = 0; i < nEntries; i++) {
tree->GetEntry(i);

if ( total_edep_layer_B >= 0) {
graph_5->SetPoint(nPoints, particleMomentum, total_edep_layer_B);
nPoints++;
}
}

graph_5->Draw("AP"); 

// Set titles                                                                         
graph_5->SetTitle("Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with absorber Muons for layer B");
graph_5->GetXaxis()->SetTitle("Particle Momentum [GeV]");
graph_5->GetYaxis()->SetTitle("Total Energy Deposition [MeV]");


canva_5->SetGrid();
canva_5->Update();
canva_5->Draw();



//                                           p_t VS  Energy Loss in layer B (TH2D) 


TCanvas* canva_bethe_mu_si_B = new TCanvas("canva_6", "Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons for layer B", 1600, 1200);
TH2D *edepmuon_abs_B = new TH2D("edepmuon_abs_B", "Muons Transverse Momentun [GeV] vs Total Deposited Energy [MeV] for Muons with absorber for layer B",250, 0.1, 2.5, 250, 0, 13);
tree->Draw("total_edep_on_B:Particle_Momentum_GeV>>edepmuon_abs_B", "total_edep_on_B > 0","COLZ");

edepmuon_abs_B->SetTitle("Transverse Momentum [GeV] vs Total Energy Deposition [MeV] with Absorber for Muons for layer B");
edepmuon_abs_B->GetXaxis()->SetTitle("Particle Momentum [GeV]");
edepmuon_abs_B->GetYaxis()->SetTitle("Total Energy Deposition [MeV]");

canva_bethe_mu_si_B->SetGrid();
edepmuon_abs_B->Draw("COLZ");
canva_bethe_mu_si_B->Update();
canva_bethe_mu_si_B->Update();
canva_bethe_mu_si_B->Draw();



//================================================  TOTAL HITS  ========================================================
TCanvas *Hit1 = new TCanvas("Hit1", "Hits on the module", 1600, 1200);
TH1D *hit1 = new TH1D("hit1", "Hits on module", 2, 0,2);
hit1->SetLineColor(kGreen);
hit1->SetFillColor(kGreen);
hit1->SetFillStyle(3001);



double_t counter_hits = 0;


    for (Long64_t i = 0; i < nEntries; ++i) {
       tree->GetEntry(i);
       hit1->Fill(Hit); 

        counter_hits += Hit; 
       
    }
    cout <<""<<  endl;
    cout<<"TOTAL HITS: "<< counter_hits <<endl;

hit1->GetXaxis()->SetTitle("Hit");
hit1->GetYaxis()->SetTitle("Entries");

hit1->Draw("HIST");

Hit1->SetGrid();
Hit1->Update();


//================================================  EFFICIENCY with HITS  ========================================================
vector<int> momentumCounts(10, 0);                
vector<int> HitCounts(10, 0);

std::vector<double> momentumEdges = { 0.10, 0.34, 0.58, 0.82, 1.06, 1.30, 1.54, 1.78, 2.02, 2.26, 2.50};


for (Long64_t i = 0; i < nEntries; ++i) {

    tree->GetEntry(i);
    int momentumRange = -1;


     for (int j = 0; j < 10; ++j) {
        if (particleMomentum >= momentumEdges[j] && particleMomentum < momentumEdges[j+1]) {
            momentumRange = j;
            break;
            }
        }
        momentumCounts[momentumRange]++;
        //GammaSum[momentumRange] += total_detected_gammas;
        //GammaGeneratedSum[momentumRange] += total_generated_gammas;
    

    if(Hit == 1){
        HitCounts[momentumRange]++;
    }
}
//===========================================
    vector<Double_t> xValues;           //===
//===========================================
    vector<Double_t> yValues;

    cout<< "\n MUONES CON ABSORBER \n" <<endl;
    cout << "Resultados: Hits por rango de momento.\n";

    for (size_t i = 0; i < HitCounts.size(); ++i) {
        Double_t momento = (i * 0.24) + 0.1 ;
        xValues.push_back(momento);

        int HITS = HitCounts[i];
        double_t EFF = ( static_cast<Double_t>(HITS)  /  momentumCounts[i] ) * 100;
        yValues.push_back(EFF);

        cout << "Rango " << i + 1 << " (" << (i * 0.24) + 0.1 << " GeV - " << ((i + 1) * 0.24) + 0.1
             << " GeV):     " << HitCounts[i] << " hits,        " << momentumCounts[i] << "  momentos  en el rango. " 
             << "       Eficiencia:     " << EFF << " %. " << endl;   
    }
cout << "\n"<< endl;


//......................................  EFFICIENCY PLOT  ...................................... 

    TCanvas* canvas_7 = new TCanvas("Efficiency_canvas_7", "Efficiency vs Transverse Momentum for muons with absorber", 1600, 1200);
    TGraph* SiAbs = new TGraph(xValues.size(), &xValues[0], &yValues[0]);
    //TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
    SiAbs->SetTitle("Efficiency vs Transverse Momentum for Muons with absorber");
    SiAbs->SetMarkerStyle(20);
    SiAbs->SetMarkerColor(kRed); 
    SiAbs->SetMarkerSize(1.5);
    SiAbs->GetXaxis()->SetTitle("Momentum [GeV/c]");
  
    SiAbs->GetYaxis()->SetTitle("Efficiency [%]");

    SiAbs->Draw("AP");

    //NoAbs->SetMarkerStyle(21);
    //NoAbs->SetMarkerColor(kBlue);
    
    //NoAbs->Draw("P SAME");

    // Añadir leyenda
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(SiAbs, "WITH ABSORBER", "lp");
    //legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
    legend->Draw();

    canvas_7->SetGrid();
    canvas_7->Update();



//================================================  MEAN ENERGY  ========================================================
 vector<Double_t> edep(10, 0); 
 vector<Double_t> edep_A(10, 0);
 vector<Double_t> edep_B(10, 0);

 for (Long64_t i = 0; i < nEntries; ++i) {

    tree->GetEntry(i);
    int momentumRange = -1;

    for (int j = 0; j < 10; ++j) {
        if (particleMomentum >= momentumEdges[j] && particleMomentum < momentumEdges[j+1]) {
            momentumRange = j;
            break;
            }
        }
        if(Hit == 1){
        edep[momentumRange]+= total_energy_deposition;
        }  
        edep_A[momentumRange] += total_edep_layer_A;
        edep_B[momentumRange] += total_edep_layer_B;
    }

vector<Double_t> yValues_mean_energy;
for (int i = 0; i < 10; ++i) {
        double avgEnergy = (HitCounts[i] > 0) ? edep[i] / HitCounts[i] : 0.0;
        yValues_mean_energy.push_back(avgEnergy);
         
         cout << "Rango " << i + 1 << " (" << (i * 0.24) + 0.1 << " GeV - " << ((i + 1) * 0.24) + 0.1
             << " GeV):     " << HitCounts[i] << " hits,        " << edep[i] << "  suma de energia depositada. " 
             << "       Promedio:     " << avgEnergy << " MeV " << endl; 
    }
cout << "\n"<< endl; 

//------------------------------------------------------------------------------------------------------------------
vector<Double_t> yValues_mean_energy_A;
for (int i = 0; i < 10; ++i) {
        double avgEnergy_A = (HitCounts[i] > 0) ? edep_A[i] / HitCounts[i] : 0.0;
        yValues_mean_energy_A.push_back(avgEnergy_A);
         
         cout << "Rango " << i + 1 << " (" << (i * 0.24) + 0.1 << " GeV - " << ((i + 1) * 0.24) + 0.1
             << " GeV):     " << HitCounts[i] << " hits,        " << edep_A[i] << "  suma de energia depositada en A. " 
             << "       Promedio:     " << avgEnergy_A << " MeV " << endl; 
    }
cout << "\n"<< endl; 

//------------------------------------------------------------------------------------------------------------------
vector<Double_t> yValues_mean_energy_B;
for (int i = 0; i < 10; ++i) {
        double avgEnergy_B = (HitCounts[i] > 0) ? edep_B[i] / HitCounts[i] : 0.0;
        yValues_mean_energy_B.push_back(avgEnergy_B);
        
         
         cout << "Rango " << i + 1 << " (" << (i * 0.24) + 0.1 << " GeV - " << ((i + 1) * 0.24) + 0.1
             << " GeV):     " << HitCounts[i] << " hits,        " << edep_B[i] << "  suma de energia depositada en B. " 
             << "       Promedio:     " << avgEnergy_B << " MeV " << endl; 
    }
cout << "\n"<< endl; 



//...................................|  MEAN ENERGY PLOTS (TOTAL, for A and for B) |...................................

TCanvas* canvas_8 = new TCanvas("average_edep_canvas_8", "Average Energy Loss vs Transverse Momentum for muons with absorber", 1600, 1200);
    TGraph* graph_8 = new TGraph(xValues.size(), &xValues[0], &yValues_mean_energy[0]);
    //TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
    graph_8->SetTitle("Average Energy Loss vs Transverse Momentum for Muons with absorber");
    graph_8->SetMarkerStyle(22);
    graph_8->SetMarkerColor(kBlue-2); 
    graph_8->SetMarkerSize(1.5);
    graph_8->GetXaxis()->SetTitle("Momentum [GeV]");
  
    graph_8->GetYaxis()->SetTitle("Average Energy Loss [MeV]");

    graph_8->Draw("AP");

    //NoAbs->SetMarkerStyle(21);
    //NoAbs->SetMarkerColor(kBlue);
    
    //NoAbs->Draw("P SAME");

    // Añadir leyenda
    TLegend* legend_8 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_8->AddEntry(graph_8, "WITH ABSORBER", "lp");
    //legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
    legend_8->Draw();

    canvas_8->SetGrid();
    canvas_8->Update();



TCanvas* canvas_9 = new TCanvas("average_edep_canvas_9", "Average Energy Loss vs Transverse Momentum for muons with absorber in Layer A", 1600, 1200);
    TGraph* graph_9 = new TGraph(xValues.size(), &xValues[0], &yValues_mean_energy_A[0]);
    //TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
    graph_9->SetTitle("Average Energy Loss vs Transverse Momentum for Muons with absorber  in Layer A");
    graph_9->SetMarkerStyle(22);
    graph_9->SetMarkerColor(kRed+2); 
    graph_9->SetMarkerSize(1.5);
    graph_9->GetXaxis()->SetTitle("Momentum [GeV]");
  
    graph_9->GetYaxis()->SetTitle("Average Energy Loss in Layer A [MeV]");

    graph_9->Draw("AP");

    //NoAbs->SetMarkerStyle(21);
    //NoAbs->SetMarkerColor(kBlue);
    
    //NoAbs->Draw("P SAME");

    // Añadir leyenda
    TLegend* legend_9 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_9->AddEntry(graph_9, "WITH ABSORBER", "lp");
    //legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
    legend_9->Draw();

    canvas_9->SetGrid();
    canvas_9->Update();



TCanvas* canvas_10 = new TCanvas("average_edep_canvas_10", "Average Energy Loss vs Transverse Momentum for muons with absorber  in Layer B", 1600, 1200);
    TGraph* graph_10 = new TGraph(xValues.size(), &xValues[0], &yValues_mean_energy_B[0]);
    //TGraph* NoAbs = new TGraph(xValues2.size(), &xValues2[0], &yValues2[0]);
    graph_10->SetTitle("Average Energy Loss vs Transverse Momentum for Muons with absorber in Layer B");
    graph_10->SetMarkerStyle(22);
    graph_10->SetMarkerColor(kViolet-6); 
    graph_10->SetMarkerSize(1.5);
    graph_10->GetXaxis()->SetTitle("Momentum [GeV]");
  
    graph_10->GetYaxis()->SetTitle("Average Energy Loss in Layer B [MeV]");

    graph_10->Draw("AP");

    //NoAbs->SetMarkerStyle(21);
    //NoAbs->SetMarkerColor(kBlue);
    
    //NoAbs->Draw("P SAME");

    // Añadir leyenda
    TLegend* legend_10 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend_10->AddEntry(graph_10, "WITH ABSORBER", "lp");
    //legend->AddEntry(NoAbs, "NO ABSORBER", "lp");
    legend_10->Draw();

    canvas_10->SetGrid();
    canvas_10->Update();


//===================================================== TOTAL DETECTED PHOTONS  ======================================================

TCanvas *canvas_11 = new TCanvas("canvas_11", "Histogram Total Detected Photons for Muons with absorber", 1600, 1200);
TH1D *photons_detected = new TH1D("Total_detected_photons", "Total detected photons for muons with absorber", 97, 0,5000);

photons_detected->SetLineColor(kGreen);
photons_detected->SetFillColor(kGreen);
photons_detected->SetFillStyle(3001);



for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (total_detected_gammas > 0 ){
    photons_detected->Fill(total_detected_gammas);  
    }   
}
photons_detected->SetTitle("Total detected photons with absorber 70 cm thick  Muons");
photons_detected->GetXaxis()->SetTitle("Photons");
photons_detected->GetYaxis()->SetTitle("Entries");

photons_detected->Draw("HIST");

canvas_11->SetGrid();
canvas_11->Update();


//================================ TOTAL GENERATED PHOTONS ================================

TCanvas *canvas_12 = new TCanvas("canvas_12", "Histogram Total Generated Photons Absorber 70 cm muons", 1600, 1200);
TH1D *photons_generated = new TH1D("Total_generated_photons", "Total generated photons muons", 92, 0,220000);

photons_generated->SetLineColor(kBlue);
photons_generated->SetFillColor(kBlue);
photons_generated->SetFillStyle(3001);


for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (total_generated_gammas > 0 ){
    photons_generated->Fill(total_generated_gammas);  
    }   
}

photons_generated->SetTitle("Total generated photons with absorber 70 cm thick  Muons");
photons_generated->GetXaxis()->SetTitle("Photons");
photons_generated->GetYaxis()->SetTitle("Entries");

photons_generated->Draw("HIST");

canvas_12->SetGrid();
canvas_12->Update();



//=====================================================  [A] TOTAL DETECTED PHOTONS  ======================================================

TCanvas *canvas_13 = new TCanvas("canvas_13", "Histogram Total Detected Photons for Muons with absorber in layer A", 1600, 1200);
TH1D *photons_detected_A = new TH1D("Total_detected_photons", "Total detected photons for muons with absorber in layer A", 92, 0,3500);

photons_detected_A->SetLineColor(kCyan+3);
photons_detected_A->SetFillColor(kCyan+3);
photons_detected_A->SetFillStyle(3004);



for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (total_detected_photons_layer_A > 0 ){
    photons_detected_A->Fill(total_detected_photons_layer_A);  
    }   
}
photons_detected_A->SetTitle("Total detected photons with absorber 70 cm thick Muons in layer A");
photons_detected_A->GetXaxis()->SetTitle("Photons");
photons_detected_A->GetYaxis()->SetTitle("Entries");

photons_detected_A->Draw("HIST");

canvas_13->SetGrid();
canvas_13->Update();


//================================================  [A] TOTAL GENERATED PHOTONS  ========================================================

TCanvas *canvas_14 = new TCanvas("canvas_14", "Histogram Total Generated Photons Absorber 70 cm muons  in layer A", 1600, 1200);
TH1D *photons_generated_A = new TH1D("Total_generated_photons", "Total generated photons muons  in layer A", 92, 0,140000);

photons_generated_A->SetLineColor(kMagenta+2);
photons_generated_A->SetFillColor(kMagenta+2);
photons_generated_A->SetFillStyle(3005);


for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (total_generated_photons_layer_A > 0 ){
    photons_generated_A->Fill(total_generated_photons_layer_A);  
    }   
}

photons_generated_A->SetTitle("Total generated photons with absorber 70 cm thick Muons in layer A");
photons_generated_A->GetXaxis()->SetTitle("Photons");
photons_generated_A->GetYaxis()->SetTitle("Entries");

photons_generated_A->Draw("HIST");

canvas_14->SetGrid();
canvas_14->Update();


//=====================================================  [B] TOTAL DETECTED PHOTONS  ======================================================

TCanvas *canvas_15 = new TCanvas("canvas_15", "Histogram Total Detected Photons for Muons with absorber in layer B", 1600, 1200);
TH1D *photons_detected_B = new TH1D("Total_detected_photons", "Total detected photons for muons with absorber in layer B", 92, 0,4600);

photons_detected_B->SetLineColor(kOrange+7);
photons_detected_B->SetFillColor(kOrange+7);
photons_detected_B->SetFillStyle(3009);



for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (total_detected_photons_layer_B > 0 ){
    photons_detected_B->Fill(total_detected_photons_layer_B);  
    }   
}
photons_detected_B->SetTitle("Total detected photons with absorber 70 cm thick Muons in layer B");
photons_detected_B->GetXaxis()->SetTitle("Photons");
photons_detected_B->GetYaxis()->SetTitle("Entries");

photons_detected_B->Draw("HIST");

canvas_15->SetGrid();
canvas_15->Update();


//================================  [B] TOTAL GENERATED PHOTONS  ================================

TCanvas *canvas_16 = new TCanvas("canvas_16", "Histogram Generated Photons Absorber 70 cm muons in layer B", 1600, 1200);
TH1D *photons_generated_B = new TH1D("Total_generated_photons", "Generated photons in layer B for muons", 91, 0,160000);

photons_generated_B->SetLineColor(kAzure-5);
photons_generated_B->SetFillColor(kAzure-5);
photons_generated_B->SetFillStyle(3010);


for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (total_generated_photons_layer_B > 0 ){
    photons_generated_B->Fill(total_generated_photons_layer_B);  
    }   
}

photons_generated_B->SetTitle("Total generated photons with absorber 70 cm thick Muons in layer B");
photons_generated_B->GetXaxis()->SetTitle("Photons");
photons_generated_B->GetYaxis()->SetTitle("Entries");

photons_generated_B->Draw("HIST");

canvas_16->SetGrid();
canvas_16->Update();



//======================================  TOTAL / AVERAGE PHOTONS DETECTION / GENERATION VS MOMENTUM  =========================================





//======================================  DISPERSSION ANGLE RESPECT TO AXIS Y  ======================================


// In layer B

TCanvas *canvas_36 = new TCanvas("canvas_36", "Histogram angle dispersion Abs 70 cm for layer B", 1600, 1200);
TH1D *graph_36 = new TH1D("dispersion_angle_B", "Angle dispersion in layer B", 93, 0, 6.2);

graph_36->SetLineColor(kSpring+2);
graph_36->SetFillColor(kSpring+2);
graph_36->SetFillStyle(3011);



for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (angle_b >= 0 ){
        graph_36->Fill(angle_b);  
   }   
}

graph_36->SetTitle("Dispersion angle by hits in with absorber 70 cm thick Muons in layer B");
graph_36->GetXaxis()->SetTitle("Degrees °");
graph_36->GetYaxis()->SetTitle("Entries");

graph_36->Draw("HIST");

canvas_36->SetGrid();
canvas_36->Update();


//IN layer A

TCanvas *canvas_37 = new TCanvas("canvas_37", "Histogram angle dispersion Abs 70 cm for layer A", 1600, 1200);
TH1D *graph_37 = new TH1D("dispersion_angle_B", "Angle dispersion in layer A", 92, 0, 6.2);

graph_37->SetLineColor(kAzure+7);
graph_37->SetFillColor(kAzure+7);
graph_37->SetFillStyle(3012);



for (Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (angle_a >= 0 ){
        graph_37->Fill(angle_a);  
   }   
}

graph_37->SetTitle("Dispersion angle by hits in with absorber 70 cm thick Muons in layer A");
graph_37->GetXaxis()->SetTitle("Degrees °");
graph_37->GetYaxis()->SetTitle("Entries");

graph_37->Draw("HIST");

canvas_37->SetGrid();
canvas_37->Update();

}





