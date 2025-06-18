//    DETECTOR CONSTRUCTION HEADER

#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalPhysics.hh"
#include "SensitiveDetector.hh"
#include "G4GenericMessenger.hh" 
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh" 
#include "G4LogicalBorderSurface.hh"
#include "G4SubtractionSolid.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4MagneticField.hh"

#include "G4Cache.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;

class G4UniformMagField;
class F01FieldSetup;
class G4VisAttributes;




class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    ~DetectorConstruction();
  
    virtual G4VPhysicalVolume *Construct(); 

    void ConstructMIDModule();     


    // TO RETRIEVE THE fScoringVolumes FOR STEPPING ACTION
    const std::vector<G4LogicalVolume*>& GetScoringVolumesA() const {
    return fScoringVolumes_A;
    }

    const std::vector<G4LogicalVolume*>& GetScoringVolumesB() const {
    return fScoringVolumes_B;
    }


   std::vector<G4LogicalVolume*> LogicBars_A;
   std::vector<G4LogicalVolume*> LogicSiPMs_A;
   std::vector<G4LogicalVolume*> fScoringVolumes_A;
      std::vector<G4LogicalVolume*> Logic_Fibers_A;
      std::vector<G4LogicalVolume*> Logic_claddings_A;

   std::vector<G4LogicalVolume*> LogicBars_B; 
   std::vector<G4LogicalVolume*> LogicSiPMs_B;
   std::vector<G4LogicalVolume*> fScoringVolumes_B;
      std::vector<G4LogicalVolume*> Logic_Fibers_B;
      std::vector<G4LogicalVolume*> Logic_claddings_B;


    G4double wavelength, lightOutput;
    G4bool MID_Module; 
    G4double distance_modules; 
    
  
  private: 

    virtual void ConstructSDandField();

    G4double env_sizeX, env_sizeY, env_sizeZ; 

    G4GenericMessenger *fMessenger;  


    G4Tubs *Solid_fiber, *Solid_cladding;
    G4Box *window;

    G4Box  *SolidWorld, *Solidbar, *Solidsipm, *SolidCube,
           *SolidSA, *Solidmylar, *SolidInnermylar;

    G4LogicalVolume *LogicWorld, *Logicbar, *Logicmylar, *Logicsipm, *LogicalSA,
                    *LogicCube, *Logicbar_A, *Logicbar_B, *Logicsipm_A, *Logicsipm_B,
                    *Logic_Fiber_A, *Logic_Fiber_B,
                    *Logic_cladding_A, *Logic_cladding_B;;

    G4VPhysicalVolume *PhysicalWorld, *Physicalbar, *Physicalsipm, *Physicalmylar, *PhysicalSA, *PhysicalCube,
                      *Physical_MID_A, *Physical_SiPM_MID_A, *Physical_Mylar_MID_A, 
                      *Physical_MID_B, *Physical_SiPM_MID_B, *Physical_Mylar_MID_B,
                      *Physical_Fiber_A, *Physical_Fiber_B,
                      *Physical_cladding_A, *Physical_cladding_B;

    
    void DefineMaterials();

    G4Material *plastic, *worldMaterial, *steel, *mylarMaterial, *fiber_core, *cladding;
    G4Element *C, *H, *O;
    G4OpticalSurface *mirrorsurface, *op_surface_bar_clad, *op_surface_core_clad, *op_surface_clad_air, *op_surface_bar_air; 


    G4Cache<F01FieldSetup*>    fEmFieldSetup;

        
    //vector<G4VisAttributes*> fVisAttributes;
    
  
};

#endif

