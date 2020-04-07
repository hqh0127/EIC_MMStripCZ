#include "PHG4CylinderStripDetector.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4Detector.h>       // for PHG4Detector
#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4String.hh>  // for G4String
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4ThreeVector.hh>  // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Element.hh>
#include <Geant4/G4Material.hh>

#include <TSystem.h>

#include <cmath>
#include <iostream>  // for operator<<, endl, basic_ost...
#include <sstream>

class G4VSolid;
class PHCompositeNode;

using namespace std;

//_______________________________________________________________
PHG4CylinderStripDetector::PHG4CylinderStripDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam, const int lyr)
  : PHG4Detector(subsys, Node, dnam)
  , m_Params(parameters)
  , m_CylinderPhysicalVolume(nullptr)
  , m_Layer(lyr)
{
}

//_______________________________________________________________
bool PHG4CylinderStripDetector::IsInDetector(G4VPhysicalVolume *volume) const
{
  set<G4VPhysicalVolume *>::const_iterator iter =
      m_PhysicalVolumesSet.find(volume);
  if (iter != m_PhysicalVolumesSet.end())
  {
    return true;
  }

  return false;
}

//_______________________________________________________________

bool PHG4CylinderStripDetector::IsInTileC(G4VPhysicalVolume *volume) const
{
  if (volume == m_CylinderCPhysicalVolume)
  {
    return true;
  }

  return false;
}

//_______________________________________________________________

bool PHG4CylinderStripDetector::IsInTileZ(G4VPhysicalVolume *volume) const
{
  if (volume == m_CylinderZPhysicalVolume)
  {
    return true;
  }

  return false;
}

//_______________________________________________________________
void PHG4CylinderStripDetector::ConstructMe(G4LogicalVolume *logicWorld)
{
  BuildMaterials();
  G4Material *TrackerMaterial = G4Material::GetMaterial(m_Params->get_string_param("gas"));

  if (!TrackerMaterial)
  {
    std::cout << "Error: Can not set Micromegas Gas" << std::endl;
    gSystem->Exit(1);
  }

  // determine length of cylinder using PHENIX's rapidity coverage if flag is true
  double radius = m_Params->get_double_param("radius") * cm;
  double thickness = m_Params->get_double_param("thickness") * cm;
  double gap = m_Params->get_double_param("gap") * cm;
  G4VSolid *cylinder_solid = new G4Tubs(G4String(GetName()),
                                        radius,
                                        radius + thickness*2 + gap,
                                        m_Params->get_double_param("length") * cm / 2., 0, twopi);
  G4VSolid *cylinder_solid_C = new G4Tubs(G4String(GetName()),
                                        radius,
                                        radius + thickness,
                                        m_Params->get_double_param("length") * cm / 2., 0, twopi);
  G4VSolid *cylinder_solid_Z = new G4Tubs(G4String(GetName()),
                                        radius + thickness + gap,
                                        radius + thickness*2 + gap,
                                        m_Params->get_double_param("length") * cm / 2., 0, twopi);
  double steplimits = m_Params->get_double_param("steplimits") * cm;
  G4UserLimits *g4userlimits = nullptr;
  if (isfinite(steplimits))
  {
    g4userlimits = new G4UserLimits(steplimits);
  }

  G4LogicalVolume *cylinder_logic = new G4LogicalVolume(cylinder_solid,
                                                        G4Material::GetMaterial("myAir"),
                                                        G4String(GetName())
                                                        );
  G4LogicalVolume *cylinder_logic_C = new G4LogicalVolume(cylinder_solid_C,
                                                        TrackerMaterial,
                                                        G4String(GetName())+"CTileLogic",
                                                        nullptr, nullptr, g4userlimits);
  G4LogicalVolume *cylinder_logic_Z = new G4LogicalVolume(cylinder_solid_Z,
                                                        TrackerMaterial,
                                                        G4String(GetName())+"ZTileLogic",
                                                        nullptr, nullptr, g4userlimits);
  G4VisAttributes *vis = new G4VisAttributes(G4Color(G4Colour::Grey())); // grey is good to see the tracks in the display
  vis->SetForceSolid(true);
  cylinder_logic_C->SetVisAttributes(vis);
  cylinder_logic_Z->SetVisAttributes(vis);
  
  vis = new G4VisAttributes(G4Color(G4Colour::Grey())); // grey is good to see the tracks in the display
  vis->SetForceSolid(true);
  vis->SetVisibility(false);
  cylinder_logic->SetVisAttributes(vis);

  PHG4Subsystem *mysys = GetMySubsystem();
  mysys->SetLogicalVolume(cylinder_logic);
  m_CylinderCPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0,0,0),
                                               cylinder_logic_C,
                                               G4String(GetName())+"CTilePhys",
                                               cylinder_logic, 0, false, OverlapCheck());
  m_PhysicalVolumesSet.insert(m_CylinderCPhysicalVolume);
  m_CylinderZPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0,0,0),
                                               cylinder_logic_Z,
                                               G4String(GetName())+"ZTilePhys",
                                               cylinder_logic, 0, false, OverlapCheck());
  m_PhysicalVolumesSet.insert(m_CylinderZPhysicalVolume);
  m_CylinderPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(m_Params->get_double_param("place_x") * cm, m_Params->get_double_param("place_y") * cm, m_Params->get_double_param("place_z") * cm),
                                               cylinder_logic,
                                               G4String(GetName()),
                                               logicWorld, 0, false, OverlapCheck());
}

void PHG4CylinderStripDetector::BuildMaterials(){
  // get the list of chemical elements 
  // ---------------------------------
  
  //G4Element *N  = G4Element::GetElement( "N" );
  //G4Element *O  = G4Element::GetElement( "O" );
  //G4Element *H  = G4Element::GetElement( "H" );
  //G4Element *C  = G4Element::GetElement( "C" );
  
  // get the list of NIST materials 
  // ---------------------------------
  G4Material *G4_N = G4Material::GetMaterial("G4_N");
  G4Material *G4_O = G4Material::GetMaterial("G4_O");
  G4Material *G4_C = G4Material::GetMaterial("G4_C");
  G4Material *G4_H = G4Material::GetMaterial("G4_H");
  G4Material *G4_Si = G4Material::GetMaterial("G4_Si");
  G4Material *G4_Ar = G4Material::GetMaterial("G4_Ar");
  G4Material *G4_Cr = G4Material::GetMaterial("G4_Cr");
  G4Material *G4_Fe = G4Material::GetMaterial("G4_Fe");
  G4Material *G4_Mn = G4Material::GetMaterial("G4_Mn");
  G4Material *G4_Ni = G4Material::GetMaterial("G4_Ni");
  G4Material *G4_Cu = G4Material::GetMaterial("G4_Cu");

  // combine elements
  // ----------------
  G4int ncomponents;
  G4double fraction;
  G4double temperature = 298.15*kelvin;
  G4double pressure = 1.*atmosphere;
  

  // air
  if (!G4Material::GetMaterial("myAir")){
    G4Material *myAir = new G4Material( "myAir", 0.001205, ncomponents=2, kStateGas, temperature, pressure);
    myAir->AddMaterial( G4_N, fraction = 0.77 );
    myAir->AddMaterial( G4_O, fraction = 0.23 );
  }

  // FR4
  if (!G4Material::GetMaterial("myFR4")){
    G4Material *myFR4 = new G4Material( "myFR4", 1.860, ncomponents=4);
    myFR4->AddMaterial( G4_C,  fraction = 0.43550 );
    myFR4->AddMaterial( G4_H,  fraction = 0.03650 );
    myFR4->AddMaterial( G4_O,  fraction = 0.28120 );
    myFR4->AddMaterial( G4_Si, fraction = 0.24680 );
  }

  // Kapton
  if (!G4Material::GetMaterial("myKapton")){
    G4Material *myKapton = new G4Material( "myKapton", 1.420, ncomponents=4);
    myKapton->AddMaterial( G4_C, 0.6911330 );
    myKapton->AddMaterial( G4_H, 0.0263620 );
    myKapton->AddMaterial( G4_N, 0.0732700 );
    myKapton->AddMaterial( G4_O, 0.2092350);
  }

  // MMgas
  if (!G4Material::GetMaterial("myMMGas")){
    G4Material *myMMGas = new G4Material( "myMMGas", 0.00170335, ncomponents=3);
    myMMGas->AddMaterial( G4_Ar, 0.900 );
    myMMGas->AddMaterial( G4_C,  0.0826586 );
    myMMGas->AddMaterial( G4_H,  0.0173414 );
  }

  // MMMesh
  if (!G4Material::GetMaterial("myMMMesh")){
    G4Material *myMMMesh = new G4Material( "myMMMesh", 2.8548, ncomponents=5);
    myMMMesh->AddMaterial( G4_Cr, 0.1900 );
    myMMMesh->AddMaterial( G4_Fe, 0.6800 );
    myMMMesh->AddMaterial( G4_Mn, 0.0200 );
    myMMMesh->AddMaterial( G4_Ni, 0.1000 );
    myMMMesh->AddMaterial( G4_Si, 0.0100 );
  }

  // MMStrips
  if (!G4Material::GetMaterial("myMMStrips")){
    G4Material *myMMStrips = new G4Material( "myMMStrips", 5.248414, G4_Cu);
    cout << myMMStrips->GetName() << endl;
  }

  // MMResistivePaste
  if (!G4Material::GetMaterial("myMMResistivePaste")){
    G4Material *myMMResistivePaste = new G4Material( "myMMResistivePaste", 0.77906, G4_C);
    cout << myMMResistivePaste->GetName() << endl;
  }

  // Copper
  if (!G4Material::GetMaterial("myCopper")){
    G4Material *myCopper = new G4Material("myCopper", 8.9600, G4_Cu);
    cout << myCopper->GetName() << endl;
  }
}
