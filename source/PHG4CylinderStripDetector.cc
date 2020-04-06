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
  G4Material *TrackerMaterial = G4Material::GetMaterial(m_Params->get_string_param("material"));

  if (!TrackerMaterial)
  {
    std::cout << "Error: Can not set material" << std::endl;
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
                                                        G4Material::GetMaterial("G4_AIR"),
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
