//==========================================================================
//  Code to extraxt fluences directly from G4 steps
//--------------------------------------------------------------------------
// 
// Author: Alex Jentsch
//
//
// Based on AIDASoft "TestSteppingAction" code 
//
//
//==========================================================================

// Framework include files
#include "DDG4/Geant4SteppingAction.h"
#include "DDG4/Geant4StepHandler.h"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim {
    
    /// Class to count steps and suspend tracks every 5 steps
    /** Class to count steps and suspens tracks every 5 steps
     * Count steps and suspend
     *
     *  \version 1.0
     *  \ingroup DD4HEP_SIMULATION
     */
    class FluenceFromG4Step : public Geant4SteppingAction {
      std::size_t m_calls_steps { 0UL };
      std::size_t m_calls_suspended { 0UL };
      std::size_t m_calls_kill { 0UL };
    
    public:
      /// Standard constructor
      FluenceFromG4Step(Geant4Context* context, const std::string& nam) : Geant4SteppingAction(context, nam) 
      {
      }
      /// Default destructor
      virtual ~FluenceFromG4Step()   {
		info("+++ Track Calls Steps: %ld", m_calls_steps);
		info("+++ Track Calls Suspended: %ld", m_calls_suspended);
		info("+++ Track Calls Killed: %ld", m_calls_kill);
      }
      /// stepping callback
      /*virtual void operator()(const G4Step* step, G4SteppingManager*) {
        if(m_calls_steps % 5 == 0 ) {
          ++m_calls_suspended;
          step->GetTrack()->SetTrackStatus(fSuspend);
        } else if((m_calls_steps + 1) % 30 == 0 ) {
          ++m_calls_kill;
          step->GetTrack()->SetTrackStatus(fStopAndKill);
        }
		++m_calls_steps;
      }*/
    };
  }    // End namespace sim
}      // End namespace dd4hep

#include "DDG4/Factories.h"
DECLARE_GEANT4ACTION_NS(dd4hep::sim, FluenceFromG4Step)
