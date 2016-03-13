#ifndef HZZ4LeptonsMCGenFilter_h
#define HZZ4LeptonsMCGenFilter_h

/* \class HZZ4LeptonsMCGenFilter
 *
 *
 * Filter to select 4 lepton events (4e, 4mu, 4leptons) 
 *
 * \author Nicola De Filippis
 *
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


class HZZ4LeptonsMCGenFilter : public edm::EDFilter {
  
 public:
  // Constructor
  explicit HZZ4LeptonsMCGenFilter(const edm::ParameterSet&);

  // Destructor
  ~HZZ4LeptonsMCGenFilter();

  /// Get event properties to send to builder to fill seed collection
  virtual bool filter(edm::Event&, const edm::EventSetup& );

  void endJob();

 private:
  int evt, ikept, taus;
  // edm::InputTag gen;  
  edm::EDGetTokenT<std::vector<reco::GenParticle> > gen;     
  bool debug;
  int leptonFlavour;
  double acceptance;
};

#endif
