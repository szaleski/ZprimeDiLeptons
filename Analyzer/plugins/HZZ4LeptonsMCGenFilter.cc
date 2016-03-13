/* \class HZZ4leptonsMCGenFilter
 *
 *  Filter of Z->tau tau leptons channel at the level of generation
 *  author:  Nicola De Filippis
 *
 */

// system include files
#include "ZprimeDiLeptons/Analyzer/interface/HZZ4LeptonsMCGenFilter.h"

// User include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// C++
#include <iostream>
#include <vector>

// Constructor
HZZ4LeptonsMCGenFilter::HZZ4LeptonsMCGenFilter(const edm::ParameterSet& pset) {
  
  // LeptonFlavour for Z final states
  // 16 = z-> tautau

  // Local Debug flag
  // gen                = pset.getParameter<edm::InputTag>("genParticles"); 
  gen                = consumes<std::vector<reco::GenParticle> >(pset.getParameter<edm::InputTag>("genParticles"));
  debug              = pset.getParameter<bool>("DebugHZZ4LeptonsMCGenFilter");
  leptonFlavour      = pset.getParameter<int>("HZZ4LeptonsMCFilterLeptonFlavour");
  acceptance         = pset.getParameter<double>("acceptance");
  
  ikept   =0;
  evt     =0;
  
}

// Destructor
HZZ4LeptonsMCGenFilter::~HZZ4LeptonsMCGenFilter() {
  
  std::cout << "number of events processed: " << evt << std::endl;
  std::cout << "number of events kept: " << ikept << std::endl;
  std::cout << "expected number of taus = " << taus << std::endl;
}

// Filter event
bool HZZ4LeptonsMCGenFilter::filter(edm::Event& event, const edm::EventSetup& setup ) {
  
  bool keepEvent   = false;
  evt++;
  
  bool ztautau  = false;
  
  // get gen particle candidates 
  edm::Handle<reco::GenParticleCollection> genCandidates;	    
  // event.getByLabel(gen, genCandidates);
  event.getByToken(gen, genCandidates);
  
  int nTau = 0;
  
  for ( reco::GenParticleCollection::const_iterator mcIter=genCandidates->begin(); mcIter!=genCandidates->end(); ++mcIter ) {

    // Taus:
    if ( mcIter->pdgId() == 15 || mcIter->pdgId() == -15) {
      // Mother is a Z
      if ( mcIter->mother()->pdgId() == 23) {
        nTau++;
      }
    }
    
  }
  
  if (nTau  == 2) ztautau = true;
  if ( leptonFlavour == 16 && ztautau ) keepEvent = true;  
  
  if (keepEvent ) ikept++;
  
  return keepEvent;
  
}


void HZZ4LeptonsMCGenFilter::endJob(){
  
  std::cout << "Number of events filtered in the acceptance= " << ikept << std::endl;
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(HZZ4LeptonsMCGenFilter);

