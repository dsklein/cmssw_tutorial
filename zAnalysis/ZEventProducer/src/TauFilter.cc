// -*- C++ -*-
//
// Package:    TauFilter
// Class:      TauFilter
// 
/**\class TauFilter TauFilter.cc zAnalysis/TauFilter/src/TauFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu Jan 30 10:29:08 PST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <iostream>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

//
// class declaration
//

class TauFilter : public edm::EDFilter {
   public:
      explicit TauFilter(const edm::ParameterSet&);
      ~TauFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

  typedef math::XYZTLorentzVector LorentzVector;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TauFilter::TauFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


TauFilter::~TauFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TauFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   Handle<vector<GenParticle> > genPs;
   iEvent.getByLabel( "genParticles", genPs );

   //iterate over the gen particles, find which one is a Z
   vector<GenParticle>::const_iterator pit = genPs->begin();
   for( pit = genPs->begin(); pit != genPs->end(); ++pit) {

	 if( pit->pdgId() != 23 ) continue;
	 
	 int ndaughters = pit->numberOfDaughters();            //better to use size_type; will change if anything complains
	 int daughterId = 0;

	 //Find the daughters of the Z; if you meet a tau, return false
	 for( int i=0; i<ndaughters; i++) {
	   daughterId = abs( pit->daughter(i)->pdgId() );
	   if( daughterId == 23 ) continue;
	   if( daughterId == 15 ) return false;
	 }
   }

   //cout << "Keeping this event" << endl;
   return true;

}

// ------------ method called once each job just before starting event loop  ------------
void 
TauFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
TauFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
TauFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
TauFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
TauFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TauFilter);
