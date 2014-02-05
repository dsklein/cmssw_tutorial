// -*- C++ -*-
//
// Package:    ZEventProducer
// Class:      ZEventProducer
// 
/**\class ZEventProducer ZEventProducer.cc zAnalysis/ZEventProducer/src/ZEventProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu Jan 30 10:27:36 PST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

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

class ZEventProducer : public edm::EDProducer {
   public:
      explicit ZEventProducer(const edm::ParameterSet&);
      ~ZEventProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

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
ZEventProducer::ZEventProducer(const edm::ParameterSet& iConfig)
{
   //register your products

  //Gen-level information
  produces<LorentzVector>( "parentZp4" ).setBranchAlias( "z_p4" );
  produces<std::vector<LorentzVector> >( "daughterp4s" ).setBranchAlias( "lep_p4s" );
  produces<std::vector<int> >( "daughterIds" ).setBranchAlias( "lep_ids" );

  //Reco-level information
  produces<std::vector<LorentzVector> >( "recoLeptonp4s" ).setBranchAlias( "recolep_p4s" );
  produces<std::vector<int> >( "recoLeptonIds" ).setBranchAlias( "recolep_ids" );

   //now do what ever other initialization is needed
  
}


ZEventProducer::~ZEventProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ZEventProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   // Retrieve the info we want
   Handle<vector<GenParticle> > genPs;
   iEvent.getByLabel( "genParticles", genPs );

   Handle<vector<Muon> > recomuons;
   iEvent.getByLabel( "muons", recomuons );

   Handle<vector<GsfElectron> > recoelectrons;
   iEvent.getByLabel( "gsfElectrons", recoelectrons );

   // Make our vectors, using auto_ptrs
   auto_ptr<LorentzVector>           z_p4( new LorentzVector );
   auto_ptr<vector<LorentzVector> >  lep_p4s( new vector<LorentzVector> );
   auto_ptr<vector<LorentzVector> >  recolep_p4s( new vector<LorentzVector> );
   auto_ptr<vector<int> >            lep_ids( new vector<int> );
   auto_ptr<vector<int> >            recolep_ids( new vector<int> );

   // Take the info from our file and store it in our own containers

   //iterate over the gen particles, find which one is a Z, get its p4 and store it
   vector<GenParticle>::const_iterator pit = genPs->begin();
   for( pit = genPs->begin(); pit != genPs->end(); ++pit) {

	 if( pit->pdgId() != 23 ) continue;
	 if( pit->status() != 3 ) continue;
	 *z_p4 = LorentzVector( pit->p4() );
	 
	 //Find all daughters of the Z; get and store their p4s and pdgIds
	 int ndaughters = pit->numberOfDaughters();            //better to use size_type; will change if anything complains
	 for( int i=0; i<ndaughters; i++) {
	   if( pit->daughter(i)->pdgId() == 23 ) continue;
	   lep_p4s->push_back( LorentzVector(pit->daughter(i)->p4() ) );
	   lep_ids->push_back( pit->daughter(i)->pdgId() );
	 }
   }

   //iterate over the reco electrons, get each one's p4 and pdgId, and push_back into our vector
   vector<GsfElectron>::const_iterator elit = recoelectrons->begin();
   for( elit = recoelectrons->begin(); elit != recoelectrons->end(); ++elit) {
	 recolep_p4s->push_back( LorentzVector( elit->p4() ) );
	 recolep_ids->push_back( elit->pdgId() );
   }

   //iterate over the reco muons, get each one's p4 and pdgId, and push_back into our vector
   vector<Muon>::const_iterator muit = recomuons->begin();
   for( muit = recomuons->begin(); muit != recomuons->end(); ++muit) {
	 recolep_p4s->push_back( LorentzVector(muit->p4() ) );
	 recolep_ids->push_back( muit->pdgId() );
   }

   //cout << "Writing this event" << endl;

   // Save our containers
   iEvent.put( z_p4,           "parentZp4" );
   iEvent.put( lep_p4s,         "daughterp4s" );
   iEvent.put( recolep_p4s,     "recoLeptonp4s" );
   iEvent.put( lep_ids,      "daughterIds" );
   iEvent.put( recolep_ids, "recoLeptonIds" );

 
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZEventProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZEventProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
ZEventProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZEventProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZEventProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZEventProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZEventProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZEventProducer);
