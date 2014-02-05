// -*- C++ -*-
//
// Package:    TrackAndPointsProducer
// Class:      TrackAndPointsProducer
// 
/**\class TrackAndPointsProducer TrackAndPointsProducer.cc prodTutorial/TrackAndPointsProducer/src/TrackAndPointsProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Tue Jan 28 14:03:17 PST 2014
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
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

//
// class declaration
//

class TrackAndPointsProducer : public edm::EDProducer {
   public:
      explicit TrackAndPointsProducer(const edm::ParameterSet&);
      ~TrackAndPointsProducer();

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

  edm::InputTag src_;
  typedef math::XYZPointD Point;
  typedef std::vector<Point> PointCollection;

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
TrackAndPointsProducer::TrackAndPointsProducer(const edm::ParameterSet& iConfig)
{
  //register your products
  
  src_ = iConfig.getParameter<edm::InputTag>( "src" );
  produces<PointCollection>( "innerPoint" ).setBranchAlias( "innerPoints");
  produces<PointCollection>( "outerPoint" ).setBranchAlias( "outerPoints");

  //now do what ever other initialization is needed

}


TrackAndPointsProducer::~TrackAndPointsProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TrackAndPointsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;


   // retrieve the tracks
   Handle<TrackCollection> tracks;
   iEvent.getByLabel( src_, tracks );

   // create the vectors. Use auto_ptr, as these pointers will automatically
   // delete when they go out of scope, a very efficient way to reduce memory leaks.
   auto_ptr<PointCollection> innerPoints( new PointCollection );
   auto_ptr<PointCollection> outerPoints( new PointCollection );

   // and already reserve some space for the new data, to control the size
   // of your executible's memory use.
   const int size = tracks->size();
   innerPoints->reserve( size );
   outerPoints->reserve( size );

   // loop over the tracks:
   for( TrackCollection::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {
	 // fill the points in the vectors
	 innerPoints->push_back( track->innerPosition() );
	 outerPoints->push_back( track->outerPosition() );
   }

   // and save the vectors
   iEvent.put( innerPoints, "innerPoint" );
   iEvent.put( outerPoints, "outerPoint" );

}

// ------------ method called once each job just before starting event loop  ------------
void 
TrackAndPointsProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrackAndPointsProducer::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
TrackAndPointsProducer::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TrackAndPointsProducer::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TrackAndPointsProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TrackAndPointsProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackAndPointsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAndPointsProducer);
