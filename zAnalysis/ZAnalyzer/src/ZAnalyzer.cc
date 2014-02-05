// -*- C++ -*-
//
// Package:    ZAnalyzer
// Class:      ZAnalyzer
// 
/**\class ZAnalyzer ZAnalyzer.cc zAnalysis/ZAnalyzer/src/ZAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu Jan 30 16:38:05 PST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/Math/interface/deltaR.h>
#include "TH1.h"
#include "TH2.h"


//
// class declaration
//

class ZAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ZAnalyzer(const edm::ParameterSet&);
      ~ZAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------

  typedef math::XYZTLorentzVector LorentzVector;

  TH1D* h_zMass;

  TH1D* h_dilMass_gen;
  TH1D* h_dilMass_reco;
  TH1D* h_eeMass_gen;
  TH1D* h_eeMass_reco;
  TH1D* h_mumuMass_gen;
  TH1D* h_mumuMass_reco;

  TH1D* h_lepPt_gen;
  TH1D* h_lepPt_reco;
  TH1D* h_lepEta_gen;
  TH1D* h_lepEta_reco;

  TH1D* h_effEE;
  TH1D* h_effMu;
  TH1D* h_effEE_numer;
  TH1D* h_effMu_numer;
  TH1D* h_effEE_denom;
  TH1D* h_effMu_denom;

  // TH1D* h_resolution_zmass;
  // TH1D* h_resolution_elpt;
  // TH1D* h_resolution_mupt;

  TH2D* h_residual_zmass;
  TH2D* h_residual_elpt;
  TH2D* h_residual_mupt;

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
ZAnalyzer::ZAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  TH1::SetDefaultSumw2();

  edm::Service<TFileService> fs;
  h_zMass =          fs->make<TH1D>("zmass",        "Mass of the Z", 30, 0, 150);
  h_dilMass_gen =    fs->make<TH1D>("dilmassgen",   "Dilepton mass (gen)", 30, 0, 150);
  h_dilMass_reco =   fs->make<TH1D>("dilmassreco",  "Dilepton mass (reco)", 30, 0, 150);
  h_eeMass_gen =     fs->make<TH1D>("eemassgen",    "M_{ee} (gen)", 30, 0, 150);
  h_eeMass_reco =    fs->make<TH1D>("eemassreco",   "M_{ee} (reco)", 30, 0, 150);
  h_mumuMass_gen =   fs->make<TH1D>("mumumassgen",  "M_{#mu#mu} (gen)", 30, 0, 150);
  h_mumuMass_reco =  fs->make<TH1D>("mumumassreco", "M_{#mu#mu} (reco)", 30, 0, 150);  

  h_lepPt_gen =      fs->make<TH1D>("lepptgen",     "Lepton p_{T} (gen)", 40, 0, 200);
  h_lepPt_reco =     fs->make<TH1D>("lepptreco",    "Lepton p_{T} (reco)", 40, 0, 200);
  h_lepEta_gen =     fs->make<TH1D>("lepetagen",    "Lepton #eta (gen)", 30, -3, 3);
  h_lepEta_reco =    fs->make<TH1D>("lepetareco",   "Lepton #eta (reco)", 30, -3, 3);

  h_effEE =          fs->make<TH1D>("electroneff",     "Electron efficiency", 40, 0, 200);
  h_effMu =          fs->make<TH1D>("muoneff",         "Muon efficiency", 40, 0, 200);
  h_effEE_numer =    fs->make<TH1D>("electroneffnum",  "Electron efficiency numerator", 40, 0, 200);
  h_effMu_numer =    fs->make<TH1D>("muoneffnum",      "Muon efficiency numerator", 40, 0, 200);
  h_effEE_denom =    fs->make<TH1D>("electroneffden",  "Electron efficiency denominator", 40, 0, 200);
  h_effMu_denom =    fs->make<TH1D>("muoneffden",      "Muon efficiency denominator", 40, 0, 200);

  // h_resolution_zmass = fs->make<TH1D>("resozmass",  "Z mass resolution", 30, 0, 150);
  // h_resolution_elpt  = fs->make<TH1D>("resoelpt",   "Electron pt resolution", 40, 0, 200);
  // h_resolution_mupt  = fs->make<TH1D>("resomupt",   "Muon pt resolution", 40, 0, 200);

  h_residual_zmass = fs->make<TH2D>("residzmass",   "Z mass residual plot", 30, 0, 150, 20, -100, 100);
  h_residual_elpt  = fs->make<TH2D>("residelpt",    "Electron pt residual plot", 40, 0, 200, 20, -100, 100);
  h_residual_mupt  = fs->make<TH2D>("residmupt",    "Muon pt residual plot", 40, 0, 200, 20, -100, 100);
}


ZAnalyzer::~ZAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   //Get our event information from the stored data file
   Handle<LorentzVector> z_p4;
   iEvent.getByLabel("zProducer", "parentZp4", z_p4);

   Handle<vector<LorentzVector> > lep_p4s;
   iEvent.getByLabel("zProducer", "daughterp4s", lep_p4s);

   Handle<vector<LorentzVector> > recolep_p4s;
   iEvent.getByLabel("zProducer", "recoLeptonp4s", recolep_p4s);

   Handle<vector<int> > lep_ids;
   iEvent.getByLabel("zProducer", "daughterIds", lep_ids);

   Handle<vector<int> > recolep_ids;
   iEvent.getByLabel("zProducer", "recoLeptonIds", recolep_ids);


   // Okay, let's get to work!
   int lep_flavor = 99;
   lep_flavor = abs( lep_ids->at(0) );

   //First, fill the histograms that don't require any gen-reco matching
   h_zMass->Fill( z_p4->mass() );

   vector<LorentzVector>::const_iterator lepit = lep_p4s->begin();
   for( lepit = lep_p4s->begin(); lepit != lep_p4s->end(); ++lepit ) {
	 h_lepPt_gen->Fill( lepit->Pt() );
	 h_lepEta_gen->Fill( lepit->Eta() );
   }

   for( lepit = recolep_p4s->begin(); lepit != recolep_p4s->end(); ++lepit ) {
	 h_lepPt_reco->Fill( lepit->Pt() );
	 h_lepEta_reco->Fill( lepit->Eta() );
   }


   ///////////////////////////// Now, let's do some matching of reco leptons to gen leptons
   vector<LorentzVector> recoleps_matched;
   LorentzVector bestmatch;

   //Loop through all gen leptons
   for( unsigned int i=0; i<lep_p4s->size(); i++ ) {
	 double dr_min = 0.4;
	 double dr_test = 99.9;
	 bool matchfound = false;

	 //Loop through all reco leptons
	 for( unsigned int j=0; j<recolep_p4s->size(); j++ ) {

	   if(recolep_ids->at(j) != lep_ids->at(i) ) continue;   //ignore reco leptons of a different flavor

	   //Check if gen and reco are well-matched in dR
	   dr_test = reco::deltaR( lep_p4s->at(i), recolep_p4s->at(j) );
	   if( dr_test < dr_min ) {
		 dr_min = dr_test;
		 bestmatch = recolep_p4s->at(j);
		 matchfound = true;
	   }
	 }

	 //If there's a good match, store its p4
	 if( matchfound==true ) recoleps_matched.push_back( bestmatch );

	 //Fill efficiency numerator/denominator, for leptons with eta<=2.5
	 if( fabs( lep_p4s->at(i).Eta() ) <= 2.5 )
	   {
		 if( lep_flavor==11 )       h_effEE_denom->Fill(lep_p4s->at(i).Pt() );
		 else if( lep_flavor==13 )  h_effMu_denom->Fill(lep_p4s->at(i).Pt() );

		 if( matchfound==true ) {
		   if( lep_flavor==11 )       h_effEE_numer->Fill(lep_p4s->at(i).Pt() );
		   else if( lep_flavor==13 )  h_effMu_numer->Fill(lep_p4s->at(i).Pt() );
		 }
	   }

   }


   ////////////////////////////// From here on, only proceed if we've got two good matches (i.e. a good reco dilepton)
   if( recoleps_matched.size() < 2 ) return;


   //Calculate dilepton invariant masses
   double dilMass_gen = -9.;
   double dilMass_reco = -9.;

   dilMass_gen  =  ( lep_p4s->at(0)         + lep_p4s->at(1)         ).M();
   dilMass_reco =  ( recoleps_matched.at(0) + recoleps_matched.at(1) ).M();

   //Fill dilepton mass histograms
   h_dilMass_gen->Fill( dilMass_gen );
   h_dilMass_reco->Fill( dilMass_reco );

   if( lep_flavor == 11) {
	 h_eeMass_gen->Fill( dilMass_gen );
	 h_eeMass_reco->Fill( dilMass_reco );
   }
   else if( lep_flavor == 13) {
	 h_mumuMass_gen->Fill( dilMass_gen );
	 h_mumuMass_reco->Fill( dilMass_reco );
   }

   //Fill residual plots
   h_residual_zmass->Fill( z_p4->mass(), z_p4->mass() - dilMass_reco );

   for( unsigned int i=0; i<2; i++) {
	 if( lep_flavor == 11 ) h_residual_elpt->Fill( lep_p4s->at(i).Pt(), lep_p4s->at(i).Pt() - recoleps_matched.at(i).Pt() );
	 if( lep_flavor == 13 ) h_residual_mupt->Fill( lep_p4s->at(i).Pt(), lep_p4s->at(i).Pt() - recoleps_matched.at(i).Pt() );
   }


}


// ------------ method called once each job just before starting event loop  ------------
void 
ZAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZAnalyzer::endJob() 
{
  using namespace std;

  //Make efficiency plots
  h_effEE->Divide( h_effEE_numer, h_effEE_denom, 1, 1, "B" );
  h_effMu->Divide( h_effMu_numer, h_effMu_denom, 1, 1, "B" );

  //Make resolution plots
  const bool olddir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);

  h_residual_zmass->FitSlicesY();
  TH1D* h_resolution_zmass = dynamic_cast<TH1D*>(gDirectory->Get("residzmass_2"));
  h_resolution_zmass->SetNameTitle( "resozmass", "Z mass resolution" );

  h_residual_elpt->FitSlicesY();
  TH1D* h_resolution_elpt = dynamic_cast<TH1D*>(gDirectory->Get("residelpt_2"));
  h_resolution_elpt->SetNameTitle( "resoelpt", "Electron pt resolution" );

  h_residual_mupt->FitSlicesY();
  TH1D* h_resolution_mupt = dynamic_cast<TH1D*>(gDirectory->Get("residmupt_2"));
  h_resolution_mupt->SetNameTitle( "resomupt", "Muon pt resolution" );

  TH1::AddDirectory(olddir);
}

// ------------ method called when starting to processes a run  ------------
void 
ZAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZAnalyzer);
