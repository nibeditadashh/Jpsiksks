// -*- C++ -*-
//
// Package:    JPsiKs0Ks0
// Class:      JPsiKs0Ks0
// 
//=================================================
// Original by: Chandiprasad Kar                  |
//<chandiprasad.kar@cern.ch>                      |
//Major variables are added 05 Sept, 2020    
// started on 11.10.21     |
//=================================================

// system include files
#include <memory>

// user include files
#include "jpsiksks/JPsiKsKsPAT/src/JPsiKs0Ks0.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include <vector>
#include <utility>
#include <string>
#include <iostream>
//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//
const static unsigned kNBR_CLOSED_TRACKS = 20;
const double PI = 3.141592653589793;
//
// constructors and destructor
//
JPsiKs0Ks0::JPsiKs0Ks0(const edm::ParameterSet& iConfig):

  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  //  trakCollection_label(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  triggerPrescales_ (consumes<pat::PackedTriggerPrescales> (iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),  
  trigTable_( iConfig.getParameter<std::vector<std::string> >("TriggerNames")),   
  puToken_(consumes<std::vector< PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PuInfoTag"))), 

  v0PtrCollection_(consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secondaryVerticesPtr"))),	       

  genParticles_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  doMC_ ( iConfig.getUntrackedParameter<bool>("doMC",false) ),
  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),
  dr0(0),dr1(0), dr01(0),dr11(0), //added for Ks01
  dpt0(0),dpt1(0), dpt01(0),dpt11(0), //added for Ks01
  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
  
  //************************ ****************************************************
  
  //*******************************************************
  nB(0), nMu(0),
  B_mass(0), B_px(0), B_py(0), B_pz(0),
  
  B_Ks0_mass(0), B_Ks0_px(0), B_Ks0_py(0), B_Ks0_pz(0),
  B_Ks0_pt1(0), B_Ks0_px1(0), B_Ks0_py1(0), B_Ks0_pz1(0), 
  B_Ks0_pt2(0), B_Ks0_px2(0), B_Ks0_py2(0), B_Ks0_pz2(0), 

  B_Ks0_px1_track(0), B_Ks0_py1_track(0), B_Ks0_pz1_track(0), 
  B_Ks0_px2_track(0), B_Ks0_py2_track(0), B_Ks0_pz2_track(0), 

  pi1dxy(0), pi2dxy(0), pi1dz(0), pi2dz(0),
  pi1dxy_e(0), pi2dxy_e(0), pi1dz_e(0), pi2dz_e(0),
  B_Ks0_charge1(0), B_Ks0_charge2(0),

 //added for Ks1
 
  B_Ks01_mass(0), B_Ks01_px(0), B_Ks01_py(0), B_Ks01_pz(0),
  B_Ks01_pt1(0), B_Ks01_px1(0), B_Ks01_py1(0), B_Ks01_pz1(0), 
  B_Ks01_pt2(0), B_Ks01_px2(0), B_Ks01_py2(0), B_Ks01_pz2(0), 

  B_Ks01_px1_track(0), B_Ks01_py1_track(0), B_Ks01_pz1_track(0), 
  B_Ks01_px2_track(0), B_Ks01_py2_track(0), B_Ks01_pz2_track(0), 
  
  pi11dxy(0), pi21dxy(0), pi11dz(0), pi21dz(0),
  pi11dxy_e(0), pi21dxy_e(0), pi11dz_e(0), pi21dz_e(0),
  B_Ks01_charge1(0), B_Ks01_charge2(0),

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0), 
  B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),

  B_Ks0_chi2(0), B_J_chi2(0), B_chi2(0), B_chi2dof(0), B_Ks01_chi2(0), //added Ks1
  B_Prob(0), B_J_Prob(0), B_ks0_Prob(0),
  bDecayVtxX(0), bDecayVtxY(0), bDecayVtxZ(0), bDecayVtxXE(0), bDecayVtxYE(0), bDecayVtxZE(0),
  bDecayVtxXYE(0), bDecayVtxXZE(0), bDecayVtxYZE(0),
  
  VDecayVtxX(0), VDecayVtxY(0), VDecayVtxZ(0), VDecayVtxXE(0), VDecayVtxYE(0), VDecayVtxZE(0),
  VDecayVtxXYE(0), VDecayVtxXZE(0), VDecayVtxYZE(0), 
  //added for Ks1
  VDecay1VtxX(0), VDecay1VtxY(0), VDecay1VtxZ(0), VDecay1VtxXE(0), VDecay1VtxYE(0), VDecay1VtxZE(0),
  VDecay1VtxXYE(0), VDecay1VtxXZE(0), VDecay1VtxYZE(0),
  pVtxIPX(0),  pVtxIPY(0),  pVtxIPZ(0),
  pVtxIPXE(0),  pVtxIPYE(0),  pVtxIPZE(0),  pVtxIPCL(0),
  pVtxIPXYE(0),  pVtxIPXZE(0),  pVtxIPYZE(0),
  
  B_l3d(0),  B_l3dE(0),  B_lxy(0), B_lxyE(0),
  B_cosalpha(0),   B_cosalphaxy(0), alpha(0),  B_treco(0),   B_trecoe(0),  B_trecoxy(0), B_trecoxye(0),
  B_pvip(0), B_pviperr(0), B_pvips(0), B_pvlzip(0), B_pvlziperr(0), B_pvlzips(0),
  B_pv2ip(0), B_pv2iperr(0), B_pv2ips(0), B_pv2lzip(0), B_pv2lziperr(0), B_pv2lzips(0),

  B_l3d_pv2(0),  B_l3dE_pv2(0),
  B_iso(0), B_mum_iso(0), B_mup_iso(0), B_pi1_iso(0),B_pi2_iso(0),B_pi3_iso(0),B_pi4_iso(0),
  
  istruemum(0), istruemup(0), istruekp(0), istruekm(0), istruekp1(0), istruekm1(0), istruebs(0),
  bunchXingMC(0), numInteractionsMC(0), trueNumInteractionsMC(0),
  run(0), event(0),
  lumiblock(0)
  
{
   //now do what ever initialization is needed
}

JPsiKs0Ks0::~JPsiKs0Ks0()
{

}

// ------------ method called to for each event  ------------
void JPsiKs0Ks0::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

  edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> theV0PtrHandle;
  iEvent.getByToken(v0PtrCollection_,  theV0PtrHandle);

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  //edm::Handle<std::vector<pat::GenericParticle> > thePATTrackHandle;  
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  fMagneticField = fieldHandle.product();

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genParticles_, pruned);

  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);
  if(isMC_ ){
    edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(puToken_, PupInfo);  
    for (std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); PVI++)
      {
	bunchXingMC->push_back(PVI->getBunchCrossing());
	numInteractionsMC->push_back(PVI->getPU_NumInteractions());
	trueNumInteractionsMC->push_back(PVI->getTrueNumInteractions());
      }
  }
  //std::cout << "Inside analyze " << std::endl;
  //***************************************
  // MC Gen information
  //***************************************
  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_ks_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion4_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_b_ct = -9999.;
  
  if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == 531) ) { //&& (dau->status() == 2) ) {
	foundit++;
	gen_b_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	gen_b_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
	for (size_t k=0; k<dau->numberOfDaughters(); k++) {
	  const reco::Candidate *gdau = dau->daughter(k);
	  if (gdau->pdgId()==443 ) { //&& gdau->status()==2) {
	    foundit++;
	    gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
	    gen_b_ct = GetLifetime(gen_b_p4,gen_b_vtx,gen_jpsi_vtx);
	    int nm=0;
	    for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
	      const reco::Candidate *mm = gdau->daughter(l);
	      if (mm->pdgId()==13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	      if (mm->pdgId()==-13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	    }
	    if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
	    else foundit-=nm;
	  }
	  for (size_t lk=0; lk<packed->size(); lk++) {
	    const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(0);
	    int stable_id = (*packed)[lk].pdgId();
	    if (dauInPrunedColl != nullptr && isAncestor(gdau,dauInPrunedColl)) {
	      if(stable_id == 211) {foundit++;
		gen_pion1_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	      }
	      if(stable_id == -211){ foundit++;
		gen_pion2_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	      }
	    }
	  }
	for (size_t lk=0; lk<packed->size(); lk++) {
	    const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(1);
	    int stable_id = (*packed)[lk].pdgId();
	    if (dauInPrunedColl != nullptr && isAncestor(gdau,dauInPrunedColl)) {
	      if(stable_id == 211) {foundit++;
		gen_pion3_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	      }
	      if(stable_id == -211){ foundit++;
		gen_pion4_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	      }
	    }
	  }

        } // for (size_t k                                                                                                                                                          
      }   // if (abs(dau->pdgId())==531 )                                                                                                                                           
      if (foundit>=9) break;
    } // for i                                                                                                                                                                      
    if (foundit!=7) {
      gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_ks_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_b_vtx.SetXYZ(0.,0.,0.);
      gen_jpsi_vtx.SetXYZ(0.,0.,0.);
      gen_b_ct = -9999.;
      //std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
  }

  nB = 0; nMu = 0;
  if ( OnlyGen_ ) { 
    tree_->Fill();
    return;
  }
  //std::cout << " found the given run " << iEvent.id().run() << "," << iEvent.id().event() << std::endl; // sanity check 

  //****************************************************************
  //Triggers

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales>            triggerPrescales;
  edm::Handle<edm::TriggerResults> hltTriggerResults;

  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  iEvent.getByToken(triggerResults_Label, hltTriggerResults ); 
 
  std::string cctk_0 = "hltJpsiTkAllConeTracksIter";
  std::string cctk_1 = "hltPsiPrimeTkAllConeTracksIter";
  std::string cctk_2 = "hltLowMassNonResonantTkAllConeTracksIter";
  std::stringstream myString;
  std::vector<float> obj_eta, obj_phi, obj_pt;
  std::vector<int> obj_charge;

  //bool foundOneTrig = false;
  const edm::TriggerNames &names = iEvent.triggerNames(*hltTriggerResults);
  for (unsigned int i = 0, n = hltTriggerResults->size(); i < n; ++i) {
    for (unsigned int it = 0; it < trigTable_.size(); it++){
      if (names.triggerName(i).find(trigTable_[it]) != std::string::npos && hltTriggerResults->accept(i))
	{ 
	  //triggernames->push_back(names.triggerName(i));
	  //triggerprescales->push_back(triggerPrescales->getPrescaleForIndex(i));
	  // cout<<"Trigger name "<<names.triggerName(i)<<" Prescale "<<triggerPrescales->getPrescaleForIndex(i)<<endl;
	  //foundOneTrig = true;
	}
    }
  }
  //if ( iEvent.isRealData() && !foundOneTrig) return;
     
  if (triggerObjects.isValid()) {
    //std::cout << "will try to match trigger object with track " << triggerObjects->size() << endl;    
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
      obj.unpackFilterLabels(iEvent,*hltTriggerResults);
      std::string cc1 = obj.collection();
      
      for (unsigned int i = 0; i < trigTable_.size(); i++)
        {
	  myString.clear(); myString.str(""); myString << trigTable_[i] << "*";
	  if ( obj.hasPathName(myString.str().c_str(), true, true) )
            {

	      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	      std::vector<std::string> pathNamesLast = obj.pathNames(true);
	            
	      //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
	            
	      if ( (cc1.find(cctk_0) != std::string::npos ) ){//|| (cc1.find(cctk_1) != std::string::npos) || ( cc1.find(cctk_2) != std::string::npos)) {
		
		obj_eta.push_back(obj.eta());
		obj_phi.push_back(obj.phi());
		obj_pt.push_back(obj.pt());
		obj_charge.push_back(obj.charge());
		
		// std::cout << myString.str().c_str() << std::endl;
		// std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
		// std::cout << "\t   Charge:   "<<obj.charge()<<std::endl;
		// std::cout << "\t   Collection: " << obj.collection() << std::endl;
		// std::cout << "\t   Type IDs:   ";
		// for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
		// std::cout << "\t   Filters:    ";
		// for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
		// std::cout << std::endl;
	      }
	    }
	}
    }
  }
    

  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  reco::Vertex bestVtxBS;


  // get primary vertex
  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_Label, recVtxs);
  edm::Handle<reco::VertexCollection> pvHandle_;
  iEvent.getByToken(primaryVertices_Label,pvHandle_);
  bestVtx = *(recVtxs->begin());
  
  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);
  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 

  nVtx = recVtxs->size();

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event(); 
  
  //std::cout<<"After Primary vertex check "<< std::endl;
  //Let's begin by looking for J/psi->mu+mu-
  std::cout<<"Patmuon size  "<<thePATMuonHandle->size()<< std::endl;
  unsigned int nMu_tmp = thePATMuonHandle->size();
 
  for(edm::View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {      
      for(edm::View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue; // <-------------------------------------
	  
	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) {glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }
	  
	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;
	  
	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
	  
	  //Let's check the vertex and mass
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	  // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************

	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;

	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }

	  KinematicParticleVertexFitter fitter;   

	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }

	  if (!psiVertexFitTree->isValid()) 
	    {
	      std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }

	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();//masa del J/psi
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();//vertice del J/psi
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }

	  double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(J_Prob_tmp<0.01)
	    {
	      continue;
	    }
	  
	   //some loose cuts go here
	  
	  if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	  if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;
	  
	  //  ***************  added two Ks0 Ks01 loop
	  std::cout<<" code is working till here 1 "<<std::endl;    
           if ( theV0PtrHandle->size()>2 && thePATMuonHandle->size()>=2 )
            { 
	      for ( vector<VertexCompositePtrCandidate>::const_iterator iVee = theV0PtrHandle->begin();   iVee != theV0PtrHandle->end(); ++iVee )
		{

		for ( vector<VertexCompositePtrCandidate>::const_iterator jVee = iVee+1;   iVee != theV0PtrHandle->end(); ++jVee )
		{
std::cout<<" code is working till here 2 "<<std::endl;
		//added for 2 ks01
		 if(iVee==jVee) continue;
	  
	  

		  //get the tracks from V0 candidate
		  vector<pat::PackedCandidate> v0daughters;
		  vector<Track> theDaughterTracks;
		  v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(0))) );
		  v0daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(iVee->daughter(1))) );

                  //added for Ks01
		  vector<pat::PackedCandidate> v01daughters;
		  vector<Track> theDaughter1Tracks;
		  v01daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(jVee->daughter(0))) );
		  v01daughters.push_back( *(dynamic_cast<const pat::PackedCandidate *>(jVee->daughter(1))) );
		  
		  for(unsigned int j = 0; j < v0daughters.size(); ++j)
		    {
		      theDaughterTracks.push_back(v0daughters[j].pseudoTrack());
		    }

 		//added for Ks01
     		for(unsigned int t = 0; t < v01daughters.size(); ++t)
		    {
		      theDaughter1Tracks.push_back(v01daughters[t].pseudoTrack());
		    }
		  std::cout<<" code is working till here 3 "<<std::endl;
		   // it does not have sences here. 
		   //if ( IsTheSame(*theDaughterTracks[0],*iMuon1) || IsTheSame(*theDaughterTracks[0],*iMuon2) ) continue;
		   //if ( IsTheSame(*theDaughterTracks[1],*iMuon1) || IsTheSame(*theDaughterTracks[1],*iMuon2) ) continue;
		  
		   //Now let's see if these two tracks make a vertex for Ks0
		   reco::TransientTrack pion1TT((*theB).build(theDaughterTracks[0]));
		   reco::TransientTrack pion2TT((*theB).build(theDaughterTracks[1]));

		// //Now let's see if these two tracks make a vertex for Ks01
			
		   reco::TransientTrack pion3TT((*theB).build(theDaughter1Tracks[0]));
		   reco::TransientTrack pion4TT((*theB).build(theDaughter1Tracks[1]));	     
		   
		   ParticleMass pion_mass = 0.13957018;
		   ParticleMass Ks0_mass = 0.497614;
		   float pion_sigma = pion_mass*1.e-6;
		   float Ks0_sigma = Ks0_mass*1.e-6;
		   
		   //initial chi2 and ndf before kinematic fits.
		   float chi = 0.;
		   float ndf = 0.;
		   vector<RefCountedKinematicParticle> pionParticles;
                   vector<RefCountedKinematicParticle> pion1Particles;//Ks01
		   // vector<RefCountedKinematicParticle> muonParticles;
		   try {
		     pionParticles.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
		     pionParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));
//added for Ks01
pion1Particles.push_back(pFactory.particle(pion3TT,pion_mass,chi,ndf,pion_sigma));
		     pion1Particles.push_back(pFactory.particle(pion4TT,pion_mass,chi,ndf,pion_sigma));
		   }
		   catch(...) {
		     std::cout<<" Exception caught ... continuing 3 "<<std::endl;
		     continue;
		   }
		   
		   RefCountedKinematicTree Ks0VertexFitTree;
                   RefCountedKinematicTree Ks01VertexFitTree;//Ks01
		   try{
		     Ks0VertexFitTree = fitter.fit(pionParticles); 
                     Ks01VertexFitTree = fitter.fit(pion1Particles); 
		   }
		   catch(...) {
		     std::cout<<" Exception caught ... continuing 4 "<<std::endl;                   
		     continue;
		   }
		   if (!Ks0VertexFitTree->isValid() && !Ks01VertexFitTree->isValid()) //added Ks01 info 
		     {
		       std::cout << "invalid vertex from the Ks0 && Ks01 vertex fit" << std::endl;
		       continue; 
		     }
		   Ks0VertexFitTree->movePointerToTheTop();
		   Ks01VertexFitTree->movePointerToTheTop();//added Ks1
		   
		   RefCountedKinematicParticle Ks0_vFit_noMC = Ks0VertexFitTree->currentParticle();
		   RefCountedKinematicVertex Ks0_vFit_vertex_noMC = Ks0VertexFitTree->currentDecayVertex();
		   //added Ks1
                   RefCountedKinematicParticle Ks01_vFit_noMC = Ks01VertexFitTree->currentParticle();
		   RefCountedKinematicVertex Ks01_vFit_vertex_noMC = Ks01VertexFitTree->currentDecayVertex();
		  
		   if( Ks0_vFit_vertex_noMC->chiSquared() < 0 && Ks01_vFit_vertex_noMC->chiSquared() < 0 ) //added Ks1
		     { 
		       std::cout << "negative chisq from ks fit" << endl;
		       continue;
		     }
		  std::cout<<" code is working till here 31 "<<std::endl;
		   //some loose cuts go here
		   
		   if(Ks0_vFit_vertex_noMC->chiSquared()>50) continue;
		   if(Ks01_vFit_vertex_noMC->chiSquared()>50) continue;//added for Ks1
		   if(Ks0_vFit_noMC->currentState().mass()<0.45 || Ks0_vFit_noMC->currentState().mass()>0.55) continue;
		   if(Ks01_vFit_noMC->currentState().mass()<0.45 || Ks01_vFit_noMC->currentState().mass()>0.55) continue;//added for Ks1

		   Ks0VertexFitTree->movePointerToTheFirstChild();
		   RefCountedKinematicParticle T1CandMC = Ks0VertexFitTree->currentParticle();
		   
		   Ks0VertexFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle T2CandMC = Ks0VertexFitTree->currentParticle();
		   

		//added for Ks1
                 
 		   Ks01VertexFitTree->movePointerToTheFirstChild();
		   RefCountedKinematicParticle T3CandMC = Ks01VertexFitTree->currentParticle();
		   
		   Ks01VertexFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle T4CandMC = Ks01VertexFitTree->currentParticle();
		   

		   //  Ks0  mass constrain
		   // do mass constrained vertex fit
		   // creating the constraint with a small sigma to put in the resulting covariance 
		   // matrix in order to avoid singularities
		   // JPsi mass constraint is applied in the final B fit
		   
		   KinematicParticleFitter csFitterKs;
		   KinematicConstraint * ks_c = new MassKinematicConstraint(Ks0_mass,Ks0_sigma);
		   // add mass constraint to the ks0 fit to do a constrained fit:  
		   
		   Ks0VertexFitTree = csFitterKs.fit(ks_c,Ks0VertexFitTree);
		   if (!Ks0VertexFitTree->isValid()){
		     std::cout << "caught an exception in the ks mass constraint fit" << std::endl;
		     continue; 
		   }
		   
		   Ks0VertexFitTree->movePointerToTheTop();
		   RefCountedKinematicParticle ks0_vFit_withMC = Ks0VertexFitTree->currentParticle();
		   
		//added Ks1
		   KinematicParticleFitter csFitterKs1;
		   KinematicConstraint * ks1_c = new MassKinematicConstraint(Ks0_mass,Ks0_sigma);
		   // add mass constraint to the ks0 fit to do a constrained fit:  
		   
		   Ks01VertexFitTree = csFitterKs1.fit(ks1_c,Ks01VertexFitTree);
		   if (!Ks01VertexFitTree->isValid()){
		     std::cout << "caught an exception in the ks1 mass constraint fit" << std::endl;
		     continue; 
		   }
		  std::cout<<" code is working till here 32 "<<std::endl; 
		   Ks01VertexFitTree->movePointerToTheTop();
		   RefCountedKinematicParticle ks01_vFit_withMC = Ks01VertexFitTree->currentParticle();

		   //Now we are ready to combine!
		   // JPsi mass constraint is applied in the final Bd fit,
		   
		   vector<RefCountedKinematicParticle> vFitMCParticles;
		   vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(ks0_vFit_withMC);
                   vFitMCParticles.push_back(ks01_vFit_withMC); //added for ks1
std::cout<<" code is working till here 33 "<<std::endl;
		   MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
		   KinematicConstrainedVertexFitter kcvFitter;
		   RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
		   if (!vertexFitTree->isValid()) {
		     std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		    continue;
		   }
		std::cout<<" code is working till here 4 "<<std::endl;
		   vertexFitTree->movePointerToTheTop();		     
		    std::cout<<" code is working till here 41 "<<std::endl;
		   RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		   RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
		   if (!bDecayVertexMC->vertexIsValid()){
		     std::cout << "B MC fit vertex is not valid" << endl;
		     continue;
		   }
		   std::cout << "B mass "<<bCandMC->currentState().mass()<< std::endl;
		   if(bCandMC->currentState().mass()<5.0 || bCandMC->currentState().mass()>6.0) continue;
		   
		   if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 ) 
		     {
		       std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
		       continue;
		     }
		   
		   double B_Prob_tmp  = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   if(B_Prob_tmp<0.01)
		     {
		       continue;
		     }		     
		   
		   // get children from final B fit
		   vertexFitTree->movePointerToTheFirstChild();
		   RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
		   vertexFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
		   
		   vertexFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle Ks0CandMC = vertexFitTree->currentParticle();

		//ADDED for Ks1
                   vertexFitTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle Ks01CandMC = vertexFitTree->currentParticle();
		
		   
		   KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMupKP;
		   KinematicParameters psiMumKP;
	       
		   if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
		   if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
		   if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
		   if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 

 		   GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
		   		       mu1CandMC->currentState().globalMomentum().y(),
 		   		       mu1CandMC->currentState().globalMomentum().z());

 		   GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
		   		       mu2CandMC->currentState().globalMomentum().y(),
 		   		       mu2CandMC->currentState().globalMomentum().z());

 		   GlobalVector Ks0p1vec(T1CandMC->currentState().globalMomentum().x(),
		   		        T1CandMC->currentState().globalMomentum().y(),
 		   		        T1CandMC->currentState().globalMomentum().z());

 		   GlobalVector Ks0p2vec(T2CandMC->currentState().globalMomentum().x(),
		   			T2CandMC->currentState().globalMomentum().y(),
		   			T2CandMC->currentState().globalMomentum().z());

		// added for Ks1

		   GlobalVector Ks01p1vec(T3CandMC->currentState().globalMomentum().x(),
		   		        T3CandMC->currentState().globalMomentum().y(),
 		   		        T3CandMC->currentState().globalMomentum().z());

 		   GlobalVector Ks01p2vec(T4CandMC->currentState().globalMomentum().x(),
		   			T4CandMC->currentState().globalMomentum().y(),
		   			T4CandMC->currentState().globalMomentum().z());


		   KinematicParameters Ks0Pi1KP = T1CandMC->currentState().kinematicParameters();
		   KinematicParameters Ks0Pi2KP = T2CandMC->currentState().kinematicParameters();
		   KinematicParameters Ks0PipKP;
		   KinematicParameters Ks0PimKP;
	       
		   if ( T1CandMC->currentState().particleCharge() > 0 ) Ks0PipKP = Ks0Pi1KP;
		   if ( T1CandMC->currentState().particleCharge() < 0 ) Ks0PimKP = Ks0Pi1KP;
		   if ( T2CandMC->currentState().particleCharge() > 0 ) Ks0PipKP = Ks0Pi2KP;
		   if ( T2CandMC->currentState().particleCharge() < 0 ) Ks0PimKP = Ks0Pi2KP;

		//added for KS1

		   KinematicParameters Ks01Pi1KP = T3CandMC->currentState().kinematicParameters();
		   KinematicParameters Ks01Pi2KP = T4CandMC->currentState().kinematicParameters();
		   KinematicParameters Ks01PipKP;
		   KinematicParameters Ks01PimKP;
	       
		   if ( T3CandMC->currentState().particleCharge() > 0 ) Ks01PipKP = Ks01Pi1KP;
		   if ( T3CandMC->currentState().particleCharge() < 0 ) Ks01PimKP = Ks01Pi1KP;
		   if ( T4CandMC->currentState().particleCharge() > 0 ) Ks01PipKP = Ks01Pi2KP;
		   if ( T4CandMC->currentState().particleCharge() < 0 ) Ks01PimKP = Ks01Pi2KP;	 
	 std::cout<<" code is working till here 5 "<<std::endl;

		   // fill candidate variables now
		   
		   if(nB==0){		    
		     nMu  = nMu_tmp;
		     cout<< "*Number of Muons : " << nMu_tmp << endl;
		   } // end nB==0		     
		
		   B_mass->push_back(bCandMC->currentState().mass());
		   B_px->push_back(bCandMC->currentState().globalMomentum().x());
		   B_py->push_back(bCandMC->currentState().globalMomentum().y());
		   B_pz->push_back(bCandMC->currentState().globalMomentum().z());
		
		   B_Ks0_mass->push_back( Ks0_vFit_noMC->currentState().mass() );
		   B_Ks0_px->push_back( Ks0_vFit_noMC->currentState().globalMomentum().x() );
		   B_Ks0_py->push_back( Ks0_vFit_noMC->currentState().globalMomentum().y() );
		   B_Ks0_pz->push_back( Ks0_vFit_noMC->currentState().globalMomentum().z() );

		//added for Ks1

		   B_Ks01_mass->push_back( Ks01_vFit_noMC->currentState().mass() );
		   B_Ks01_px->push_back( Ks01_vFit_noMC->currentState().globalMomentum().x() );
		   B_Ks01_py->push_back( Ks01_vFit_noMC->currentState().globalMomentum().y() );
		   B_Ks01_pz->push_back( Ks01_vFit_noMC->currentState().globalMomentum().z() );

		   B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
		   B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		   B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		   B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );

		   B_Ks0_pt1->push_back(Ks0p1vec.perp());
		   B_Ks0_px1->push_back(Ks0Pi1KP.momentum().x());
		   B_Ks0_py1->push_back(Ks0Pi1KP.momentum().y());
		   B_Ks0_pz1->push_back(Ks0Pi1KP.momentum().z());
		   B_Ks0_px1_track->push_back(v0daughters[0].px());
		   B_Ks0_py1_track->push_back(v0daughters[0].py());
		   B_Ks0_pz1_track->push_back(v0daughters[0].pz());
		   B_Ks0_charge1->push_back(T1CandMC->currentState().particleCharge());

		   B_Ks0_pt2->push_back(Ks0p2vec.perp());
		   B_Ks0_px2->push_back(Ks0Pi2KP.momentum().x());
		   B_Ks0_py2->push_back(Ks0Pi2KP.momentum().y());
		   B_Ks0_pz2->push_back(Ks0Pi2KP.momentum().z());
		   B_Ks0_px2_track->push_back(v0daughters[1].px());
		   B_Ks0_py2_track->push_back(v0daughters[1].py());
		   B_Ks0_pz2_track->push_back(v0daughters[1].pz());
		   B_Ks0_charge2->push_back(T2CandMC->currentState().particleCharge());

		// added for Ks1

		   B_Ks01_pt1->push_back(Ks01p1vec.perp());
		   B_Ks01_px1->push_back(Ks01Pi1KP.momentum().x());
		   B_Ks01_py1->push_back(Ks01Pi1KP.momentum().y());
		   B_Ks01_pz1->push_back(Ks01Pi1KP.momentum().z());
		   B_Ks01_px1_track->push_back(v01daughters[0].px());
		   B_Ks01_py1_track->push_back(v01daughters[0].py());
		   B_Ks01_pz1_track->push_back(v01daughters[0].pz());
		   B_Ks01_charge1->push_back(T3CandMC->currentState().particleCharge());

		   B_Ks01_pt2->push_back(Ks01p2vec.perp());
		   B_Ks01_px2->push_back(Ks01Pi2KP.momentum().x());
		   B_Ks01_py2->push_back(Ks01Pi2KP.momentum().y());
		   B_Ks01_pz2->push_back(Ks01Pi2KP.momentum().z());
		   B_Ks01_px2_track->push_back(v01daughters[1].px());
		   B_Ks01_py2_track->push_back(v01daughters[1].py());
		   B_Ks01_pz2_track->push_back(v01daughters[1].pz());
		   B_Ks01_charge2->push_back(T4CandMC->currentState().particleCharge());


		   B_J_pt1->push_back(Jp1vec.perp());
		   B_J_px1->push_back(psiMu1KP.momentum().x());
		   B_J_py1->push_back(psiMu1KP.momentum().y());
		   B_J_pz1->push_back(psiMu1KP.momentum().z());
		   B_J_charge1->push_back(mu1CandMC->currentState().particleCharge());

		   B_J_pt2->push_back(Jp2vec.perp());
		   B_J_px2->push_back(psiMu2KP.momentum().x());
		   B_J_py2->push_back(psiMu2KP.momentum().y());
		   B_J_pz2->push_back(psiMu2KP.momentum().z());
		   B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());

		   B_Ks0_chi2->push_back(Ks0_vFit_vertex_noMC->chiSquared());
		   B_Ks01_chi2->push_back(Ks01_vFit_vertex_noMC->chiSquared());//added for ks1
		   B_J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
		   B_chi2->push_back(bDecayVertexMC->chiSquared());
		   B_chi2dof->push_back(bDecayVertexMC->chiSquared()/bDecayVertexMC->degreesOfFreedom());
		   
		   //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		 //  double ks0_Prob_tmp  = TMath::Prob(Ks0_vFit_vertex_noMC->chiSquared() + Ks01_vFit_vertex_noMC->chiSquared() ,(int)Ks0_vFit_vertex_noMC->degreesOfFreedom()+Ks01_vFit_vertex_noMC->degreesOfFreedom());
		 double ks0_Prob_tmp  = TMath::Prob(Ks0_vFit_vertex_noMC->chiSquared(),(int)Ks0_vFit_vertex_noMC->degreesOfFreedom());
		   B_Prob    ->push_back(B_Prob_tmp);
		   B_J_Prob  ->push_back(J_Prob_tmp);
		   B_ks0_Prob ->push_back(ks0_Prob_tmp);

		   // ************
		   bDecayVtxX->push_back((*bDecayVertexMC).position().x());
		   bDecayVtxY->push_back((*bDecayVertexMC).position().y());
		   bDecayVtxZ->push_back((*bDecayVertexMC).position().z());
		   bDecayVtxXE->push_back(bDecayVertexMC->error().cxx());
		   bDecayVtxYE->push_back(bDecayVertexMC->error().cyy());
		   bDecayVtxZE->push_back(bDecayVertexMC->error().czz());
		   bDecayVtxXYE->push_back(bDecayVertexMC->error().cyx());
		   bDecayVtxXZE->push_back(bDecayVertexMC->error().czx());
		   bDecayVtxYZE->push_back(bDecayVertexMC->error().czy());

		   VDecayVtxX->push_back( Ks0_vFit_vertex_noMC->position().x() );
		   VDecayVtxY->push_back( Ks0_vFit_vertex_noMC->position().y() );
		   VDecayVtxZ->push_back( Ks0_vFit_vertex_noMC->position().z() );
		   VDecayVtxXE->push_back( Ks0_vFit_vertex_noMC->error().cxx() );
		   VDecayVtxYE->push_back( Ks0_vFit_vertex_noMC->error().cyy() );
		   VDecayVtxZE->push_back( Ks0_vFit_vertex_noMC->error().czz() );
		   VDecayVtxXYE->push_back( Ks0_vFit_vertex_noMC->error().cyx() );
		   VDecayVtxXZE->push_back( Ks0_vFit_vertex_noMC->error().czx() );
		   VDecayVtxYZE->push_back( Ks0_vFit_vertex_noMC->error().czy() );

		  //added for Ks1

		   VDecay1VtxX->push_back( Ks01_vFit_vertex_noMC->position().x() );
		   VDecay1VtxY->push_back( Ks01_vFit_vertex_noMC->position().y() );
		   VDecay1VtxZ->push_back( Ks01_vFit_vertex_noMC->position().z() );
		   VDecay1VtxXE->push_back( Ks01_vFit_vertex_noMC->error().cxx() );
		   VDecay1VtxYE->push_back( Ks01_vFit_vertex_noMC->error().cyy() );
		   VDecay1VtxZE->push_back( Ks01_vFit_vertex_noMC->error().czz() );
		   VDecay1VtxXYE->push_back( Ks01_vFit_vertex_noMC->error().cyx() );
		   VDecay1VtxXZE->push_back( Ks01_vFit_vertex_noMC->error().czx() );
		   VDecay1VtxYZE->push_back( Ks01_vFit_vertex_noMC->error().czy() );
std::cout<<" code is working till here 6 "<<std::endl;
		   // ********************* muon-trigger-machint**************** 
		   
		   const pat::Muon* muon1 = &(*iMuon1);
		   const pat::Muon* muon2 = &(*iMuon2);

		   int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;
		   
		   if (muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) tri_Dim25_tmp = 1;
		   if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
		   if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr) tri_JpsiTkTk_tmp = 1;
		   
		   std::cout<<"Do the trigger matching "<<std::endl;
		   tri_Dim25->push_back( tri_Dim25_tmp );	       
		   tri_JpsiTk->push_back( tri_JpsiTk_tmp );
                   tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );
		   // Match the track with the trigger object
		   float dr0_t = 99999.;
		   float dpt0_t = 99999.;
		   for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		     float dp = pion1TT.track().phi() - obj_phi[ii];
		     float de = pion1TT.track().eta() - obj_eta[ii];      
		     if (dp>float(M_PI)) dp-=float(2*M_PI);  
		     float dr = std::sqrt(de*de + dp*dp);

		     if (dr < dr0_t) dr0_t = dr;
		     float dpt = pion1TT.track().pt() - obj_pt[ii];
		     if (abs(dpt) < dpt0_t) dpt0_t = abs(dpt);
		   }

		   float dr1_t = 99999.;
		   float dpt1_t = 99999.;
		   for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		     float dp = pion2TT.track().phi() - obj_phi[ii];
		     float de = pion2TT.track().eta() - obj_eta[ii];    
		     if (dp>float(M_PI)) dp-=float(2*M_PI);  
		     float dr = std::sqrt(de*de + dp*dp);		     
		     if (dr < dr1_t) dr1_t = dr;
		     float dpt = pion2TT.track().pt() - obj_pt[ii];
		     if (abs(dpt) < dpt1_t) dpt1_t =abs(dpt);

		   }

		// added for Ks01
                  float dr01_t = 99999.;
		   float dpt01_t = 99999.;
		   for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		     float dp = pion3TT.track().phi() - obj_phi[ii];
		     float de = pion3TT.track().eta() - obj_eta[ii];      
		     if (dp>float(M_PI)) dp-=float(2*M_PI);  
		     float dr = std::sqrt(de*de + dp*dp);

		     if (dr < dr01_t) dr01_t = dr;
		     float dpt = pion3TT.track().pt() - obj_pt[ii];
		     if (abs(dpt) < dpt01_t) dpt01_t = abs(dpt);
		   }

		   float dr11_t = 99999.;
		   float dpt11_t = 99999.;
		   for (uint ii=0 ; ii < obj_eta.size(); ii++) {
		     float dp = pion4TT.track().phi() - obj_phi[ii];
		     float de = pion4TT.track().eta() - obj_eta[ii];    
		     if (dp>float(M_PI)) dp-=float(2*M_PI);  
		     float dr = std::sqrt(de*de + dp*dp);		     
		     if (dr < dr11_t) dr11_t = dr;
		     float dpt = pion4TT.track().pt() - obj_pt[ii];
		     if (abs(dpt) < dpt11_t) dpt11_t =abs(dpt);

		   }

		   
		   dr0->push_back(dr0_t);
		   dr1->push_back(dr1_t);

		   dpt0->push_back(dpt0_t);
		   dpt1->push_back(dpt1_t);

		//added for Ks01
		   dr01->push_back(dr01_t);
		   dr11->push_back(dr11_t);

		   dpt01->push_back(dpt01_t);
		   dpt11->push_back(dpt11_t);

		   // ************
		  
		   mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
		   mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
		   mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
		   mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
		   mu1PF->push_back(iMuon1->isPFMuon());
		   mu2PF->push_back(iMuon2->isPFMuon());
		   mu1loose->push_back(muon::isLooseMuon(*iMuon1));
		   mu2loose->push_back(muon::isLooseMuon(*iMuon2));

		   mumC2->push_back( glbTrackM->normalizedChi2() );
		   mumNHits->push_back( glbTrackM->numberOfValidHits() );
		   mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
		   mupC2->push_back( glbTrackP->normalizedChi2() );
		   mupNHits->push_back( glbTrackP->numberOfValidHits() );
		   mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
                   mumdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
		   mupdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
		   mumdz->push_back(glbTrackM->dz(bestVtx.position()) );
		   mupdz->push_back(glbTrackP->dz(bestVtx.position()) );
		   muon_dca->push_back(dca);

		   pi1dxy->push_back(v0daughters[0].dxy());
		   pi2dxy->push_back(v0daughters[1].dxy());
		   pi1dz->push_back(v0daughters[0].dz());
		   pi2dz->push_back(v0daughters[1].dz());

		   pi1dxy_e->push_back(v0daughters[0].dxyError());
		   pi2dxy_e->push_back(v0daughters[1].dxyError());
		   pi1dz_e->push_back(v0daughters[0].dzError());
		   pi2dz_e->push_back(v0daughters[1].dzError());

//added for Ks01

 		   pi11dxy->push_back(v01daughters[0].dxy());
		   pi21dxy->push_back(v01daughters[1].dxy());
		   pi11dz->push_back(v01daughters[0].dz());
		   pi21dz->push_back(v01daughters[1].dz());

		   pi11dxy_e->push_back(v01daughters[0].dxyError());
		   pi21dxy_e->push_back(v01daughters[1].dxyError());
		   pi11dz_e->push_back(v01daughters[0].dzError());
		   pi21dz_e->push_back(v01daughters[1].dzError());
		   
		   // // ********************* loop over all the primary vertices and we choose the one with the best pointing angle **************** 
		   // // reco::Vertex bestVtxIP;

		   // // Double_t pVtxIPX_temp = -10000.0;
		   // // Double_t pVtxIPY_temp = -10000.0;
		   // // Double_t pVtxIPZ_temp = -10000.0;
		   // // Double_t pVtxIPXE_temp = -10000.0;
		   // // Double_t pVtxIPYE_temp = -10000.0;
		   // // Double_t pVtxIPZE_temp = -10000.0;
		   // // Double_t pVtxIPXYE_temp = -10000.0;
		   // // Double_t pVtxIPXZE_temp = -10000.0;
		   // // Double_t pVtxIPYZE_temp = -10000.0;
		   // // Double_t pVtxIPCL_temp = -10000.0;
		   // // Double_t lip1 = -1000000.0;
		   // // for(size_t i = 0; i < recVtxs->size(); ++i) {
		   // //   const Vertex &vtx = (*recVtxs)[i];
		            
		   // //   Double_t dx1 = (*bDecayVertexMC).position().x() - vtx.x(); 
		   // //   Double_t dy1 = (*bDecayVertexMC).position().y() - vtx.y();
		   // //   Double_t dz1 = (*bDecayVertexMC).position().z() - vtx.z();
		   // //   float cosAlphaXYb1 = ( bCandMC->currentState().globalMomentum().x() * dx1 + bCandMC->currentState().globalMomentum().y()*dy1 + bCandMC->currentState().globalMomentum().z()*dz1  )/( sqrt(dx1*dx1+dy1*dy1+dz1*dz1)* bCandMC->currentState().globalMomentum().mag() );

		   // //   if(cosAlphaXYb1>lip1)
		   // //     {
		   // // 	 lip1 = cosAlphaXYb1 ;
		   // // 	 pVtxIPX_temp = vtx.x();
		   // // 	 pVtxIPY_temp = vtx.y();
		   // // 	 pVtxIPZ_temp = vtx.z();
		   // // 	 pVtxIPXE_temp = vtx.covariance(0, 0);
		   // // 	 pVtxIPYE_temp = vtx.covariance(1, 1);
		   // // 	 pVtxIPZE_temp = vtx.covariance(2, 2);
		   // // 	 pVtxIPXYE_temp = vtx.covariance(0, 1);
		   // // 	 pVtxIPXZE_temp = vtx.covariance(0, 2);
		   // // 	 pVtxIPYZE_temp = vtx.covariance(1, 2);
		   // // 	 pVtxIPCL_temp = (TMath::Prob(vtx.chi2(),(int)vtx.ndof()) );

		   // // 	 bestVtxIP = vtx;
			  
		   // //     }
                 
		   // // }

		   // // // // try refitting the primary without the tracks in the B reco candidate		   
		   // // // const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(iMuon1->originalObject());
		   // // // const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(iMuon2->originalObject());
		   // // // reco::TrackRef patTrack1_1 = iTrack1->track();
		        
		   // // // // first get tracks from the original primary
		   // // // vector<reco::TransientTrack> vertexTracks;
		        
		   // // // for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestVtxIP.tracks_begin();
		   // // // 	 iTrack != bestVtxIP.tracks_end(); ++iTrack) {
		   // // //   // compare primary tracks to check for matches with B cand
		   // // //   reco::TrackRef trackRef = iTrack->castTo<TrackRef>();
		            
		   // // //   // the 3 tracks in the Bs candidate are  patTrack1_1 rmu1 and rmu2 
		   // // //   if (  !( (patTrack1_1.key()==trackRef.key()) || (rmu1->track().key()==trackRef.key()) || (rmu2->track().key()==trackRef.key()) ) ) {		        
		   // // //     reco::TransientTrack tt((*theB).build(trackRef));
		   // // //     vertexTracks.push_back(tt);
		   // // //   }//else { std::cout << "found track match with primary" << endl;}
		   // // // }
		        
		   // // // // *** if no tracks in primary or no reco track included in primary then don't do anything ***
		   // // // reco::Vertex bestVtxRf = bestVtxIP;
		   // // // GlobalPoint PVRfP = GlobalPoint( bestVtxIP.x(), bestVtxIP.y(), bestVtxIP.z() );
		        
		   // // // if (  vertexTracks.size()>0 && (bestVtxIP.tracksSize()!=vertexTracks.size()) ) {
		   // // //   AdaptiveVertexFitter theFitter;
		   // // //   TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);
		   // // //   if ( v.isValid() ) {    
		   // // //     //set bestVtxRf as new best vertex to fill variables for refitting PV
		   // // //     bestVtxRf = reco::Vertex(v);
		        
		   // // //   }
		   // // // }
		   		   
		   // // // priRfVtxX->push_back( bestVtxRf.x() );
		   // // // priRfVtxY->push_back( bestVtxRf.y() );
		   // // // priRfVtxZ->push_back( bestVtxRf.z() );
		   // // // priRfVtxXE->push_back( bestVtxRf.covariance(0, 0) );
		   // // // priRfVtxYE->push_back( bestVtxRf.covariance(1, 1) );
		   // // // priRfVtxZE->push_back( bestVtxRf.covariance(2, 2) );
		   // // // priRfVtxXYE->push_back( bestVtxRf.covariance(0, 1) );
		   // // // priRfVtxXZE->push_back( bestVtxRf.covariance(0, 2) );
		   // // // priRfVtxYZE->push_back( bestVtxRf.covariance(1, 2) );  
		   // // // priRfVtxCL->push_back( ChiSquaredProbability((double)(bestVtxRf.chi2()),(double)(bestVtxRf.ndof())) );
		   // // // priRfNTrkDif->push_back( bestVtxIP.tracksSize() - vertexTracks.size() );
		   

		   nB++;	
		   
		   int pvIndex = -1;
		   saveIP(vertexFitTree,*pvHandle_.product(), pvIndex);
		   std::cout<<" pvindex "<<pvIndex<<std::endl;
		   SaveIso(vertexFitTree, fMagneticField, theB, pvIndex,
		   	   thePATTrackHandle, thePATMuonHandle,
		   	   muon1TT, muon2TT,pion1TT, pion2TT, pion3TT, pion4TT);


		   if(isMC_)saveTruthMatch(iEvent);

		   /////////////////////////////////////////////////
		   pionParticles.clear();
 		   pion1Particles.clear();
		   muonParticles.clear();
		   vFitMCParticles.clear();
		   }//V0
		   
		}//Ks01
	     }//Ks0
	}//mu2
    }//mu1    
//}
//}
   
   //fill the tree and clear the vectors
   if (nB > 0 ) 
     {
       //std::cout << "filling tree" << endl;
       tree_->Fill();
     }
   // *********

   nB = 0; nMu = 0;

   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();
   B_Ks0_mass->clear(); B_Ks0_px->clear(); B_Ks0_py->clear(); B_Ks0_pz->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();

   B_Ks0_pt1->clear(); B_Ks0_px1->clear(); B_Ks0_py1->clear(); B_Ks0_pz1->clear(); B_Ks0_charge1->clear(); 
   B_Ks0_pt2->clear(); B_Ks0_px2->clear(); B_Ks0_py2->clear(); B_Ks0_pz2->clear(); B_Ks0_charge2->clear(); 

   B_Ks0_px1_track->clear(); B_Ks0_py1_track->clear(); B_Ks0_pz1_track->clear(); 
   B_Ks0_px2_track->clear(); B_Ks0_py2_track->clear(); B_Ks0_pz2_track->clear(); 

   B_J_pt1->clear();  B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_pt2->clear();  B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   B_Ks0_chi2->clear(); B_J_chi2->clear(); B_chi2->clear(); B_chi2dof->clear();
   B_Prob->clear(); B_J_Prob->clear(); B_ks0_Prob->clear();


   bDecayVtxX->clear(); bDecayVtxY->clear(); bDecayVtxZ->clear(); 
   bDecayVtxXE->clear(); bDecayVtxYE->clear(); bDecayVtxZE->clear(); 
   bDecayVtxXYE->clear(); bDecayVtxXZE->clear(); bDecayVtxYZE->clear();  

   VDecayVtxX->clear(); VDecayVtxY->clear(); VDecayVtxZ->clear();
   VDecayVtxXE->clear(); VDecayVtxYE->clear(); VDecayVtxZE->clear();
   VDecayVtxXYE->clear(); VDecayVtxXZE->clear(); VDecayVtxYZE->clear();
   
//added for Ks01
  B_Ks01_mass->clear(); B_Ks01_px->clear(); B_Ks01_py->clear(); B_Ks01_pz->clear();
   B_Ks01_pt1->clear(); B_Ks01_px1->clear(); B_Ks01_py1->clear(); B_Ks01_pz1->clear(); B_Ks01_charge1->clear();
   B_Ks01_pt2->clear(); B_Ks01_px2->clear(); B_Ks01_py2->clear(); B_Ks01_pz2->clear(); B_Ks01_charge2->clear();

   B_Ks01_px1_track->clear(); B_Ks01_py1_track->clear(); B_Ks01_pz1_track->clear();
   B_Ks01_px2_track->clear(); B_Ks01_py2_track->clear(); B_Ks01_pz2_track->clear();
  B_Ks01_chi2->clear();
  VDecay1VtxX->clear(); VDecay1VtxY->clear(); VDecay1VtxZ->clear();
   VDecay1VtxXE->clear(); VDecay1VtxYE->clear(); VDecay1VtxZE->clear();
   VDecay1VtxXYE->clear(); VDecay1VtxXZE->clear(); VDecay1VtxYZE->clear();

   // *********

   nVtx = 0;
   priVtxX = 0; priVtxY = 0; priVtxZ = 0; 
   priVtxXE = 0; priVtxYE = 0; priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;
     
   pi1dxy->clear(); pi2dxy->clear(); pi1dz->clear(); pi2dz->clear();
   pi1dxy_e->clear(); pi2dxy_e->clear(); pi1dz_e->clear(); pi2dz_e->clear();

//added for Ks01
   pi11dxy->clear(); pi21dxy->clear(); pi11dz->clear(); pi21dz->clear();
   pi11dxy_e->clear(); pi21dxy_e->clear(); pi11dz_e->clear(); pi21dz_e->clear();

   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear();
   dr0->clear(); dr1->clear(); dpt0->clear(); dpt1->clear();
   dr01->clear(); dr11->clear(); dpt01->clear(); dpt11->clear();//added for Ks01
   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 

   // priRfVtxX->clear(); priRfVtxY->clear(); priRfVtxZ->clear(); priRfVtxXE->clear(); priRfVtxYE->clear(); 
   // priRfVtxZE->clear(); priRfVtxXYE->clear(); priRfVtxXZE->clear(); priRfVtxYZE->clear(); priRfVtxCL->clear(); 
   // priRfNTrkDif->clear(); 
   
   pVtxIPX->clear();  pVtxIPY->clear();  pVtxIPZ->clear();
   pVtxIPXE->clear();  pVtxIPYE->clear();  pVtxIPZE->clear();  pVtxIPCL->clear();
   pVtxIPXYE->clear();  pVtxIPXZE->clear();  pVtxIPYZE->clear();

   B_l3d->clear();  B_l3dE->clear();  B_lxy->clear(); B_lxyE->clear();
   B_cosalpha->clear();   B_cosalphaxy->clear(); alpha->clear();  B_treco->clear();   B_trecoe->clear();  B_trecoxy->clear(); B_trecoxye->clear();
   B_pvip->clear(); B_pviperr->clear(); B_pvips->clear(); B_pvlzip->clear(); B_pvlziperr->clear(); B_pvlzips->clear();
   B_pv2ip->clear(); B_pv2iperr->clear(); B_pv2ips->clear(); B_pv2lzip->clear(); B_pv2lziperr->clear(); B_pv2lzips->clear();
   
   B_l3d_pv2->clear();  B_l3dE_pv2->clear();
   B_iso->clear(); B_mum_iso->clear(); B_mup_iso->clear(); B_pi1_iso->clear();B_pi2_iso->clear();B_pi3_iso->clear();B_pi4_iso->clear();

   
   istruemum->clear(); istruemup->clear(); istruekp->clear(); istruekm->clear(); istruekp1->clear(); istruekm1->clear(); istruebs->clear();
   bunchXingMC->clear(); numInteractionsMC->clear(); trueNumInteractionsMC->clear();

   

}
bool JPsiKs0Ks0::IsTheSametk(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if(deltaR(tk.eta(),tk.phi(),mu.eta(),mu.phi()) < 0.00001 )return true;
  //if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool JPsiKs0Ks0::IsTheSame(const reco::Track& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool JPsiKs0Ks0::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
  if (ancestor == particle ) return true;
  for (size_t i=0; i< particle->numberOfMothers(); i++) {
    if (isAncestor(ancestor,particle->mother(i))) return true;
  }
  return false;
}

double JPsiKs0Ks0::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
  TVector3 pv_dv = decay_vtx - production_vtx;
  TVector3 b_p3  = b_p4.Vect();
  pv_dv.SetZ(0.);
  b_p3.SetZ(0.);
  Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
  return lxy*b_p4.M()/b_p3.Mag();
}

double
JPsiKs0Ks0::calEta (double Px, double Py, double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double
JPsiKs0Ks0::calPhi (double Px, double Py, double Pz)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

double
JPsiKs0Ks0::calEtaPhiDistance (double Px1, double Py1, double Pz1,
				double Px2, double Py2, double Pz2)
{
  double phi1 = calPhi (Px1,Py1,Pz1);
  double eta1 = calEta (Px1,Py1,Pz1);
  double phi2 = calPhi (Px2,Py2,Pz2);
  double eta2 = calEta (Px2,Py2,Pz2);
  return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}
void
JPsiKs0Ks0::saveTruthMatch(const edm::Event& iEvent){
  double deltaEtaPhi;

  for (vector<int>::size_type i = 0; i < B_mass->size(); i++) {//{{{
   
    //-----------------------
    // truth match with mu-
    //-----------------------
    double TruthMatchMuonMaxR_ = 0.004;
    double TruthMatchKaonMaxR_ = 0.3;
    deltaEtaPhi = calEtaPhiDistance(gen_muon1_p4.Px(), gen_muon1_p4.Py(), gen_muon1_p4.Pz(),
     				    B_J_px1->at(i), B_J_py1->at(i), B_J_pz1->at(i));
    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemum->push_back(true);
    } else {
      istruemum->push_back(false);
    }

    //-----------------------
    // truth match with mu+
    //-----------------------

    deltaEtaPhi = calEtaPhiDistance(gen_muon2_p4.Px(), gen_muon2_p4.Py(), gen_muon2_p4.Pz(),
     				    B_J_px2->at(i), B_J_py2->at(i), B_J_pz2->at(i));

    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemup->push_back(true);
    }
    else {
      istruemup->push_back(false);
    }

    //---------------------------------
    // truth match with pion+ track   
    //---------------------------------                          
    deltaEtaPhi = calEtaPhiDistance(gen_pion1_p4.Px(), gen_pion1_p4.Py(), gen_pion1_p4.Pz(),
                                    B_Ks0_px1->at(i), B_Ks0_py1->at(i), B_Ks0_pz1->at(i));
    if (deltaEtaPhi < TruthMatchKaonMaxR_){
      istruekp->push_back(true);
    } else {
      istruekp->push_back(false);
    }

    //---------------------------------                                                                                                                           
    // truth match with pion- track                                                                                                                                 
    //---------------------------------                                                                                                                           
    deltaEtaPhi = calEtaPhiDistance(gen_pion2_p4.Px(), gen_pion2_p4.Py(), gen_pion2_p4.Pz(),
                                    B_Ks0_px2->at(i), B_Ks0_py2->at(i), B_Ks0_pz2->at(i));
    if (deltaEtaPhi < TruthMatchKaonMaxR_){
      istruekm->push_back(true);
    } else {
      istruekm->push_back(false);
    }

    //---------------------------------
    //    // truth match with pion+ track   Ks01
    //        //---------------------------------   

deltaEtaPhi = calEtaPhiDistance(gen_pion3_p4.Px(), gen_pion3_p4.Py(), gen_pion3_p4.Pz(),
                                    B_Ks01_px1->at(i), B_Ks01_py1->at(i), B_Ks01_pz1->at(i));
    if (deltaEtaPhi < TruthMatchKaonMaxR_){
      istruekp1->push_back(true);
    } else {
      istruekp1->push_back(false);
    }


	//---------------------------------                                                                                                                           
    // truth match with pion- track                                                                                                     Ks01                            
        //---------------------------------       
    deltaEtaPhi = calEtaPhiDistance(gen_pion4_p4.Px(), gen_pion4_p4.Py(), gen_pion4_p4.Pz(),
                                    B_Ks01_px2->at(i), B_Ks01_py2->at(i), B_Ks01_pz2->at(i));
    if (deltaEtaPhi < TruthMatchKaonMaxR_){
      istruekm1->push_back(true);
    } else {
      istruekm1->push_back(false);
    }

    //---------------------------------------
    // truth match with Bs or Bs bar 
    //---------------------------------------                                                                                                
    if ( istruemum->back() && istruemup->back() && istruekm->back() && istruekp->back() && istruekm1->back() && istruekp1->back()) {
      istruebs->push_back(true);
    } else {
      istruebs->push_back(false);
    }



    }//}}}

}
void JPsiKs0Ks0::savePUinMC(const edm::Event& iEvent){
  // #################################
  // # Save pileup information in MC #
  // #################################
  // edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
  // iEvent.getByToken(puToken_, PupInfo);
  // //edm::LumiReWeighting lumi_weights;
  // //lumi_weights      = edm::LumiReWeighting("PileupMC_2018.root", "DataPileupHistogram2018_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  // //lumi_weights      = edm::LumiReWeighting("PileupMC_2016.root", "DataPileupHistogram2016_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  // //lumi_weights      = edm::LumiReWeighting("PileupMC_2017.root", "DataPileupHistogram2017_rereco.root", "input_Event/N_TrueInteractions", "pileup");
  // //float tnpv = -1;       // True number of primary vertices
  // //float wpu = 1;         // Pile-up re-weight factor
  
  // for (std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); PVI++)
  //   {
  //     bunchXingMC->push_back(PVI->getBunchCrossing());
  //     numInteractionsMC->push_back(PVI->getPU_NumInteractions());
  //     trueNumInteractionsMC->push_back(PVI->getTrueNumInteractions());
  //     //int bx = PVI->getBunchCrossing();
  //     //  if (bx == 0)tnpv = PVI->getTrueNumInteractions();        
  //   }

  //wpu = lumi_weights.weight(tnpv);
  //  fpuw8 = wpu ;
}


cov33_t JPsiKs0Ks0::GlobalError2SMatrix_33(GlobalError m_in) 
{
  cov33_t m_out;
  for (int i=0; i<3; i++) {
    for (int j=i; j<3; j++)  {
      m_out(i,j) = m_in.matrix()(i,j);
    }
  }
  return m_out;
}
  
cov99_t JPsiKs0Ks0::makeCovarianceMatrix(const cov33_t cov_vtx1,
			     const cov77_t cov_vtx2) 
{
  cov99_t cov;
  cov.Place_at(cov_vtx1,0,0);
  cov.Place_at(cov_vtx2.Sub<cov66_t>(0,0),3,3);
  return cov;
}
jac9_t JPsiKs0Ks0::makeJacobianVector3d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
 				     ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
				     const GlobalPoint &vtx2, const TVector3 &tv3momentum) 
{
  return makeJacobianVector3d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}

jac9_t JPsiKs0Ks0::makeJacobianVector3d(const AlgebraicVector3 &vtx1, 
				     const AlgebraicVector3 &vtx2, 
				     const AlgebraicVector3 &momentum) 
{
  jac9_t jac;
  const AlgebraicVector3 dist = vtx2 - vtx1;
  const double factor2 = 1. / ROOT::Math::Mag2(momentum);
  const double lifetime = ROOT::Math::Dot(dist, momentum) * factor2;
  jac.Place_at(-momentum*factor2,0);
  jac.Place_at( momentum*factor2,3);
  jac.Place_at( factor2*(dist-2*lifetime*momentum*factor2),6);
  return jac;
}
jac9_t JPsiKs0Ks0::makeJacobianVector2d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2,
				     const AlgebraicVector3 &momentum) {
  jac9_t jac;
  const double momentumMag = ROOT::Math::Mag(momentum);
  const AlgebraicVector3 dist = vtx2 - vtx1;
  const double distMag = ROOT::Math::Mag(dist);
  const double factorPositionComponent = 1./(distMag*momentumMag);
  const double factorMomentumComponent = 1./pow(momentumMag,3);
  jac(0)=-dist(0)*factorPositionComponent;
  jac(1)=-dist(1)*factorPositionComponent;
  jac(3)= dist(0)*factorPositionComponent;
  jac(4)= dist(1)*factorPositionComponent;
  jac(6)= momentum(0)*factorMomentumComponent;
  jac(7)= momentum(1)*factorMomentumComponent;
  return jac;
}

jac9_t JPsiKs0Ks0::makeJacobianVector2d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
			    ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
			    const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
  return makeJacobianVector2d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}


void JPsiKs0Ks0::saveIP(const RefCountedKinematicTree& vertexFitTree,
		     const reco::VertexCollection& vertices, int & pvIndex){

  vertexFitTree->movePointerToTheTop();

  RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
  RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
  
  std::cout<<" inside the save IP vertices size "<<vertices.size()<<std::endl;
  auto candTransientTrack = bCandMC->refittedTransientTrack();
  // find the first primary vertex
  const reco::Vertex* bestVertex_t(0);
  int bestVertexIndex(-1);
  double minDistance(999.);
  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);
    auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
    //std::cout<<" inside the vertex loop "<<i<<"\t"<<impactParameter3D.second.value()<<endl;
    if (impactParameter3D.first and impactParameter3D.second.value() < minDistance){
      minDistance = impactParameter3D.second.value();
      bestVertex_t = &vertex;
      bestVertexIndex = i;
      //std::cout<<" index list "<<i<< "\t"<<minDistance<<std::endl;
    }
  }
  pvIndex = bestVertexIndex;
  //std::cout<<" inside the save IP index of vtx "<<bestVertexIndex<<std::endl;
  pVtxIPX->push_back( bestVertex_t->x());
  pVtxIPY->push_back( bestVertex_t->y());
  pVtxIPZ->push_back( bestVertex_t->z());
  pVtxIPXE->push_back( bestVertex_t->covariance(0, 0) );
  pVtxIPYE->push_back( bestVertex_t->covariance(1, 1) );
  pVtxIPZE->push_back( bestVertex_t->covariance(2, 2) );
  pVtxIPXYE->push_back(bestVertex_t->covariance(0, 1) );
  pVtxIPXZE->push_back(bestVertex_t->covariance(0, 2) );  
  pVtxIPYZE->push_back( bestVertex_t->covariance(1, 2) );
  pVtxIPCL->push_back(  ChiSquaredProbability((double)(bestVertex_t->chi2()),(double)(bestVertex_t->ndof())) );      

  // find second best vertex
  const reco::Vertex* bestVertex2(0);
  //int bestVertexIndex2(-1);
  double minDistance2(999.);
  for ( unsigned int i = 0; i<vertices.size(); ++i ){
    const auto & vertex = vertices.at(i);    
    auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
    if (impactParameter3D.first and impactParameter3D.second.value() < minDistance2 and impactParameter3D.second.value() > minDistance){
      minDistance2 = impactParameter3D.second.value();
      bestVertex2 = &vertex;
      //bestVertexIndex2 = i;
    }
  }
  //std::cout<<" inside the save IP second vertex "<<std::endl;
  //  if (! bestVertex) continue;
  auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex_t);
  auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex_t);
  //pv = bestVertex;
  //pvIndex = bestVertexIndex;
  double longitudinalImpactParameter(0.0), longitudinalImpactParameterErr(0.0);
  double distaceOfClosestApproach(0.0), distaceOfClosestApproachErr(0.0) ;
  if (impactParameterZ.first) {
    longitudinalImpactParameter    = impactParameterZ.second.value();
    longitudinalImpactParameterErr = impactParameterZ.second.error();
  }
  if(impactParameter3D.first) {
    distaceOfClosestApproach       = impactParameter3D.second.value();
    distaceOfClosestApproachErr    = impactParameter3D.second.error();
  }
  //  std::cout<<" inside dca "<<distaceOfClosestApproach<<"\t "<<distaceOfClosestApproachErr<<"\t"<<distaceOfClosestApproachErr/distaceOfClosestApproach<<std::endl;
  //std::cout<<" inside long "<<longitudinalImpactParameter<<"\t"<<longitudinalImpactParameterErr<<"\t"<<longitudinalImpactParameterErr/longitudinalImpactParameter<<std::endl;

  B_pvip->push_back(distaceOfClosestApproach);
  B_pviperr->push_back(distaceOfClosestApproachErr);
  B_pvips->push_back(distaceOfClosestApproachErr/distaceOfClosestApproach);
  B_pvlzip->push_back(longitudinalImpactParameter);
  B_pvlziperr->push_back(longitudinalImpactParameterErr);
  B_pvlzips->push_back(longitudinalImpactParameterErr/longitudinalImpactParameter);
  
  //std::cout<<" inside pvip "<<longitudinalImpactParameter<<"\t"<<longitudinalImpactParameterErr<<"\t"<<longitudinalImpactParameterErr/longitudinalImpactParameter<<std::endl;
  
  VertexDistance3D distance3D;
  auto dist = distance3D.distance(*bestVertex_t, bDecayVertexMC->vertexState() );
  double decayLength(-1.), decayLengthErr(0);
  decayLength    = dist.value();
  decayLengthErr = dist.error();
  
  VertexDistanceXY distanceXY;
  auto distXY = distanceXY.distance(*bestVertex_t, bDecayVertexMC->vertexState() );

  B_l3d ->push_back(decayLength);
  B_l3dE ->push_back(decayLengthErr);
  B_lxy ->push_back(distXY.value());
  B_lxyE ->push_back(distXY.error());
  
  if (bestVertex2){
    double longitudinalImpactParameter2(0.0), longitudinalImpactParameter2Err(0.0);
    double distaceOfClosestApproach2(0.0), distaceOfClosestApproach2Err(0.0) ;
    auto impactParameter3D2 = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex2);
    auto impactParameterZ2  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex2);
    //pv2 = bestVertex2;
    //pv2Index = bestVertexIndex2;
    if (impactParameterZ2.first) {
      longitudinalImpactParameter2    = impactParameterZ2.second.value();
      longitudinalImpactParameter2Err = impactParameterZ2.second.error();
    }
    if (impactParameter3D2.first) {
      distaceOfClosestApproach2       = impactParameter3D2.second.value();
      distaceOfClosestApproach2Err    = impactParameter3D2.second.error();
    }
    B_pv2ip->push_back(distaceOfClosestApproach2);
    B_pv2iperr->push_back(distaceOfClosestApproach2Err);
    B_pv2ips->push_back(distaceOfClosestApproach2Err/distaceOfClosestApproach2);
    B_pv2lzip->push_back(longitudinalImpactParameter2);
    B_pv2lziperr->push_back(longitudinalImpactParameter2Err);
    B_pv2lzips->push_back(longitudinalImpactParameter2Err/longitudinalImpactParameter2);
    
    // compute decay length
    VertexDistance3D distance3D;
    auto dist = distance3D.distance(*bestVertex2, bDecayVertexMC->vertexState() );
    B_l3d_pv2 ->push_back(dist.value());
    B_l3dE_pv2 ->push_back(dist.error());
  }
  TVector3 plab(bCandMC->currentState().globalMomentum().x(),
		bCandMC->currentState().globalMomentum().y(),
		bCandMC->currentState().globalMomentum().z());
  TVector3 p1(bestVertex_t->x(), bestVertex_t->y(), bestVertex_t->z());
  TVector3 p2(bDecayVertexMC->vertexState().position().x(), 
	      bDecayVertexMC->vertexState().position().y(), 
	      bDecayVertexMC->vertexState().position().z());
  TVector3 pDiff = p2-p1;
  TVector3 pDiffXY = TVector3(pDiff.X(), pDiff.Y(), 0.);
  TVector3 ptrans  = TVector3(plab.X(), plab.Y(), 0.);
  double cosAlpha(-999.),cosAlphaXY(-999.),decayTime(-999.),decayTimeError(-999.), decayTimeXY(-999.),decayTimeXYError(-999.);
  cosAlpha  = plab.Dot(pDiff) / (plab.Mag() * pDiff.Mag());
  cosAlphaXY  = ptrans.Dot(pDiffXY) / (ptrans.Mag() * pDiffXY.Mag());
  
  B_cosalpha -> push_back(cosAlpha);
  B_cosalphaxy -> push_back(cosAlphaXY);
  alpha->push_back(TMath::ACos(cosAlpha));
  // compute decayTime information

  const double massOverC =  bCandMC->currentState().mass()/TMath::Ccgs();

  // get covariance matrix for error propagation in decayTime calculation
  auto vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(bestVertex_t->error()),
					     bCandMC->currentState().kinematicParametersError().matrix());
  auto vtxDistanceJac3d = makeJacobianVector3d(bestVertex_t->position(), bDecayVertexMC->vertexState().position(), plab);
  auto vtxDistanceJac2d = makeJacobianVector2d(bestVertex_t->position(), bDecayVertexMC->vertexState().position(), plab);

  decayTime = dist.value() / plab.Mag() * cosAlpha * massOverC;
  decayTimeError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac3d)) * massOverC;

  decayTimeXY = distXY.value() / plab.Perp() * cosAlphaXY * massOverC;
  decayTimeXYError = TMath::Sqrt(ROOT::Math::Similarity(vtxDistanceCov, vtxDistanceJac2d)) * massOverC;

  B_treco -> push_back(decayTime);
  B_trecoe -> push_back(decayTimeError);
  B_trecoxy -> push_back(decayTimeXY);
  B_trecoxye -> push_back(decayTimeXYError);
  
}
// Isolation
void JPsiKs0Ks0::SaveIso(const RefCountedKinematicTree& vertexFitTree, const MagneticField *fMagneticField, 
		      edm::ESHandle<TransientTrackBuilder> theB,
		      unsigned int pvIndex,
		      edm::Handle< edm::View<pat::PackedCandidate> > thePATTrackHandle,
		      //edm::Handle<std::vector<pat::GenericParticle> > thePATTrackHandle,
		      edm::Handle< edm::View<pat::Muon> > thePATMuonHandle,
		      const reco::TransientTrack muon1TT, 
		      const reco::TransientTrack muon2TT,
		      const reco::TransientTrack pion1TT,
		      const reco::TransientTrack pion2TT,		     
		     const reco::TransientTrack pion3TT,
                      const reco::TransientTrack pion4TT
 

		      ){

  vertexFitTree->movePointerToTheTop(); // Bs --> Jpsi Kss
  RefCountedKinematicParticle b_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex vertex = vertexFitTree->currentDecayVertex();
  
  std::cout<<"pv index in iso "<<pvIndex<<std::endl;
  
  ClosestApproachInRPhi ClosestApp;
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  
  double  bpx = b_KP->currentState().globalMomentum().x();
  double  bpy = b_KP->currentState().globalMomentum().y();
  double  bpz = b_KP->currentState().globalMomentum().z();
  double masstmp = b_KP->currentState().mass();
  TLorentzVector tmp_bmeson_lv;
  tmp_bmeson_lv.SetXYZM(bpx,bpy, bpz, masstmp);
  
  double sumppt(0),summpt(0),sumtrkmpt(0),sumtrkppt(0);
  double sumBpt(0);
  //std::cout<<" Inside the Isolation folder "<<std::endl;
  //std::vector<std::pair<int,std::pair<float,float> > > fNstTracks;
  for(edm::View<pat::PackedCandidate>::const_iterator iTrack = thePATTrackHandle->begin();                                                                 
      iTrack != thePATTrackHandle->end(); ++iTrack )   {
    
    if(iTrack->charge()==0) continue;                                                                                                           
    if(!iTrack->trackHighPurity()) continue;                                                                                                  
    if(!iTrack->hasTrackDetails()) continue;
    bool skip_this_track = false; 
    for(edm::View<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin(); iMuon != thePATMuonHandle->end(); ++iMuon)
      {
	if ( IsTheSametk(*iTrack,*iMuon) ) skip_this_track = true;
	
      }	
    
    
    if (iTrack->vertexRef().key()!=pvIndex) continue;//{ std::cout<<"passing our pvindex "<<std::endl;}// continue;
    
    if(skip_this_track)continue;    
    reco::TransientTrack TrackIsoTT((*theB).build(iTrack->pseudoTrack()));
    //std::cout<<" Inside the Isolation and checking track details "<<iTrack->pt()<<std::endl;

    if( !(TrackIsoTT.isValid()) )continue;    
    if(iTrack->pt()<0.8) continue;
    bool skipping_track= false;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), muon1TT.track().eta(),  muon1TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==muon1TT.track().charge())) skipping_track = true;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), muon2TT.track().eta(),  muon2TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==muon2TT.track().charge())) skipping_track = true;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), pion1TT.track().eta(),  pion1TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==pion1TT.track().charge())) skipping_track = true;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), pion2TT.track().eta(),  pion2TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==pion2TT.track().charge())) skipping_track = true;
 if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), pion3TT.track().eta(),  pion3TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==pion3TT.track().charge())) skipping_track = true;
    if(deltaR(TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), pion4TT.track().eta(),  pion4TT.track().phi())<0.0001 && (TrackIsoTT.track().charge()==pion4TT.track().charge())) skipping_track = true;
    if(skipping_track) continue;
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), muon1TT.initialFreeState());
    if (ClosestApp.status() != false)
      {
	if ( ClosestApp.distance() < 0.1 ) {
	  double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), muon1TT.track().eta(),  muon1TT.track().phi());
	  if(deltaR_tmp<0.5) sumppt += iTrack->pt();
	}                 
      }
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), muon2TT.initialFreeState());
    if (ClosestApp.status() != false)
      {
	if ( ClosestApp.distance() < 0.1 ) {
	  double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), muon2TT.track().eta(),  muon2TT.track().phi());
	  if(deltaR_tmp<0.5) summpt += iTrack->pt();
	}
      }
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), pion1TT.initialFreeState());
    if (ClosestApp.status() != false)
      {
	if ( ClosestApp.distance() < 0.1 ) {
	  double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), pion1TT.track().eta(),  pion1TT.track().phi());
	    if(deltaR_tmp<0.5) sumtrkmpt += iTrack->pt();
	}                 
	}
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), pion2TT.initialFreeState());
      if (ClosestApp.status() != false)
	{
	  if ( ClosestApp.distance() < 0.1 ) {
	    double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), pion2TT.track().eta(),  pion2TT.track().phi());
	    if(deltaR_tmp<0.5) sumtrkppt += iTrack->pt();
	  }
	}

	ClosestApp.calculate(TrackIsoTT.initialFreeState(), pion3TT.initialFreeState());
    if (ClosestApp.status() != false)
      { 
        if ( ClosestApp.distance() < 0.1 ) {
          double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), pion3TT.track().eta(),  pion3TT.track().phi());
            if(deltaR_tmp<0.5) sumtrkmpt += iTrack->pt();
        }
        }
    ClosestApp.calculate(TrackIsoTT.initialFreeState(), pion4TT.initialFreeState());
      if (ClosestApp.status() != false)
        {
          if ( ClosestApp.distance() < 0.1 ) {
            double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), pion4TT.track().eta(),  pion4TT.track().phi());
            if(deltaR_tmp<0.5) sumtrkppt += iTrack->pt();
          }
        }
  
      
      VertexDistance3D distance3D;      
      const GlobalPoint BVP = GlobalPoint( vertex->position() );
      if(vertex->vertexIsValid()){       
	TrajectoryStateOnSurface  tsos = extrapolator.extrapolate(TrackIsoTT.initialFreeState(), BVP);
	if (tsos.isValid()) {      			   
	  Measurement1D doca = distance3D.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), vertex->vertexState());
	  double svDoca = doca.value();	    
	  //docatrks.push_back(svDoca);
	  // check to ensure the goodness of the track
	  if( svDoca<0.05){
	    double deltaR_tmp = deltaR(iTrack->eta(), iTrack->phi(), tmp_bmeson_lv.Eta(), tmp_bmeson_lv.Phi());			   
	    if(deltaR_tmp<0.7) sumBpt += iTrack->pt();			     
	  }
	  
	}
      }
    }

  float BIso_ =tmp_bmeson_lv.Pt()/(sumBpt +tmp_bmeson_lv.Pt());
  float mupIso_ = -99., mumIso_ = -99.0;
  if(muon1TT.track().charge() == 1){
    mupIso_ = muon1TT.track().pt()/(sumppt +muon1TT.track().pt() );
    mumIso_ = muon2TT.track().pt()/(summpt +muon2TT.track().pt() );
  }else {
    mupIso_ = muon2TT.track().pt()/(sumppt +muon2TT.track().pt() );
    mumIso_ = muon1TT.track().pt()/(summpt +muon1TT.track().pt() );
  }

  float trkpIso_ = pion1TT.track().pt()/(sumtrkppt +pion1TT.track().pt() );
  float trkmIso_ = pion2TT.track().pt()/(sumtrkmpt +pion2TT.track().pt() );
  float trk1pIso_ = pion3TT.track().pt()/(sumtrkppt +pion3TT.track().pt() );
  float trk1mIso_ = pion4TT.track().pt()/(sumtrkmpt +pion4TT.track().pt() );

  B_iso->push_back(BIso_); 
  B_mum_iso->push_back(mumIso_); 
  B_mup_iso->push_back(mupIso_); 
  B_pi1_iso->push_back(trkpIso_);
  B_pi2_iso->push_back(trkmIso_);
  B_pi3_iso->push_back(trk1pIso_);
  B_pi4_iso->push_back(trk1mIso_);

  
}
// ------------ method called once each job just before starting event loop  ------------

void 
JPsiKs0Ks0::beginJob()
{

  std::cout << "Beginning analyzer job with value of doMC = " << doMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bs->J/psi Ks0 ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_Ks0_mass", &B_Ks0_mass);
  tree_->Branch("B_Ks0_px", &B_Ks0_px);
  tree_->Branch("B_Ks0_py", &B_Ks0_py);
  tree_->Branch("B_Ks0_pz", &B_Ks0_pz);
 
  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_Ks0_pt1", &B_Ks0_pt1);
  tree_->Branch("B_Ks0_px1", &B_Ks0_px1);
  tree_->Branch("B_Ks0_py1", &B_Ks0_py1);
  tree_->Branch("B_Ks0_pz1", &B_Ks0_pz1);
  tree_->Branch("B_Ks0_px1_track", &B_Ks0_px1_track);
  tree_->Branch("B_Ks0_py1_track", &B_Ks0_py1_track);
  tree_->Branch("B_Ks0_pz1_track", &B_Ks0_pz1_track);
  tree_->Branch("B_Ks0_charge1", &B_Ks0_charge1); 
 
  tree_->Branch("B_Ks0_pt2", &B_Ks0_pt2);
  tree_->Branch("B_Ks0_px2", &B_Ks0_px2);
  tree_->Branch("B_Ks0_py2", &B_Ks0_py2);
  tree_->Branch("B_Ks0_pz2", &B_Ks0_pz2);
  tree_->Branch("B_Ks0_px2_track", &B_Ks0_px2_track);
  tree_->Branch("B_Ks0_py2_track", &B_Ks0_py2_track);
  tree_->Branch("B_Ks0_pz2_track", &B_Ks0_pz2_track);
  tree_->Branch("B_Ks0_charge2", &B_Ks0_charge2);

  tree_->Branch("B_Ks01_pt1", &B_Ks01_pt1);
  tree_->Branch("B_Ks01_px1", &B_Ks01_px1);
  tree_->Branch("B_Ks01_py1", &B_Ks01_py1);
  tree_->Branch("B_Ks01_pz1", &B_Ks01_pz1);
  tree_->Branch("B_Ks01_px1_track", &B_Ks01_px1_track);
  tree_->Branch("B_Ks01_py1_track", &B_Ks01_py1_track);
  tree_->Branch("B_Ks01_pz1_track", &B_Ks01_pz1_track);
  tree_->Branch("B_Ks01_charge1", &B_Ks01_charge1);

  tree_->Branch("B_Ks01_pt2", &B_Ks01_pt2);
  tree_->Branch("B_Ks01_px2", &B_Ks01_px2);
  tree_->Branch("B_Ks01_py2", &B_Ks01_py2);
  tree_->Branch("B_Ks01_pz2", &B_Ks01_pz2);
  tree_->Branch("B_Ks01_px2_track", &B_Ks01_px2_track);
  tree_->Branch("B_Ks01_py2_track", &B_Ks01_py2_track);
  tree_->Branch("B_Ks01_pz2_track", &B_Ks01_pz2_track);
  tree_->Branch("B_Ks01_charge2", &B_Ks01_charge2);


  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("B_chi2", &B_chi2);
  tree_->Branch("B_chi2dof", &B_chi2dof);
  tree_->Branch("B_Ks0_chi2", &B_Ks0_chi2);
  tree_->Branch("B_Ks01_chi2", &B_Ks01_chi2);
  tree_->Branch("B_J_chi2", &B_J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_ks0_Prob", &B_ks0_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
       
  // *************************

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_->Branch("bDecayVtxX",&bDecayVtxX);
  tree_->Branch("bDecayVtxY",&bDecayVtxY);
  tree_->Branch("bDecayVtxZ",&bDecayVtxZ);
  tree_->Branch("bDecayVtxXE",&bDecayVtxXE);
  tree_->Branch("bDecayVtxYE",&bDecayVtxYE);
  tree_->Branch("bDecayVtxZE",&bDecayVtxZE);
  tree_->Branch("bDecayVtxXYE",&bDecayVtxXYE);
  tree_->Branch("bDecayVtxXZE",&bDecayVtxXZE);
  tree_->Branch("bDecayVtxYZE",&bDecayVtxYZE);

  tree_->Branch("VDecayVtxX",&VDecayVtxX);
  tree_->Branch("VDecayVtxY",&VDecayVtxY);
  tree_->Branch("VDecayVtxZ",&VDecayVtxZ);
  tree_->Branch("VDecayVtxXE",&VDecayVtxXE);
  tree_->Branch("VDecayVtxYE",&VDecayVtxYE);
  tree_->Branch("VDecayVtxZE",&VDecayVtxZE);
  tree_->Branch("VDecayVtxXYE",&VDecayVtxXYE);
  tree_->Branch("VDecayVtxXZE",&VDecayVtxXZE);
  tree_->Branch("VDecayVtxYZE",&VDecayVtxYZE);

//added for ks01
  tree_->Branch("VDecay1VtxX",&VDecayVtxX);
  tree_->Branch("VDecay1VtxY",&VDecayVtxY);
  tree_->Branch("VDecay1VtxZ",&VDecayVtxZ);
  tree_->Branch("VDecay1VtxXE",&VDecayVtxXE);
  tree_->Branch("VDecay1VtxYE",&VDecayVtxYE);
  tree_->Branch("VDecay1VtxZE",&VDecayVtxZE);
  tree_->Branch("VDecay1VtxXYE",&VDecayVtxXYE);
  tree_->Branch("VDecay1VtxXZE",&VDecayVtxXZE);
  tree_->Branch("VDecay1VtxYZE",&VDecayVtxYZE);

  tree_->Branch("pVtxIPX",     &pVtxIPX);
  tree_->Branch("pVtxIPY",     &pVtxIPY);
  tree_->Branch("pVtxIPZ",     &pVtxIPZ);
  tree_->Branch("pVtxIPXE",     &pVtxIPXE);
  tree_->Branch("pVtxIPYE",     &pVtxIPYE);
  tree_->Branch("pVtxIPZE",     &pVtxIPZE);
  tree_->Branch("pVtxIPXYE",     &pVtxIPXYE);
  tree_->Branch("pVtxIPXZE",     &pVtxIPXZE);
  tree_->Branch("pVtxIPYZE",     &pVtxIPYZE);
  tree_->Branch("pVtxIPCL",     &pVtxIPCL);

  // tree_->Branch("priRfVtxX",&priRfVtxX);
  // tree_->Branch("priRfVtxY",&priRfVtxY);
  // tree_->Branch("priRfVtxZ",&priRfVtxZ);
  // tree_->Branch("priRfVtxXE",&priRfVtxXE);
  // tree_->Branch("priRfVtxYE",&priRfVtxYE);
  // tree_->Branch("priRfVtxZE",&priRfVtxZE);
  // tree_->Branch("priRfVtxXYE",&priRfVtxXYE);
  // tree_->Branch("priRfVtxXZE",&priRfVtxXZE);
  // tree_->Branch("priRfVtxYZE",&priRfVtxYZE);
  // tree_->Branch("priRfVtxCL",&priRfVtxCL);
  // tree_->Branch("priRfNTrkDif",&priRfNTrkDif);

  tree_->Branch("pi1dxy",&pi1dxy);
  tree_->Branch("pi2dxy",&pi2dxy);
  tree_->Branch("pi1dz",&pi1dz);
  tree_->Branch("pi2dz",&pi2dz);

  tree_->Branch("pi1dxy_e",&pi1dxy_e);
  tree_->Branch("pi2dxy_e",&pi2dxy_e);
  tree_->Branch("pi1dz_e",&pi1dz_e);
  tree_->Branch("pi2dz_e",&pi2dz_e);
//added for Ks01
 tree_->Branch("pi11dxy",&pi11dxy);
  tree_->Branch("pi21dxy",&pi21dxy);
  tree_->Branch("pi11dz",&pi11dz);
  tree_->Branch("pi21dz",&pi21dz);

  tree_->Branch("pi11dxy_e",&pi11dxy_e);
  tree_->Branch("pi21dxy_e",&pi21dxy_e);
  tree_->Branch("pi11dz_e",&pi11dz_e);
  tree_->Branch("pi21dz_e",&pi21dz_e);

  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
  tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk); 
  tree_->Branch("dr0",&dr0);
  tree_->Branch("dr1",&dr1);
  tree_->Branch("dr01",&dr01);
  tree_->Branch("dr11",&dr11);
  tree_->Branch("dpt0",&dpt0);
  tree_->Branch("dpt1",&dpt1);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);



  tree_->Branch("B_pvip",&B_pvip);
  tree_->Branch("B_pviperr",&B_pviperr);
  tree_->Branch("B_pvips",&B_pvips);
  tree_->Branch("B_pvlzip",&B_pvlzip);
  tree_->Branch("B_pvlziperr",&B_pvlziperr);
  tree_->Branch("B_pvlzips",&B_pvlzips);
  tree_->Branch("B_pv2ip",&B_pv2ip);
  tree_->Branch("B_pv2iperr",&B_pv2iperr);
  tree_->Branch("B_pv2ips",&B_pv2ips);
  tree_->Branch("B_pv2lzip",&B_pv2lzip);
  tree_->Branch("B_pv2lziperr",&B_pv2lziperr);
  tree_->Branch("B_pv2lzips",&B_pv2lzips);
  tree_->Branch("B_l3d_pv2",&B_l3d_pv2);
  tree_->Branch("B_l3dE_pv2",&B_l3dE_pv2);
  tree_->Branch("B_l3d",&B_l3d);
  tree_->Branch("B_l3dE",&B_l3dE);
  tree_->Branch("B_lxy", &B_lxy);
  tree_->Branch("B_lxyE",&B_lxyE);
  tree_->Branch("B_cosalpha",&B_cosalpha);  
  tree_->Branch("B_cosalphaxy",&B_cosalphaxy);
  tree_->Branch("alpha",&alpha);
  tree_->Branch("B_treco",&B_treco);
  tree_->Branch("B_trecoe",&B_trecoe);
  tree_->Branch("B_trecoxy",&B_trecoxy);
  tree_->Branch("B_trecoxye",&B_trecoxye);
  tree_->Branch("B_iso",&B_iso);
  tree_->Branch("B_mum_iso",&B_mum_iso);
  tree_->Branch("B_mup_iso",&B_mup_iso);
  tree_->Branch("B_pi1_iso",&B_pi1_iso);
  tree_->Branch("B_pi2_iso",&B_pi2_iso);
  tree_->Branch("B_pi3_iso",&B_pi3_iso);
  tree_->Branch("B_pi4_iso",&B_pi4_iso);

  
  
  //gen information
  if (isMC_) {
    tree_->Branch("gen_b_p4",     "TLorentzVector",  &gen_b_p4);
    tree_->Branch("gen_ks_p4",   "TLorentzVector",  &gen_ks_p4);
    tree_->Branch("gen_pion1_p4",  "TLorentzVector",  &gen_pion1_p4);
    tree_->Branch("gen_pion2_p4",  "TLorentzVector",  &gen_pion2_p4);
    tree_->Branch("gen_pion3_p4",  "TLorentzVector",  &gen_pion3_p4);
    tree_->Branch("gen_pion4_p4",  "TLorentzVector",  &gen_pion4_p4);//added for ks01
    tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
    tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
    tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
    tree_->Branch("gen_b_vtx",    "TVector3",        &gen_b_vtx);
    tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
    tree_->Branch("gen_b_ct",     &gen_b_ct,        "gen_b_ct/F");
  }
  tree_->Branch("istruemum",  &istruemum );
  tree_->Branch("istruemup",  &istruemup );
  tree_->Branch("istruekp",   &istruekp  );
  tree_->Branch("istruekm",   &istruekm  );
  tree_->Branch("istruekp1",   &istruekp1  );
  tree_->Branch("istruekm1",   &istruekm1  );

  tree_->Branch("istruebs",   &istruebs  );
  
  tree_->Branch("bunchXingMC",&bunchXingMC);
  tree_->Branch("numInteractionsMC",&numInteractionsMC);
  tree_->Branch("trueNumInteractionsMC",&trueNumInteractionsMC);
    //}


}
// ------------ method called once each job just after ending the event loop  ------------
void JPsiKs0Ks0::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKs0Ks0);

