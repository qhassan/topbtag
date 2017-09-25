// -*- C++ -*-
//
// Package:    test_fi/Testfi
// Class:      Testfi
// 
/**\class Testfi Testfi.cc test_fi/Testfi/plugins/Testfi.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Qamar Ul Hassan
//         Created:  Fri, 22 Sep 2017 16:10:52 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Eff : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Eff(const edm::ParameterSet&);
      ~Eff();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

      std::unordered_map<std::string,TH1*> histContainer_;
      std::unordered_map<std::string,TH2*> histContainer2d_;

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
Eff::Eff(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))


{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


Eff::~Eff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Eff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//   using namespace edm;
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &primVtx = vertices->front();

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  float rho=*rhoH;
  rho=rho;


    histContainer_["cutflow"]->Fill(0);
    histContainer_["ecutflow"]->Fill(0);
    histContainer_["mucutflow"]->Fill(0);

//
  // MUONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId  
//      
//
  
//  float leptonpt(0), leptonphi(0);
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  histContainer_["nrecomuons"]->Fill(muons->size());
  std::vector<const pat::Muon *> selectedMuons,vetoMuons;

  for (const pat::Muon &mu : *muons) {

        //kinematics
        bool passPt( mu.pt() > 26 );
        bool passVetoPt( mu.pt()>10 );
        bool passEta(fabs(mu.eta()) < 2.1 );

        //distance to the PV
        float dz(fabs( mu.vertex().z() - primVtx.z()));
        bool passDB( mu.dB()<0.2 && dz<0.5 );

        //isolation
        float relchIso((mu.chargedHadronIso())/mu.pt());
        bool passIso( relchIso<0.2 );

        if(mu.isPFMuon()
           && mu.isGlobalMuon()
           && mu.normChi2() < 10
           && mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
           && mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
           && mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
           && mu.numberOfMatchedStations() > 1 )
              {
                if( passPt && passEta && passDB && passIso)
                  {
                    selectedMuons.push_back( &mu );
                  }
                else if(passVetoPt && passEta)
                   {
                     vetoMuons.push_back( &mu );
                   }
                  }

              }
  

histContainer_["nselmuons"]->Fill(selectedMuons.size());

  //
  // ELECTRONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification  
  // 
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  histContainer_["nrecoelectrons"]->Fill(electrons->size());
  std::vector<const pat::Electron *> selectedElectrons,vetoElectrons;
  for (const pat::Electron &el : *electrons) {

      //kinematics cuts
      bool passPt(el.pt()>30);
      bool passVetoPt(el.pt()>20);
      bool passEta(fabs(el.eta()) < 2.5 && (fabs(el.superCluster()->eta()) < 1.4442 || fabs(el.superCluster()->eta()) > 1.5660));

      //isolation
      float relchIso((el.chargedHadronIso())/el.pt());
      bool passIso( relchIso<0.2 );
      bool passVetoIso( relchIso<0.15 );
      if(//el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <=0
        el.dB() < 0.02
        && el.passConversionVeto() == true
        )
              {

      if(passPt && passEta && passIso)
        {
        selectedElectrons.push_back(&el);
        }

      else if(passVetoPt && passEta && passVetoIso)
        {
          vetoElectrons.push_back(&el);
        }

            }

  }
  histContainer_["nselelectrons"]->Fill(selectedElectrons.size());

      //require only 1 tight lepton in the event
      int nSelectedLeptons(selectedElectrons.size()+selectedMuons.size());
      if(nSelectedLeptons>1 || nSelectedLeptons==0) return;
      histContainer_["cutflow"]->Fill(1);
      if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(1);
      if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(1);

      //require no other leptons in the event
      int nVetoLeptons(vetoElectrons.size()+vetoMuons.size());
      if(nVetoLeptons>0) return;
      histContainer_["cutflow"]->Fill(2);
      if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(2);
      if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(2);

      //
      // JETS
      //

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  std::vector<const pat::Jet *> selectedJets;
  histContainer_["nrecojets"]->Fill(jets->size());
  int njets30(0), nj(0), nBtags(0);
  for (const pat::Jet &j : *jets) {

     float dR2lepton= selectedMuons.size()==1 ?
      deltaR(j,*(selectedMuons[0])) :
      deltaR(j,*(selectedElectrons[0]));
     float rawEnergy(j.energy()*j.jecFactor("Uncorrected"));
     if ( j.numberOfDaughters() > 1
        && (j.neutralHadronEnergy() + j.HFHadronEnergy())/rawEnergy < 0.99
        && j.neutralEmEnergyFraction() < 0.99
        && (j.chargedEmEnergyFraction() < 0.99 || fabs(j.eta()) >= 2.4)
        && (j.chargedHadronEnergyFraction() > 0. || fabs(j.eta()) >= 2.4)
        && (j.chargedMultiplicity() > 0 || fabs(j.eta()) >= 2.4)
        && dR2lepton>0.4) 
           {
             //parton matched to the jet
             const reco::Candidate *genParton = j.genParton();
             bool isLightFromW(false),isB(false);
             if(genParton)
               {
                isB=(abs(genParton->pdgId())==5);
                if(genParton->mother())
                isLightFromW=(abs(genParton->mother()->pdgId())==24);
               }
               //loop over jet charged constituents
               int ncharged(0);
               TLorentzVector chargedJet(0,0,0,0);
               float sumptcharged(0);
               for(size_t ipf=0; ipf<j.numberOfDaughters(); ipf++)
                  {
                   const pat::PackedCandidate *pfConst=dynamic_cast<const pat::PackedCandidate *>(j.daughter(ipf));
                   if(pfConst==0) continue;
                   if(pfConst->charge()==0 || pfConst->fromPV()==0) continue;
                   sumptcharged += pfConst->pt();
                   chargedJet += TLorentzVector(pfConst->px(),pfConst->py(),pfConst->pz(),pfConst->energy());
                   ncharged++;
                  }
              //N-1 plots
              if(fabs(j.eta()) < 2.5)
                {
                 histContainer_["jetpt"]->Fill(j.pt());
                 histContainer_["chjetpt"]->Fill(sumptcharged);
                 if(isB)
                   {
                    histContainer_["bjetpt"]->Fill(j.pt());
                    histContainer_["bchjetpt"]->Fill(sumptcharged);
                   }
                 else if(isLightFromW)
                   {
                    histContainer_["lightjetpt"]->Fill(j.pt());
                    histContainer_["lightchjetpt"]->Fill(sumptcharged);
                   }
                 else
                   {
                    histContainer_["otherjetpt"]->Fill(j.pt());
                    histContainer_["otherchjetpt"]->Fill(sumptcharged);
                   }

              }
  if(sumptcharged>15 /*j.pt() > 30*/)         histContainer_["jeteta"]->Fill(fabs(j.eta()));
  if( fabs(j.eta()) < 2.5 && sumptcharged>15 /*j.pt() > 30*/)
    {
      if(j.pt()>30) njets30++;
      selectedJets.push_back( &j );
//      float csv=j.bDiscriminator("combinedSecondaryVertexBJetTags");
      float csv=j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      if(csv>0.800) nBtags++;
      histContainer_["jetcsv"]->Fill(csv);
      histContainer_["jetpileupid"]->Fill(j.userFloat("pileupJetId:fullDiscriminant"));
      nj++;
    }
//  histContainer_["jetcsv"]->Fill(csv);
  histContainer_["nseljets"]->Fill(selectedJets.size());
  histContainer_["nseljetsfull"]->Fill(njets30);

           }
  }

   //
   // MET
   //
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
//  float metpt = mets->at(0).pt();
//  float metphi = mets->at(0).phi();
//  float dphi_met_lepton = deltaPhi(leptonphi, metphi); // use the function to restrict to the 0,pi range
//  float mt=sqrt(2*leptonpt*metpt*(1-cos(dphi_met_lepton)));

  //charged met
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  TLorentzVector chMet(0,0,0,0);
  float chHT(0);
  for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
  //require not to be associated to other PVs and to be charged
  const pat::PackedCandidate &pf = (*pfs)[i];
  if (pf.fromPV() == 0 || pf.charge()== 0) continue;
  chMet -= TLorentzVector(pf.px(),pf.py(),0,pf.pt());
  chHT += pf.pt();
  }
//  float dphi_chmet_lepton = deltaPhi(leptonphi, chMet.Phi()); // use the function to restrict to the 0,pi range
//  float chmt=sqrt(2*leptonpt*chMet.Pt()*(1-cos(dphi_chmet_lepton)));

  //
  //  FINAL SELECTION PLOTS
  // 
  if(selectedJets.size()>=2 && nBtags>=2)
    {
      histContainer_["cutflow"]->Fill(3);
      if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(3);
      if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(3);
         }
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
Eff::beginJob()
{
  edm::Service<TFileService> fs;

    histContainer_["cutflow"]   = fs->make<TH1F>("cutflow",    ";Selection cut;Events", 8, 0., 8.);
    histContainer_["ecutflow"]  = fs->make<TH1F>("ecutflow",   ";Selection cut;Events", 8, 0., 8.);
    histContainer_["mucutflow"] = fs->make<TH1F>("mucutflow",  ";Selection cut;Events", 8, 0., 8.);
    TString steps[]={"reco","=1 good lepton","=0 loose leptons","#geq2 jets","#geq3 jets","#geq4 jets","#geq1 b-tag","#geq2 b-tags"};
    for(size_t i=0; i<sizeof(steps)/sizeof(TString); i++)
       {
         histContainer_["cutflow"]   -> GetXaxis()->SetBinLabel(i+1,steps[i]);
         histContainer_["ecutflow"]  -> GetXaxis()->SetBinLabel(i+1,steps[i]);
         histContainer_["mucutflow"] -> GetXaxis()->SetBinLabel(i+1,steps[i]);
       }

    histContainer_["nrecomuons"] = fs->make<TH1F>("nrecomuons", ";# reconstructed muons; Events", 10, 0., 10.);
    histContainer_["nselmuons"]      = fs->make<TH1F>("nselmuons",       ";# reconstructed muons;Events",              5, 0.,5.);

    histContainer_["nrecoelectrons"]    = fs->make<TH1F>("nrecoelectrons",  ";# reconstructed electrons;Events",         5, 0., 5.);
    histContainer_["nselelectrons"]     = fs->make<TH1F>("nselelectrons",   ";# selected electrons;Events", 5, 0., 5.);

    histContainer_["nrecojets"]      = fs->make<TH1F>("nrecojets",   ";#reconstructed jets;Events", 50, 0., 50.);

    for(size_t i=0; i<4; i++)
       {
        TString pf(""); if(i==1) pf="b"; if(i==2) pf="light"; if(i==3) pf="other";
        histContainer_[(pf+"jetpt").Data()]     = fs->make<TH1F>(pf+"jetpt",       ";Transverse momentum [GeV];# jets", 100, 0., 250.);
        histContainer_[(pf+"chjetpt").Data()]   = fs->make<TH1F>(pf+"chjetpt",     ";Charged transverse momentum [GeV];# jets", 100, 0., 200.);
       }
    histContainer_["jeteta"]         = fs->make<TH1F>("jeteta",      ";Pseudo-rapidity;# jets", 100, 0., 3.);
    histContainer_["jetcsv"]         = fs->make<TH1F>("jetcsv",      ";Combined secondary vertes;# jets", 100, -1.2, 1.2);
    histContainer_["jetpileupid"]    = fs->make<TH1F>("jetpileupid", ";Pileup jet id;#jets", 100, -1.2, 1.2);
    histContainer_["nseljets"]       = fs->make<TH1F>("nseljets",    ";#selected jets;Events", 6, 3., 10.);
    histContainer_["nseljetsfull"]   = fs->make<TH1F>("nseljetsfull",    ";#selected jets;Events", 6, 3., 10.);

    histContainer_["nsvtx"]          = fs->make<TH1F>("nsvtx",    ";# secondary vertices;Events",5, 0., 5.);
    histContainer_["nvertices"] = fs->make<TH1F>("nvertices",    ";# vertices;Events", 100, 0., 100.);

    for(size_t ijet=2; ijet<=4; ijet++)
      for(size_t imet=0; imet<2; imet++)
         {
          TString prefix(imet==0 ? "": "ch");
          TString postfix(""); postfix += ijet;
          histContainer_[(prefix+"metpt"+postfix).Data()] = fs->make<TH1F>(prefix+"metpt"+postfix,    ";"+prefix+" Missing transverse energy [GeV];Events", 100, 0., 300.);
          histContainer_[(prefix+"metphi"+postfix).Data()] = fs->make<TH1F>(prefix+"metphi"+postfix,    ";"+prefix+" Missing transverse energy #phi [rad];Events", 50, -3.2, 3.2);
          histContainer_[(prefix+"mt"+postfix).Data()] = fs->make<TH1F>(prefix+"mt"+postfix,    ";"+prefix+" Transverse mass [GeV]; Events", 100, 0., 200.);
         }

    return;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Eff::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Eff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Eff);
