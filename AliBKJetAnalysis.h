#ifndef ALIBKJETANALYSIS_H
#define ALIBKJETANALYSIS_H
class TH1;
class TH1D;
class TH2;
class TH3;
class TF1;
class THnSparse;
class TRandom3;
class TGraphErrors;
class TProfile;
class THistManager;
class TLorentzVector;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliEmcalTrackSelection;
class AliAnalysisUtils;
class AliCalorimeterUtils;
class AliMultSelection;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"
#include <vector>

using namespace std;

class AliBKJetAnalysis : public AliAnalysisTaskEmcalJet {
public:
  typedef std::vector<TLorentzVector>   TLorentzVector1D;
  typedef std::vector<TLorentzVector1D> TLorentzVector2D;
  typedef std::vector<Double_t> Double1D; 
  typedef std::vector<Bool_t> Bool1D; 
  enum { kBTypeBegin=1,kBData=1, kBMixingSimple, kBTypeEnd };
  enum { kBDiJetSelBegin=1,kNoCut=1,kB1,kB2,kB3,kC1,kC2,kC3,kD1,kD2,kD3,kM1,kM2,kM3,kM4,kM5,kM6,kM7,kBDiJetSelEnd }; 
  AliBKJetAnalysis();
  AliBKJetAnalysis(const char *name);
  AliBKJetAnalysis(const char *name, const char *option);
  AliBKJetAnalysis(const AliBKJetAnalysis& ap);
  AliBKJetAnalysis& operator =(const AliBKJetAnalysis& ap);
  virtual void    FinishTaskOutput();
  ~AliBKJetAnalysis();
  virtual void UserCreateOutputObjects();
  void SetIsAA(bool is=true){ fIsAA = is; }
  void SetLeadingParticlePtMin(Double_t m){ fLeadingParticlePtMin=m; };
  void SetIsMC(Bool_t ismc){fIsMC = ismc;};
  void SetScalingFactorHist(TH1D* sfhist) {fScalingFactorHist = (TH1D*)sfhist->Clone();};
  void SetPtHardBin(double b) {pthardbin = b;};
  void RhoSparse(AliJetContainer *ktContainer, AliJetContainer *aktContainer, Int_t numberofexcludingjets, Bool_t isfulljet);
  Double1D &GetDijetPtPair() { return fDijetPtPair; };
  Double1D& GetDijetInvM(){return fDijetInvM;};
  Bool1D& GetDijetSelectionCut(){return fDijetSelectionCut;};
  TLorentzVector1D& GetJets(){return fJets;};
  Bool_t    GetIsGenGoodVtx(){return IsGenGoodVtx;};
protected:
  Bool_t       Run();
  TAxis AxisFix( TString name, int nbin, Double_t xmin, Double_t xmax);
  TAxis AxisVar( TString name, std::vector<Double_t> bin );
  TAxis AxisLog( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0);
  TAxis AxisStr( TString name, std::vector<TString> bin );
  THnSparse * CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt="");
  THnSparse * CreateTHnSparse(TString name, TString title, TString templ, Option_t * opt="");
  Long64_t FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w=1. );
  Long64_t FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w=1. );
	Bool_t MeasurePtHardBinScalingFactor ();
  Bool_t isOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);
  Bool_t MeasureJets(AliJetContainer *jetContainer, TLorentzVector1D &Jets, TLorentzVector1D &JetsBeforeCorr, Bool_t istruth, Bool_t isfulljet);
  void CheckDijetSelections(TLorentzVector1D Jets, TLorentzVector2D &sj, Bool1D &disel);
  Double_t getDiffR(double phi1, double phi2, double eta1, double eta2);
 
private:
  TList*                          fOutput = nullptr; //!
  TString                         fOption=""; 
  AliJetContainer*                fJetsCont = nullptr; //!
  UInt_t                          fFilterBit=768;
  Bool_t                          IsFirstEvent=1;
  THistManager*                   fHistos=nullptr; //!
  Bool_t                          fIsAA=0; 
  Int_t                           fNDiJetSelection=0; 
  Int_t                           fNType=0;
  Bool_t                          fIsMC=kFALSE; 
  Double_t                        fLeadingParticlePtMin=0; 
  TRandom3*                       fBSRandom=nullptr; //!
  Double1D                        fDijetPtPair;
  Double1D                        fDijetInvM;
  Bool1D                          fDijetSelectionCut;
  TLorentzVector1D                fJets;
  TF1                             *invmscale=nullptr;//!
  AliAnalysisUtils                *fUtils=nullptr;//!
  TString                         filename="";
  Double_t                        pyxsechistsf = 1;
  TH1D*                           fScalingFactorHist = nullptr;
  Bool_t                          IsGenGoodVtx = false;
  TAxis                           binJetTPtCut; //!
  TAxis                           binJetPtCut; //!
  TAxis                           binningpthard; //!
  AliVEvent *fEvt = nullptr;                           //!
  AliCalorimeterUtils*            fCaloUtils = nullptr;//!
  TClonesArray*                   fMCArray=nullptr; //!
  AliMultSelection*               sel=nullptr;//!
  TF1*                            tsf=nullptr;//!
  TF1*                            tsfl=nullptr;//!
  TF1*                            tsfh=nullptr;//!
  Double_t                        Rhopt = 0; 
  Double_t                        Rhom = 0;
  Double_t                        Rhoptfulljet = 0; 
  TLorentzVector                  RHOCHJET;
  TLorentzVector                  RHOFULLJET;
  Double_t                        sf = 1;
  Double_t                        genzvtx = -30;
  Int_t                           NTrials = -1;
  Double_t                        XSection = -1;
  TLorentzVector                  p6;
  TLorentzVector                  p7;
  Double_t                        vertex[3];
  Bool_t                          IsGoodVertex = false;
  Double_t                        pthardbin = 0.5; //first bin
  double                          fEventsPassed = 0;
  Double_t                        ptlead = 0;
  ClassDef(AliBKJetAnalysis, 1)
};
#endif
