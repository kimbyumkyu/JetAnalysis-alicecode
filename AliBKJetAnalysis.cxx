/*************************************************************************
 *
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Comment describing what this class does needed!

//==================================================================
// Class for di-charged and full jet analyses.
// by Beomkyu KIM
//==================================================================
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TList.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <THnSparse.h>
#include <TRandom.h>
#include "THistManager.h"
#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"
#include "AliCentrality.h"
#include "AliBKJetAnalysis.h"
#include <TRandom3.h>
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliEmcalPythiaInfo.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliAnalysisUtils.h"
//#include "AliAnalysisTaskCounter.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliCalorimeterUtils.h"
#include "AliEMCALGeometry.h"
#include "AliMultSelection.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskRhoSparse.h"
const Double_t pionmass = AliPID::ParticleMass(AliPID::kPion);
const Double_t pi = TMath::Pi();

AliBKJetAnalysis::AliBKJetAnalysis()
	: AliAnalysisTaskEmcalJet("AliBKJetAnalysis", kTRUE), fOutput(0)
{
	DefineOutput(1, TList::Class());
}

AliBKJetAnalysis::AliBKJetAnalysis(const char *name)
	: AliAnalysisTaskEmcalJet(name, kTRUE), fOutput(0)
{
	DefineOutput(1, TList::Class());
}
AliBKJetAnalysis::AliBKJetAnalysis(const char *name, const char *option)
	: AliAnalysisTaskEmcalJet(name, kTRUE), fOutput(0), fOption(option)
{
	DefineOutput(1, TList::Class());
}

AliBKJetAnalysis::AliBKJetAnalysis(const AliBKJetAnalysis &ap)
	: AliAnalysisTaskEmcalJet(ap.fName, kTRUE), fOutput(ap.fOutput)
{
}
AliBKJetAnalysis &AliBKJetAnalysis::operator=(const AliBKJetAnalysis &ap)
{

	this->~AliBKJetAnalysis();
	new (this) AliBKJetAnalysis(ap);
	return *this;
}

AliBKJetAnalysis::~AliBKJetAnalysis()
{
	delete fOutput;
	delete fBSRandom;
}

void AliBKJetAnalysis::UserCreateOutputObjects()
{
	fBSRandom = new TRandom3;
	fBSRandom->SetSeed();
	//===============================
	// BINS
	//===============================
	auto binJetMass = AxisFix("JetMass", 50, 0, 100);
	auto binInvM = AxisVar("InvM", {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 160, 170, 180, 190, 200, 220, 240, 270, 300, 500, 1000, 5000});
	auto binCent = AxisFix("Cent", 1, 0, 100);
	auto binDiJetSel = AxisStr("DiJetSel", {"NoCut", "B1", "B2", "B3", "C1", "C2", "C3", "D1", "D2", "D3", "M1", "M2", "M3", "M4", "M5", "M6", "M7"});
	fNDiJetSelection = binDiJetSel.GetNbins();
	if (fIsAA)
	{
		binCent = AxisFix("Cent", 20, 0, 100);
	}
	else
	{
		binCent = AxisVar("Cent", {0, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100});
	}
	auto binLog1k = AxisLog("Log1k", 500, 0.1, 1000, 0);
	//auto binLog3c     = AxisVar("Log3c",{0,5,7,9,12,16,21,28,36,45,57,70,85,99,115,132,150,169,190,212,235,500});
	auto binLog3c = AxisVar("Log3c", {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 160, 170, 180, 190, 200, 220, 240, 270, 300, 500, 1000, 5000});
	auto bintpt = AxisVar("JetTPt", {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 160, 170, 180, 190, 200, 220, 240, 270, 300, 500, 1000, 5000});
	auto bin1c = AxisFix("Fix1c", 100, 0, 100);
	auto binAsim = AxisFix("Asim", 100, 0, 1);
	auto binM = AxisFix("BinM", 600, -300, 300);
	auto binpthardbin = AxisFix("pthardbin", fOption.Contains("LHC13")? 10 : 20, 0, fOption.Contains("LHC13")? 10 : 20);
	auto binjetpt = AxisVar("binjetpt", {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 65, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 5000});
	auto binjt = AxisLog("JtBin", 40, 0.01, 10, 0);
	auto binz = AxisFix("zbin", 5, 0, 1);
	auto binjetmultiplicity = AxisFix("mbin", 50, 0, 50);
	auto binjes = AxisFix("jesbin", 100, -1, 1);

	//https://alimonitor.cern.ch/users/download.jsp?view=true&path=/alice/cern.ch/user/a/aliprod/bin/simrun.sh
	if (fOption.Contains("LHC13"))
		binningpthard = AxisVar("binpthard", {0, 5, 11, 21, 36, 57, 84, 117, 152, 191, 234});
	else 
		binningpthard = AxisVar("binprhard", {5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235, 10000});

	//===============================
	// HISTOGRAMS
	//===============================
	fHistos = new THistManager("jethists");

	const int nbins = 100;
	Double_t logbins[nbins + 1];
	Double_t low = 0.1;
	Double_t high = 300;
	Double_t logbw = (log(high) - log(low)) / nbins;
	for (int ij = 0; ij <= nbins; ij++)
		logbins[ij] = low * exp(ij * logbw);

	fHistos->CreateTH1("zvtx", "zvtx", 60, -30, 30, "s");
	fHistos->CreateTH1("mczvtx", "mczvtx", 60, -30, 30, "s");
	fHistos->CreateTH1("hCent", "histcent", 100, 0, 100, "s");
	fHistos->CreateTH1("nCentrality", "nCentrality", 10, 0, 100, "s");
	fHistos->CreateTH2("trketaphi", "trketaphi", 20, -1, 1, 90, 0, TMath::TwoPi(), "s");
	fHistos->CreateTH2("jetetaphi", "jetetaphi", 20, -1, 1, 90, 0, TMath::TwoPi(), "s");
	fHistos->CreateTH1("fulljetpt", "fulljetpt", 300, 0, 300, "s");
	fHistos->CreateTH1("fulljettruept", "fulljettruept", 300, 0, 300, "s");
	fHistos->CreateTH1("fulljeteta", "fulljeteta", 40, -2, 2, "s");
	fHistos->CreateTH1("z4060", "z4060", 2000, 0, 1, "s");
	fHistos->CreateTH1("z6080", "z6080", 2000, 0, 1, "s");
	fHistos->CreateTH1("z80100", "z80100", 2000, 0, 1, "s");
	fHistos->CreateTH1("z100150", "z100150", 2000, 0, 1, "s");
	//vector<TString> ent = {"All", "PassPileUp", "PassPileUpGoodz", "GoodzNTrials", "GoodzNX", "NXoNTrials"};
	vector<TString> ent = {"All", "PassPileUp", "PassPileUpGoodz", "GoodzNTrials", "GoodzNX", "NXoNTrials", "AllBeforePassPileUp"};
	auto h = fHistos->CreateTH1("hEventNumbers", "", ent.size(), 0, ent.size());
	for (auto i = 0u; i < ent.size(); i++)
		h->GetXaxis()->SetBinLabel(i + 1, ent.at(i).Data());

	CreateTHnSparse("hJetPtLeading", "", 4, {binDiJetSel, binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hRho", "Rho dist", 3, {binCent, bin1c, binpthardbin}, "s");
	CreateTHnSparse("hDiJetDPhi_0_2pi", "DiJet #Delta#Phi", 5, {binDiJetSel, binCent, binjetpt, AxisFix("", 100, 0, 2 * pi), binpthardbin}, "s");
	CreateTHnSparse("hDiJetDPhi_0_2piTruth", "DiJet #Delta#Phi", 5, {binDiJetSel, binCent, binjetpt, AxisFix("", 100, 0, 2 * pi), binpthardbin}, "s");
	CreateTHnSparse("hDiJetInvMPtPair", "DiJet PtPair", 5, {binDiJetSel, binCent, binjetpt, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hDiJetInvMPtPairTruth", "DiJet PtPair Truth", 5, {binDiJetSel, binCent, binjetpt, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hDiJetInvMPtPairRes", "DiJet InvM PtPair Res Matrix", 7, {binDiJetSel, binCent, binjetpt, binjetpt, binjetpt, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hDiJetInvMPtPairMiss", "DiJet InvM PtPair missing tracks", 5, {binDiJetSel, binCent, binjetpt, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hDiJetInvMPtPairFake", "DiJet InvM PtPair fake tracks", 5, {binDiJetSel, binCent, binjetpt, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hPtHardBin", "pthardbin", 2, {binCent, AxisFix("", 20, 0, 20)}, "s");

	CreateTHnSparse("hJetPt", "Inclusive jet pt", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetPtTest", "Inclusive jet pt", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetPtBeforeCorr", "Inclusive jet pt", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetPtBeforeMatching", "Inclusive matched jet pt", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetPtTruth", "Inclusive Gen jet pt", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetPtRes", "Inclusive jet pt Res Matrix", 4, {binCent, binjetpt, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetPtFake", "Inclusive jet pt fake tracks", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetPtMiss", "Inclusive jet pt  missing tracks", 3, {binCent, binjetpt, binpthardbin}, "s");

	CreateTHnSparse("hFullJetPt", "Inclusive jet pt", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hTrueFullJetPt", "Inclusive jet pt", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hTrueFullJetPtMatched", "Inclusive jet pt", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hFullJetMultiplicity", "Inclusive jet pt", 3, {binCent, binjetmultiplicity, binpthardbin}, "s");
	CreateTHnSparse("hTrueFullJetMultiplicity", "Inclusive jet pt", 3, {binCent, binjetmultiplicity, binpthardbin}, "s");
	CreateTHnSparse("hFullJetJES", "Full jet pt jes", 4, {binCent, binjetpt, binjes,binpthardbin}, "s");

	CreateTHnSparse("hJetJt", "Inclusive jet jt", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");
	CreateTHnSparse("hJetJtWeight", "Inclusive jet jt over jt", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");
	CreateTHnSparse("hPerpJtWeight", "Perpendicular jt over jt", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");
	CreateTHnSparse("hJetJtBeforeCorr", "Inclusive jet jt before correction", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");
	CreateTHnSparse("hJetJtBeforeMatching", "Jet jt before mathcing", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");
	CreateTHnSparse("hJetJtTruth", "Inclusive Gen jet jt", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");
	CreateTHnSparse("hJetJtTruthWeight", "Inclusive Gen jet jt ove jt", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");
	CreateTHnSparse("hJetJtRes", "Inclusive jet jt Res Matrix", 7, {binCent, binjetpt, binjetpt, binz, binjt, binjt, binpthardbin}, "s");
	CreateTHnSparse("hJetJtFake", "Inclusive jet jt fake tracks", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");
	CreateTHnSparse("hJetJtMiss", "Inclusive jet jt missing tracks", 5, {binCent, binjetpt, binz, binjt, binpthardbin}, "s");

	//test purpose
	CreateTHnSparse("hJetMass", "Inclusive jet mass", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetMassRes", "Inclusive jet mass  Res Matrix", 4, {binCent, binjetpt, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetMassTest", "Inclusive jet mass", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetMassBeforeCorr", "Inclusive jet mass", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetMassBeforeMatching", "Inclusive matched jet mass", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetMassTruth", "Inclusive Gen jet mass", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetMassMiss", "Inclusive Gen jet mass missed", 3, {binCent, binjetpt, binpthardbin}, "s");
	CreateTHnSparse("hJetMassFake", "Inclusive Gen jet mass faked", 3, {binCent, binjetpt, binpthardbin}, "s");

	PostData(1, fHistos->GetListOfHistograms());

	//===============================
	// For Sure
	//===============================
	//std::cout<< "DEBUG4 IsAA?"<< (fIsAA?"AA":"pp")<<std::endl;
	//std::cout<<"NBins of Cent : "<<binCent.GetNbins()<<"\t"<<binCent.GetXmin()<<"\t"<<binCent.GetXmax()<<endl;
	if (fNDiJetSelection != kBDiJetSelEnd - 1)
	{
		cout << "fNDiJetSelection(" << fNDiJetSelection
			 << ") is not match with kBDiJetSelEnd(" << kBDiJetSelEnd << ")" << endl;
		gSystem->Exit(1);
	}
	for (auto i = 0u; i < kBDiJetSelEnd; i++)
	{
		fDijetSelectionCut.push_back(false);
		fDijetPtPair.push_back(0.);
		fDijetInvM.push_back(0.);
	}
	fUtils = new AliAnalysisUtils();
	fUtils->SetMaxVtxZ(10);
}

//________________________________________________________________________

Bool_t AliBKJetAnalysis::Run()
{
	using TMath::Abs;
	fCent = -1;
	//Cent = InputEvent()->GetCentrality();
	sel = (AliMultSelection *)InputEvent()->FindListObject("MultSelection");
	if (sel)
	{
		if (fOption.Contains("13f"))
			fCent = sel->GetMultiplicityPercentile("V0C");
		else if (
			fOption.Contains("13b") ||
			fOption.Contains("13c") ||
			fOption.Contains("13d") ||
			fOption.Contains("13e"))
			fCent = sel->GetMultiplicityPercentile("V0A");
		else
			fCent = sel->GetMultiplicityPercentile("V0M");
	}
	fHistos->FillTH1("hCent", fCent);

	AliInputEventHandler *inputHandler = (AliInputEventHandler *)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();

	Bool_t IsMinimumBias = kFALSE;
	if (fOption.Contains("LHC13d") || fOption.Contains("LHC13e"))
		IsMinimumBias = inputHandler->IsEventSelected() & AliVEvent::kEMCEJE;
	else
		IsMinimumBias = inputHandler->IsEventSelected() & AliVEvent::kINT7;

	if (fOption.Contains("MC"))
		IsMinimumBias = inputHandler->IsEventSelected() & AliVEvent::kAny;
	//IsMinimumBias = inputHandler->IsEventSelected() & AliVEvent::kEMCEJE;
	//if (fOption.Contains("MC") || fOption.Contains("Emb"))
	//	IsMinimumBias = inputHandler->IsEventSelected() & AliVEvent::kAny;
	fHistos->FillTH1("hEventNumbers", "All", 1);
	if (!IsMinimumBias)
	{
		PostData(1, fHistos->GetListOfHistograms());
		return false;
	}

	fHistos->FillTH1("hEventNumbers", "AllBeforePassPileUp", 1);
	//pt hard bin scaling-----------------------------------------------------------
	//https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JetMCProductionsCrossSections
	sf = 1.;
	pthardbin = 0.5;
	if (fIsMC)
	{
		Bool_t rejecttailsOK = this->MeasurePtHardBinScalingFactor();
		if (!rejecttailsOK)
			return false;
	}
	//cout<<"Scaling factor = "<<sf<<endl;

	AliVEvent *event = InputEvent();
	event->IsA() == AliESDEvent::Class()
		? fEvt = dynamic_cast<AliESDEvent *>(event)
		: fEvt = dynamic_cast<AliAODEvent *>(event);
	if (!fEvt)
		return false;

	if (!fIsMC && !fOption.Contains("LHC15o") && fUtils->IsPileUpSPD(InputEvent()))
	{
		PostData(1, fHistos->GetListOfHistograms());
		return false;
	}

	fHistos->FillTH1("hEventNumbers", "PassPileUp", 1);

	// Fill z_vertex and cut z_vertex range
	IsGoodVertex = false;
	IsGenGoodVtx = false;
	if (fOption.Contains("Emb") || fOption.Contains("MC"))
	{
		fHistos->FillTH1("zvtx", genzvtx);
		if (fabs(genzvtx) <= 10)
		{
			IsGenGoodVtx = true;
		}
	}

	const AliVVertex *vtx = InputEvent()->GetPrimaryVertex();
	if (fUtils->IsVertexSelected2013pA(InputEvent()))
	{
		//if (1){
		if (vtx->GetNContributors() > 1)
		{
			double zvtx = vtx->GetZ();
			vtx->GetXYZ(vertex);
			fHistos->FillTH1("zvtx", zvtx);
			if (fabs(zvtx) <= 10)
				IsGoodVertex = true;
		}
		//if (IsGoodVertex) cout<<"bkkim zvtx : "<<zvtx<<endl;
	}

	if (IsGoodVertex)
	{
		fHistos->FillTH1("hEventNumbers", "PassPileUpGoodz", 1);
		fHistos->FillTH1("hEventNumbers", "GoodzNTrials", NTrials);
		fHistos->FillTH1("hEventNumbers", "GoodzNX", XSection);
		fHistos->FillTH1("hEventNumbers", "NXoNTrials", XSection / NTrials);
		fHistos->FillTH1("nCentrality", fCent);
		FillTHnSparse("hPtHardBin", {fCent, pthardbin}, 1);
	}

	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	//===============================
	// RESUM JETS
	//===============================

	auto jetContainer = GetJetContainer(0);
	//cout<<"First jet Container name : "<<jetContainer -> GetName()<<endl;
	auto trkContainer = GetTrackContainer(0); //trk container recall for QA purpose
	//cout<<"trkContainer name "<<trkContainer->GetName()<<endl;
	auto cluContainer = GetClusterContainer(0);
	auto ktContainer = GetJetContainer(1);
	auto fulljetContainer = GetJetContainer(2);
	auto ktfulljetContainer = GetJetContainer(3);

	AliJetContainer *mcContainer = nullptr;
	AliJetContainer *mcFullJetContainer = nullptr;
	AliParticleContainer *mcparticleContainer = nullptr;
	if (fOption.Contains("Emb") || fOption.Contains("MC"))
	{
		mcContainer = GetJetContainer(4);
		mcFullJetContainer = GetJetContainer(5);
		mcparticleContainer = GetParticleContainer(1);
		//cout<<"MC jet Container name : "<<mcContainer -> GetName()<<endl;
	}

	//AliAnalysisTaskRhoSparse *fRho = (AliAnalysisTaskRhoSparse *)(AliAnalysisManager::GetAnalysisManager()->GetTask("AliAnalysisTaskRhoSparse"));
	//if (jetContainer->GetRhoParameter())
	//cout << "Rho name " << ktContainer->GetRhoName() << endl;
	//AliRhoParameter *rhoParam = dynamic_cast<AliRhoParameter *>(InputEvent()->FindListObject("Rho"));
	//cout << "RhoValue = " << rhoParam->GetVal() << endl;
	//cout << "RhoValue2 = " << jetContainer->GetRhoVal() << endl;
	cout << "Rho value = " << RHO << endl;
	FillTHnSparse("hRho", {fCent, RHO, pthardbin}, sf);
	// =============================
	// QA plots
	// ============================

	//simple full jet test start

	//


	for (auto trk : trkContainer->all())
	{
		//fHistos->FillTH2("trketaphi", trk->Eta(), trk->Phi());
		//cout<<"trk -> Eta : "<<trk->Eta()<<endl;
		//cout<<"trk -> Label : "<<trk->GetLabel()<<endl;
	}
	for (int itrack = 0; itrack < trkContainer->GetNParticles(); itrack++)
	{
		AliAODTrack *track = static_cast<AliAODTrack *>(trkContainer->GetParticle(itrack));
		if (!track)
			continue;
		//fHistos->FillTH2("trketaphi", track->Eta(), track->Phi());
	}

	//this->MeasureBgDensity(ktContainer);
	this->RhoSparse(ktContainer, jetContainer, 2);

	TLorentzVector1D RecJets;
	TLorentzVector1D RecJetsBeforeCorr;

	Bool_t isjetok = this->MeasureJets(jetContainer, RecJets, RecJetsBeforeCorr, false, false);
	if (!isjetok)
	{
		PostData(1, fHistos->GetListOfHistograms());
		return false;
	}
	std::sort(RecJets.begin(), RecJets.end(), [&](const TLorentzVector &x, const TLorentzVector &y) { return x.Pt() > y.Pt(); });

	TLorentzVector1D TrueJets;
	TLorentzVector1D TrueJetsBeforeCorr;
	TLorentzVector1D TrueFullJets;
	TLorentzVector1D TrueFullJetsBeforeCorr;

	if (fOption.Contains("Emb") || fOption.Contains("MC"))
	{
		this->MeasureJets(mcContainer, TrueJets, TrueJetsBeforeCorr, true, false);
		std::sort(TrueJets.begin(), TrueJets.end(), [&](const TLorentzVector &x, const TLorentzVector &y) { return x.Pt() > y.Pt(); });
		this->MeasureJets(mcFullJetContainer, TrueFullJets, TrueFullJetsBeforeCorr, true, true);
		std::sort(TrueFullJets.begin(), TrueFullJets.end(), [&](const TLorentzVector &x, const TLorentzVector &y) { return x.Pt() > y.Pt(); });
	}

	TLorentzVector1D matchedjets;
	if ((fOption.Contains("Emb") || fOption.Contains("MC")))
	{ // Inclusive jet pt response matrix
		if (IsGenGoodVtx)
		{
			TLorentzVector1D pjets = TrueJets;
			TLorentzVector1D rjets = RecJets;
			for (auto pj : pjets)
			{
				FillTHnSparse("hJetPtTruth", {fCent, pj.Pt(), pthardbin}, sf);
				FillTHnSparse("hJetMassTruth", {fCent, pj.M(), pthardbin}, sf);
			}

			if (IsGoodVertex)
			{
				for (auto rj : rjets)
				{
					FillTHnSparse("hJetPtBeforeMatching", {fCent, rj.Pt(), pthardbin}, sf);
					FillTHnSparse("hJetMassBeforeMatching", {fCent, rj.M(), pthardbin}, sf);
				}
				for (auto rj : RecJetsBeforeCorr)
				{
					FillTHnSparse("hJetPtBeforeCorr", {fCent, rj.Pt(), pthardbin}, sf);
					FillTHnSparse("hJetMassBeforeCorr", {fCent, rj.M(), pthardbin}, sf);
				}

				for (auto pj : pjets)
				{
					TLorentzVector maxjet(0, 0, 0, 0);
					Double_t leastdr = 10;
					for (auto rj : rjets)
					{
						if (rj.DeltaR(pj) < 0.24 && maxjet.Pt() < rj.Pt())
						{
							maxjet = rj;
						}
					}
					if (maxjet.Pt() > 0)
					{
						matchedjets.push_back(maxjet);
						FillTHnSparse("hJetPt", {fCent, maxjet.Pt(), pthardbin}, sf);
						rjets.erase(std::remove_if(rjets.begin(), rjets.end(), [&](const TLorentzVector &x) { return x.Pt() == maxjet.Pt(); }), rjets.end());
						FillTHnSparse("hJetPtRes", {fCent, maxjet.Pt(), pj.Pt(), pthardbin}, sf);
						FillTHnSparse("hJetMassRes", {fCent, maxjet.M(), pj.M(), pthardbin}, sf);
					}
					if (maxjet.Pt() == 0)
					{
						FillTHnSparse("hJetPtMiss", {fCent, pj.Pt(), pthardbin}, sf);
						FillTHnSparse("hJetMassMiss", {fCent, pj.M(), pthardbin}, sf);
					}
				}
				if (!fIsAA)
				{
					for (auto rj : rjets)
					{
						FillTHnSparse("hJetPtFake", {fCent, rj.Pt(), pthardbin}, sf);
						FillTHnSparse("hJetMassFake", {fCent, rj.M(), pthardbin}, sf);
					}
				}
				TLorentzVector1D RecJetsMatched = matchedjets;
				//cout<<"matchedjets"<<endl;
				for (auto mj : RecJetsMatched)
				{
					//FillTHnSparse("hJetPt", {fCent, mj.Pt(), pthardbin}, sf);
					FillTHnSparse("hJetMass", {fCent, mj.M(), pthardbin}, sf);
				}
				RecJets = RecJetsMatched;
			}
			else
			{
				for (auto pj : pjets)
				{
					FillTHnSparse("hJetPtMiss", {fCent, pj.Pt(), pthardbin}, sf);
					FillTHnSparse("hJetMassMiss", {fCent, pj.M(), pthardbin}, sf);
				}
			}
		}
	}
	else if (IsGoodVertex)
	{
		for (auto rj : RecJets)
		{
			//FillTHnSparse("hJetPt", {fCent, rj.Pt(), pthardbin}, sf);
			FillTHnSparse("hJetMass", {fCent, rj.M(), pthardbin}, sf);
		}
		for (auto rj : RecJetsBeforeCorr)
		{
			FillTHnSparse("hJetPtBeforeCorr", {fCent, rj.Pt(), pthardbin}, sf);
			FillTHnSparse("hJetMassBeforeCorr", {fCent, rj.M(), pthardbin}, sf);
		}
	}

	TLorentzVector2D sj(fNDiJetSelection + 1, TLorentzVector1D(2));
	Bool1D recdisel(fNDiJetSelection + 1, false);
	this->CheckDijetSelections(RecJets, sj, recdisel);
	TLorentzVector2D sjkine(fNDiJetSelection + 1, TLorentzVector1D(2));
	Bool1D truedisel(fNDiJetSelection + 1, false);
	if (fOption.Contains("Emb") || fOption.Contains("MC"))
	{
		this->CheckDijetSelections(TrueJets, sjkine, truedisel);
	}

	//cout <<"sj size : "<<sj.size()<<endl;
	//===============================
	// ptPair : leading - subleading
	//===============================

	for (int ids = kBDiJetSelBegin; ids < kBDiJetSelEnd; ids++)
	{ // NOTE : begin with 1
		auto j = sj[ids];
		//=== SKIP Empty DiJet
		//fDijetSelectionCut[ids] = false;
		if (fOption.Contains("Emb") || fOption.Contains("MC"))
		{
			auto truej = sjkine[ids];
			auto dijet = truej[0] + truej[1];
			auto invM = dijet.M();
			auto ptpair = dijet.Pt();
			Bool_t truthdijetcut = truedisel.at(ids);
			if (IsGenGoodVtx && truthdijetcut && IsGoodVertex)
			{
				if (!recdisel[ids])
					FillTHnSparse("hDiJetInvMPtPairMiss", {(double)ids, fCent, invM, ptpair, pthardbin}, sf);
			}

			if (truthdijetcut && IsGenGoodVtx)
			{
				FillTHnSparse("hDiJetInvMPtPairTruth", {double(ids), fCent, invM, ptpair, pthardbin}, sf);
				FillTHnSparse("hDiJetDPhi_0_2piTruth", {double(ids), fCent, invM, TVector2::Phi_0_2pi(truej[0].DeltaPhi(truej[1])), pthardbin}, sf);
			}
		}
		if (!recdisel[ids])
			continue;
		Double_t nphi = 0;
		auto diJetSel = Double_t(ids);
		auto dijet = j[0] + j[1];
		auto invM = dijet.M();
		auto dPhi = j[0].DeltaPhi(j[1]);
		auto dPhiA = fabs(dPhi);
		auto dPhi_0_2pi = TVector2::Phi_0_2pi(dPhi);
		auto tpt = j[0].Pt();
		auto apt = j[1].Pt();
		auto ptpair = dijet.Pt();
		auto testdphi = nphi - j[0].Phi();
		auto ptAsim = (tpt - apt) / (tpt + apt);
		auto eAsim = (j[0].E() - j[1].E()) / (j[0].E() + j[1].E());
		auto kty = j[0].Pt() * TMath::Sin(dPhi_0_2pi);
		if (ids >= kM1)
			ptpair = kty;

		double tratio = 1.;
		double tratioh = 1.;
		double tratiol = 1.;

		if (fOption.Contains("Emb") || fOption.Contains("MC"))
		{
			auto truej = sjkine[ids];
			auto truedijet = truej[0] + truej[1];
			Double_t truthptpair = truedijet.Pt();
			Double_t truthinvM = truedijet.M();
			Bool_t truthdijetcut = truedisel.at(ids);
			auto dPhi_0_2piTruth = TVector2::Phi_0_2pi(truej[0].DeltaPhi(truej[1]));

			if (IsGenGoodVtx)
			{
				if (IsGoodVertex)
				{
					if (truthdijetcut)
					{
						FillTHnSparse("hDiJetInvMPtPairRes", {diJetSel, fCent, invM, truthinvM, ptpair, truthptpair, pthardbin}, sf);
					}
					else
					{
						FillTHnSparse("hDiJetInvMPtPairFake", {diJetSel, fCent, invM, ptpair, pthardbin}, sf);
						//FillTHnSparse("hDiJetInvMPtPairMiss", {diJetSel, fCent, truthinvM, truthptpair, pthardbin}, sf);
					}
				}
				else if (truthdijetcut)
					FillTHnSparse("hDiJetInvMPtPairMiss", {diJetSel, fCent, truthinvM, truthptpair, pthardbin}, sf);
			} // end of IsGenGoodVtx
		}	  // end of Emb and MC
		if (IsGoodVertex)
		{
			FillTHnSparse("hJetPtLeading", {diJetSel, fCent, tpt, pthardbin}, sf * tratio);
			FillTHnSparse("hDiJetDPhi_0_2pi", {diJetSel, fCent, invM, dPhi_0_2pi, pthardbin}, sf);
			FillTHnSparse("hDiJetInvMPtPair", {diJetSel, fCent, invM, ptpair, pthardbin}, sf);
		}
	}

	for (auto ij : fulljetContainer->accepted_momentum())
	{
		auto j = ij.second;
		TLorentzVector sum(0, 0, 0, 0);
		double sumpt = 0;
		double lpt = 0;
		for (int it = 0; it < j->GetNumberOfTracks(); it++)
		{

			auto trk = j->Track(it);
			if (trk->Pt() > 100)
				continue;
			TLorentzVector temp;
			temp.SetXYZM(trk->Px(), trk->Py(), trk->Pz(), pionmass);
			sum += temp;
			sumpt += trk->Pt();
			if (lpt < temp.Pt())
				lpt = temp.Pt();
		}
		for (int it = 0; it < j->GetNumberOfClusters(); it++)
		{
			auto clu = j->Cluster(it);
			TLorentzVector temp;
			clu->GetMomentum(temp, vertex);
			sum += temp;
			sumpt += temp.Pt();
			if (lpt < temp.Pt())
				lpt = temp.Pt();
		}
		if (lpt < fLeadingParticlePtMin)
			continue;
		auto jetpt = sum.Pt();
		for (int it = 0; it < j->GetNumberOfTracks(); it++)
		{
			auto trk = j->Track(it);
			if (trk->Pt() > 100)
				continue;
			if (jetpt > 40 && jetpt < 60)
				fHistos->FillTH1("z4060", trk->Pt() / jetpt);
			if (jetpt > 60 && jetpt < 80)
				fHistos->FillTH1("z6080", trk->Pt() / jetpt);
			if (jetpt > 80 && jetpt < 100)
				fHistos->FillTH1("z80100", trk->Pt() / jetpt);
			if (jetpt > 100 && jetpt < 150)
				fHistos->FillTH1("z100150", trk->Pt() / jetpt);
		}
	}

	TLorentzVector1D RecFullJets;
	TLorentzVector1D RecFullJetsBeforeCorr;
	Double_t z = 0, jt = 0, deltaR = 10;
	auto radius = fulljetContainer->GetJetRadius();
	TLorentzVector T;
	this->MeasureJets(fulljetContainer, RecFullJets, RecFullJetsBeforeCorr, false, true);
	if (IsGoodVertex)
	{
		for (auto jet : RecFullJets){
			for (int itrack = 0; itrack < trkContainer->GetNParticles(); itrack++)
			{
				AliAODTrack *track = static_cast<AliAODTrack *>(trkContainer->GetParticle(itrack));
				if (!track)
					continue;
				if (fabs(track->Eta()) > (0.25+radius))
					continue;
				if	(!(((AliAODTrack*) track)->TestFilterBit(768)))
					continue; // primary charged track selection
				if (track->Pt() < 0.15)
					continue;
				T.SetXYZM(track->Px(), track->Py(), track->Pz(), pionmass);
				fHistos->FillTH2("trketaphi", track->Eta(), track->Phi());
				deltaR = getDiffR(jet.Phi(), track->Phi(), jet.Eta(), track->Eta());
				if (deltaR>radius)
					continue;
				z = (T.Vect() * jet.Vect().Unit()) / jet.P();
				jt = (T.Vect() - z * jet.Vect()).Mag();
				FillTHnSparse("hJetJt", {fCent, jet.Pt(), z, jt, pthardbin}, sf);
				FillTHnSparse("hJetJtWeight", {fCent, jet.Pt(), z, jt, pthardbin}, sf/jt);
			}
			// Perpendicular bg jt
			TLorentzVector vOrtho;
			vOrtho.SetVect(jet.Vect());
			vOrtho.SetE(jet.E());
			vOrtho.SetPhi(jet.Phi() + TMath::Pi() / 2);
			for (int itrack = 0; itrack < trkContainer->GetNParticles(); itrack++)
			{
				AliAODTrack *track = static_cast<AliAODTrack *>(trkContainer->GetParticle(itrack));
				if (!track)
					continue;
				if (fabs(track->Eta()) > 0.9)
					continue;
				if	(!(((AliAODTrack*) track)->TestFilterBit(768)))
					continue; // primary charged track selection
				if (track->Pt() < 0.15)
					continue;
				T.SetXYZM(track->Px(), track->Py(), track->Pz(), pionmass);
				deltaR = getDiffR(vOrtho.Phi(), track->Phi(), vOrtho.Eta(), track->Eta());
				if (deltaR>radius)
					continue;
				z = (T.Vect() * vOrtho.Vect().Unit()) / vOrtho.P();
				jt = (T.Vect() - z * vOrtho.Vect()).Mag();
				FillTHnSparse("hPerpJtWeight", {fCent, vOrtho.Pt(), z, jt, pthardbin}, sf/jt);
			}

			

		}
	}
	matchedjets.clear();
	Double_t mcz = 0, mcjt = 0;
	TLorentzVector MCT(0,0,0,0);
	if ((fOption.Contains("Emb") || fOption.Contains("MC")))
	{
		if (IsGenGoodVtx)
		{
			TLorentzVector1D pjets = TrueFullJets;
			TLorentzVector1D rjets = RecFullJets;
			FillTHnSparse("hTrueFullJetMultiplicity", {fCent, double(pjets.size()), pthardbin}, sf);
			for (auto pj : pjets)
			{

				FillTHnSparse("hTrueFullJetPt", {fCent, pj.Pt(), pthardbin}, sf);
				TLorentzVector maxjet(0, 0, 0, 0);
				Double_t leastdr = 10;
				for (auto rj : rjets)
				{
					if (rj.DeltaR(pj) < 0.24 && maxjet.Pt() < rj.Pt())
					{
						maxjet = rj;
					}
				}
				if (maxjet.Pt() > 0)
				{
					matchedjets.push_back(maxjet);

					FillTHnSparse("hTrueFullJetPtMatched", {fCent, pj.Pt(), pthardbin}, sf);
					FillTHnSparse("hFullJetPt", {fCent, maxjet.Pt(), pthardbin}, sf);
					FillTHnSparse("hFullJetJES", {fCent, pj.Pt(),(maxjet.Pt()-pj.Pt())/pj.Pt(), pthardbin}, sf);

					rjets.erase(std::remove_if(rjets.begin(), rjets.end(), [&](const TLorentzVector &x) { return x.Pt() == maxjet.Pt(); }), rjets.end());
					for (int itrack = 0; itrack < trkContainer->GetNParticles(); itrack++)
					{
						AliAODTrack *track = static_cast<AliAODTrack *>(trkContainer->GetParticle(itrack));
						if (!track)
							continue;
						if (fabs(track->Eta()) > (0.3 + radius))
							continue;
						if (!(((AliAODTrack *)track)->TestFilterBit(768)))
							continue; // primary charged track selection
						if (track->Pt() < 0.15)
							continue;
						T.SetXYZM(track->Px(), track->Py(), track->Pz(), pionmass);
						deltaR = getDiffR(maxjet.Phi(), track->Phi(), maxjet.Eta(), track->Eta());
						if (deltaR > radius)
							continue;
						z = (T.Vect() * maxjet.Vect().Unit()) / maxjet.P();
						jt = (T.Vect() - z * maxjet.Vect()).Mag();
						AliAODMCParticle *mcpart = dynamic_cast<AliAODMCParticle *>(mcparticleContainer->GetParticle(track->GetLabel()));
						if (!mcpart)
							continue;

						MCT.SetXYZT(mcpart->Px(), mcpart->Py(), mcpart->Pz(), mcpart->E());
						mcz = (MCT.Vect() * pj.Vect().Unit()) / pj.P();
						mcjt = (MCT.Vect() - mcz * pj.Vect()).Mag();
						FillTHnSparse("hJetJtRes", {fCent, maxjet.Pt(), pj.Pt(), z, jt, mcjt, pthardbin}, sf);
					}
				}
				FillTHnSparse("hFullJetMultiplicity", {fCent, double(matchedjets.size()), pthardbin}, sf);
			}
		}
	}

	PostData(1, fHistos->GetListOfHistograms());
	//PostData(1, fOutput);
	//cout<<"\n\n\n"<<endl;
	return kTRUE;
}
void AliBKJetAnalysis::FinishTaskOutput()
{
}
Bool_t AliBKJetAnalysis::isOverlapping(AliEmcalJet *jet1, AliEmcalJet *jet2)
{
	for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
	{
		Int_t jet1Track = jet1->TrackAt(i);
		for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
		{
			Int_t jet2Track = jet2->TrackAt(j);
			if (jet1Track == jet2Track)
				return kTRUE;
		}
	}
	return kFALSE;
}

void AliBKJetAnalysis::RhoSparse(AliJetContainer *ktContainer, AliJetContainer * aktContainer , Int_t numberofexcludingjets) {
	// Lets exclude a dijet
	AliEmcalJet *leading = nullptr;
	AliEmcalJet *subleading = nullptr;
	Int_t n = 0;
	//for (auto ij : aktContainer->accepted_momentum())
	//	{
	//	auto j = ij.second;
	//	cout << "sg jet pt : " << j->Pt() << endl;
	//	n++;
	//}

	n = 0;
	Int_t njetacc = 0;
	static Double_t rhovec[999];
 	 Double_t TotaljetAreaPhys=0;
  	Double_t TotalTPCArea=2*TMath::Pi()*0.9;
	for (auto iBg : ktContainer->accepted_momentum())
	{
		if (n < numberofexcludingjets) {
			n++;
			continue;
		}
		auto bgjet = iBg.second;

		Bool_t matched = false;
		for (auto iSg : aktContainer->accepted_momentum())
		{
			auto sgjet = iSg.second;
			matched = (isOverlapping(bgjet, sgjet)) ? true : false;
		}
		//cout << "n = "<<n<< " kt jet pt : " << bgjet->Pt() << " matched : "<< matched << " jet eta = "<<bgjet->Eta()<< endl;
		if (bgjet -> GetNumberOfTracks()>0){
			rhovec[njetacc] = bgjet->Pt() / bgjet->Area();
			TotaljetAreaPhys += bgjet->Area();
			njetacc++;
		}
			
		n++;
	}

  	Double_t OccCorr=1;
    OccCorr = TotaljetAreaPhys/TotalTPCArea;
	if (njetacc > 0 ) RHO = TMath::Median(njetacc, rhovec);
	else
		RHO = 0;

	RHO *= OccCorr;
	//cout << "jet pt end\n\n"
	//	 << endl;
}

//___________________________________________________________________
void AliBKJetAnalysis::MeasureBgDensity(AliJetContainer *ktContainer)
{
	using TMath::Abs;
	int n = 0;
	TLorentzVector1D rhoarray;
	Double1D Sumpt;
	Double1D Summ;
	Double1D Sumpx;
	Double1D Sumpy;
	Double1D Sumpz;
	Double1D Sume;
	TLorentzVector leadingkt;
	Bool_t issubleading = false;

	for (auto ij : ktContainer->accepted_momentum())
	{
		auto j = ij.second;
		if (fabs(j->Eta()) > 0.7)
			continue;
		double lpt = 0;
		double sumpt = 0;
		double summ = 0;
		TLorentzVector sumkt(0, 0, 0, 0);

		TLorentzVector avec(0, 0, 0, 0);
		avec.SetPtEtaPhiE(j->AreaPt(), j->AreaEta(), j->AreaPhi(), j->AreaE());
		for (int it = 0; it < j->GetNumberOfTracks(); it++)
		{
			auto trk = j->Track(it);
			if (trk->Pt() > 100)
				continue;
			//cout<< " isMC : " <<fIsMC<<" Bg Density TestFilterBit :  "<<(((AliAODTrack*) trk)->TestFilterBit(768))<<endl;
			//if (fOption.Contains("Emb") && !(((AliAODTrack *)trk)->TestFilterBit(768)))
			//	continue;
			//if (fOption.Contains("Emb")) cout<<"KT jet : " << " Label = " <<trk -> GetLabel()<< " FilterBit = " << ((AliAODTrack *)trk)->TestFilterBit(768)<<endl;
			TLorentzVector temp;
			temp.SetXYZM(trk->Px(), trk->Py(), trk->Pz(), pionmass);
			if (lpt < temp.Pt())
				lpt = temp.Pt();
			sumkt += temp;
			sumpt += temp.Pt();
			summ += (sqrt(pionmass * pionmass + temp.Pt() * temp.Pt()) - temp.Pt());
		}
		//if (lpt<5) continue;
		//cout<<"\n\n"<<endl;
		TLorentzVector avec2(0, 0, 0, 0);
		avec2.SetPtEtaPhiE(j->AreaPt(), j->AreaEta(), j->AreaPhi(), j->AreaE());
		sumpt = sumkt.Pt() / j->Area();
		summ = summ / j->Area();

		if (fabs(sumkt.DeltaPhi(leadingkt)) > pi / 2. && !issubleading)
		{
			issubleading = true;
		}
		//if (n !=0 || !issubleading){
		if (n > 1)
		{
			//Sumpt.push_back(sumpt);
			//Summ.push_back(summ);
			Sumpt.push_back(sumkt.Pt() / j->Area());
			Summ.push_back(sumkt.M() / j->Area());
			Sumpx.push_back(sumkt.Px() / avec2.Px());
			Sumpy.push_back(sumkt.Py() / avec2.Py());
			Sumpz.push_back(sumkt.Pz() / avec2.Pz());
			Sume.push_back(sumkt.E() / avec2.E());
		}
		n++;
	}

	Double_t rhopt[Sumpt.size()];
	Double_t rhom[Summ.size()];
	Double_t rhopx[Sumpx.size()];
	Double_t rhopy[Sumpy.size()];
	Double_t rhopz[Sumpz.size()];
	Double_t rhoe[Sume.size()];
	int i = 0;
	for (auto k : Sumpt)
	{
		rhopt[i] = k;
		i++;
	}
	i = 0;
	for (auto k : Summ)
	{
		rhom[i] = k;
		i++;
	}
	i = 0;
	for (auto k : Sumpx)
	{
		rhopx[i] = k;
		i++;
	}
	i = 0;
	for (auto k : Sumpy)
	{
		rhopy[i] = k;
		i++;
	}
	i = 0;
	for (auto k : Sumpz)
	{
		rhopz[i] = k;
		i++;
	}

	i = 0;
	for (auto k : Sume)
	{
		rhoe[i] = k;
		i++;
	}

	RHO = TMath::Median(Sumpt.size(), rhopt);
	RHOM = TMath::Median(Summ.size(), rhom);
	RHOPX = TMath::Median(Sumpx.size(), rhopx);
	RHOPY = TMath::Median(Sumpy.size(), rhopy);
	RHOPZ = TMath::Median(Sumpz.size(), rhopz);
	RHOE = TMath::Median(Sume.size(), rhoe);
	//FillTHnSparse("hRho", {fCent, RHO, pthardbin}, sf);
}

Bool_t AliBKJetAnalysis::MeasurePtHardBinScalingFactor()
{
	NTrials = -1;
	XSection = -1;
	genzvtx = -30;


	if (fIsMC)
	{
		if (fOption.Contains("AOD"))
		{
			AliVEvent *event = InputEvent();
			if (fOption.Contains("Emb"))
				event = AliAnalysisTaskEmcalEmbeddingHelper::GetInstance()->GetExternalEvent();
			if (!event)
			{
				Printf("ERROR: Could not retrieve event");
				sf = 0;
			}
			AliAODMCHeader *cHeaderAOD = dynamic_cast<AliAODMCHeader *>(event->FindListObject(AliAODMCHeader::StdBranchName()));
			genzvtx = cHeaderAOD->GetVtxZ();
			fHistos->FillTH1("mczvtx", genzvtx);
			TList *genHeaders = cHeaderAOD->GetCocktailHeaders();
			NTrials = -1;
			XSection = -1;
			AliGenEventHeader *gh = 0;
			for (Int_t i = 0; i < genHeaders->GetEntries(); i++)
			{
				gh = (AliGenEventHeader *)genHeaders->At(i);
				TString GeneratorName = gh->GetName();
				if (GeneratorName.CompareTo("Pythia") == 0)
				{
					AliGenPythiaEventHeader *gPythia = dynamic_cast<AliGenPythiaEventHeader *>(gh);
					NTrials = gPythia->Trials();
					XSection = gPythia->GetXsection();
					Double_t pthard = gPythia->GetPtHard();
					pthardbin = double(binningpthard.FindBin(pthard)) - 0.5;
					Int_t nTriggerJets = gPythia->NTriggerJets();
					TParticle *jet = nullptr;
					Float_t tmpjet[] = {0, 0, 0, 0};
					for (Int_t ijet = 0; ijet < nTriggerJets; ijet++)
					{
						gPythia->TriggerJet(ijet, tmpjet);
						jet = new TParticle(94, 21, -1, -1, -1, -1, tmpjet[0], tmpjet[1], tmpjet[2], tmpjet[3], 0, 0, 0, 0);
						if (jet->Pt() > 4.0 * pthard)
							return false;
					}
				}
			}
			sf = XSection / NTrials;

			fMCArray = (TClonesArray *)event->FindListObject("mcparticles");
			const Int_t nTracksMC = fMCArray->GetEntriesFast();
			for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++)
			{
				AliAODMCParticle *trackMC = dynamic_cast<AliAODMCParticle *>(fMCArray->At(iTracks));
				Int_t pdgCode = trackMC->PdgCode();
				if (iTracks == 4)
					p6.SetXYZT(trackMC->Px(), trackMC->Py(), trackMC->Pz(), trackMC->E());
				if (iTracks == 5)
					p7.SetXYZT(trackMC->Px(), trackMC->Py(), trackMC->Pz(), trackMC->E());
			}
			const UInt_t ntracks = event->GetNumberOfTracks();
			for (UInt_t it = 0; it < ntracks; it++)
			{
				auto t = dynamic_cast<AliAODTrack *>(event->GetTrack(it));
				if (!t)
					continue;
				//cout << t->TestFilterBit(768) << endl;
			}
			//cout<<"\n\n"<<endl;
		}
		else
		{
			AliMCEvent *mcEvent = MCEvent();
			AliStack *stack = mcEvent->Stack();
			AliGenPythiaEventHeader *gPythia = dynamic_cast<AliGenPythiaEventHeader *>(mcEvent->GenEventHeader());
			NTrials = gPythia->Trials();
			XSection = gPythia->GetXsection();
			sf = XSection / NTrials;
			Int_t nPrim = stack->GetNprimary();
			for (Int_t i = 0; i < nPrim; i++)
			{
				TParticle *p = stack->Particle(i);
				if (!p)
					continue;
			}
		}
	}
	//The method above doesn't work for MC 13b4_fix and plus AOD files
	//For these MC productions, the method below is used..
	if (fOption.Contains("LHC13"))
	//if (fOption.Contains("13plus") || fOption.Contains("12a15e") || fOption.Contains("12a15f") || fOption.Contains("13e4") || fOption.Contains("13fix"))
	{
		TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
		if (!tree)
			sf = 0;
		;
		TFile *curfile = tree->GetCurrentFile();
		if (!curfile)
			sf = 0;
		TString fCurrFileName = TString(curfile->GetName());
		if (!filename.EqualTo(fCurrFileName))
		{
			filename = fCurrFileName.Data();
			fCurrFileName.ReplaceAll(gSystem->BaseName(fCurrFileName.Data()), "");
			TFile *fxsec = TFile::Open(Form("%s%s", fCurrFileName.Data(), "pyxsec_hists.root"));
			TKey *key = (TKey *)fxsec->GetListOfKeys()->At(0);
			TList *list = dynamic_cast<TList *>(key->ReadObj());
			XSection = ((TProfile *)list->FindObject("h1Xsec"))->GetBinContent(1);
			NTrials = ((TH1F *)list->FindObject("h1Trials"))->GetBinContent(1);
			Int_t entries = ((TH1F *)list->FindObject("h1Trials"))->GetEntries();
			NTrials /= entries;
			pyxsechistsf = XSection / NTrials;
			fxsec->Close();
		}
		sf = pyxsechistsf;
	}
	return true;
}

Bool_t AliBKJetAnalysis::MeasureJets(AliJetContainer *jetContainer, TLorentzVector1D &Jets, TLorentzVector1D &JetsBeforeCorr, Bool_t istruth, Bool_t isfulljet)
{
	using TMath::Abs;

	//=== NOTE : We will not use sj[0] because bin in THnSparse is always begin with 1.
	for (auto ij : jetContainer->accepted_momentum())
	{
		auto j = ij.second;
		fHistos->FillTH2("jetetaphi", j->Eta(), j->Phi(), sf);
		TLorentzVector sum(0, 0, 0, 0);
		double sumpt = 0;
		double lpt = 0;

		for (int it = 0; it < j->GetNumberOfTracks(); it++)
		{
			auto trk = j->Track(it);

			if (istruth)
			{
				if (!((AliAODMCParticle *)trk->IsPhysicalPrimary()))
					continue;
				if (trk->Charge() == 0)
					continue;
			}
			if (!istruth && trk->Pt() > 100)
				continue;
			//if (fOption.Contains("Emb") && !istruth && !(((AliAODTrack *)trk)->TestFilterBit(768)))
			//continue;
			//if (fOption.Contains("Emb")) cout<<"Truth : "<< istruth << " Label = " <<trk -> GetLabel()<< " FilterBit = " << ((AliAODTrack *)trk)->TestFilterBit(768)<<endl;

			if (!istruth && fOption.Contains("MBTR"))
			{
				if (fBSRandom->Uniform(0, 100) < 5.)
					continue;
			}
			TLorentzVector temp;
			temp.SetXYZM(trk->Px(), trk->Py(), trk->Pz(), pionmass);
			sum += temp;
			sumpt += trk->Pt();
			if (lpt < temp.Pt())
				lpt = temp.Pt();
		}
		//cout<<"\n\n"<<endl;
		if (isfulljet)
		{
			for (int it = 0; it < j->GetNumberOfClusters(); it++)
			{
				auto clu = j->Cluster(it);
				TLorentzVector temp;
				clu->GetMomentum(temp, vertex);
				sum += temp;
				sumpt += temp.Pt();
				if (lpt < temp.Pt())
					lpt = temp.Pt();
			}
		}

		TLorentzVector avec(0, 0, 0, 0);
		avec.SetPtEtaPhiE(j->AreaPt(), j->AreaEta(), j->AreaPhi(), j->AreaE());
		TLorentzVector sumcorr(
			sum.Px() - RHOPX * avec.Px(), sum.Py() - RHOPY * avec.Py(), sum.Pz() - RHOPZ * avec.Pz(), sum.E() - RHOE * avec.E());

		if (fOption.Contains("MC") || fOption.Contains("Emb"))
		{
			//if (TMath::Abs(sumcorr.DeltaR(p6))<0.4 && (sumcorr.Pt()/p6.Pt())>4) return false;
			//if (TMath::Abs(sumcorr.DeltaR(p7))<0.4 && (sumcorr.Pt()/p7.Pt())>4) return false;
		}

		//if (lpt < fLeadingParticlePtMin)
		//	continue;

		//if (istruth)
		//{
		//	if (fabs(sum.Eta()) < 0.5)
		//Jets.push_back(sum);
		//}

		Double_t jetetacut = isfulljet ? 0.3  : 0.5 ; 
		if (istruth || !fIsAA)
		{
			if (fabs(sum.Eta()) < jetetacut)
			{
				Jets.push_back(sum);
			}
		}
		else
		{
			if (fabs(sumcorr.Eta()) < jetetacut)
			{
				Jets.push_back(sumcorr);
			}
		}

		if (fabs(sum.Eta()) < jetetacut)
		{
			JetsBeforeCorr.push_back(sum);
		}
	}
	return true;
}

void AliBKJetAnalysis::CheckDijetSelections(TLorentzVector1D Jets, TLorentzVector2D &sj, Bool1D &disel)
{
	using TMath::Abs;

	//=== DiJetSelection 1 : LS = Leading-SubLeading
	TLorentzVector1D zsj(2, TLorentzVector(0, 0, 0, 0));
	sj[0] = zsj;
	for (auto j : Jets)
	{
		if (sj[0][0].Pt() < j.Pt())
		{
			sj[0][0] = j;
		}
		else if (sj[0][1].Pt() < j.Pt())
		{
			sj[0][1] = j;
		}
	}

	//0 lorentz vector

	disel[kNoCut] = true;
	sj[kNoCut] = sj[0];
	Double_t kinecut = 20;

	//DiJetSelection 1 : TpTJET>20, ApTJET>20, pT_pair
	if (sj[0][0].Pt() > kinecut && sj[0][1].Pt() > kinecut)
	{
		sj[kB1] = sj[0];
		disel[kB1] = true;
	}
	//DiJetSelection 2 : TpTJET>20, ApTJET>20, pT_pair
	//Leading-Subreading & dPhi>pi/2
	if (sj[0][0].Pt() > kinecut && sj[0][1].Pt() > kinecut &&
		fabs(sj[0][0].DeltaPhi(sj[0][1])) > pi / 2.)
	{
		sj[kB2] = sj[0];
		disel[kB2] = true;
	}

	//DiJetSelection 3 : TpTJET>20, ApTJET>20, pT_pair
	//Leading - Subleading in opposite hemisphere
	if (sj[0][0].Pt() > kinecut)
	{
		auto tj = sj[kB3][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kB3][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(tj.DeltaPhi(jet)) < pi / 2.)
				continue;
			if (sj[kB3][1].Pt() < jet.Pt())
				sj[kB3][1] = jet;
		}
		if (sj[kB3][1].Pt() > kinecut)
			disel[kB3] = true;
	}

	if (sj[0][0].Pt() > 30 && sj[0][1].Pt() > 30)
	{
		sj[kC1] = sj[0];
		disel[kC1] = true;
	}
	//DiJetSelection 2 : TpTJET>20, ApTJET>20, pT_pair
	//Leading-Subreading & dPhi>pi/2
	if (sj[0][0].Pt() > 30 && sj[0][1].Pt() > 30 &&
		fabs(sj[0][0].DeltaPhi(sj[0][1])) > pi / 2.)
	{
		sj[kC2] = sj[0];
		disel[kC2] = true;
	}

	//DiJetSelection 3 : TpTJET>20, ApTJET>20, pT_pair
	//Leading - Subleading in opposite hemisphere
	if (sj[0][0].Pt() > 30)
	{
		auto tj = sj[kC3][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kC3][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(tj.DeltaPhi(jet)) < pi / 2.)
				continue;
			if (sj[kC3][1].Pt() < jet.Pt())
				sj[kC3][1] = jet;
		}
		if (sj[kC3][1].Pt() > 30)
			disel[kC3] = true;
	}

	if (sj[0][0].Pt() > 40 && sj[0][1].Pt() > 40)
	{
		sj[kD1] = sj[0];
		disel[kD1] = true;
	}
	//DiJetSelection 2 : TpTJET>20, ApTJET>20, pT_pair
	//Leading-Subreading & dPhi>pi/2
	if (sj[0][0].Pt() > 40 && sj[0][1].Pt() > 40 &&
		fabs(sj[0][0].DeltaPhi(sj[0][1])) > pi / 2.)
	{
		sj[kD2] = sj[0];
		disel[kD2] = true;
	}

	//DiJetSelection 3 : TpTJET>20, ApTJET>20, pT_pair
	//Leading - Subleading in opposite hemisphere
	if (sj[0][0].Pt() > 30)
	{
		auto tj = sj[kD3][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kD3][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(tj.DeltaPhi(jet)) < pi / 2.)
				continue;
			if (sj[kD3][1].Pt() < jet.Pt())
				sj[kD3][1] = jet;
		}
		if (sj[kD3][1].Pt() > 30)
			disel[kD3] = true;
	}

	//DiJetSelection 4 : 20<TpTJET<40, 20<ApTJET, kTy
	if (sj[0][0].Pt() > 20 && sj[0][0].Pt() < 40)
	{
		auto tj = sj[kM1][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kM1][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(TVector2::Phi_0_2pi(tj.DeltaPhi(jet)) - pi) > pi / 3.)
				continue;
			if (sj[kM1][1].Pt() < jet.Pt())
				sj[kM1][1] = jet;
		}
		if (sj[kM1][1].Pt() > 20)
			disel[kM1] = true;
	}
	else
		sj[kM1] = zsj;

	//DiJetSelection 5 : 40<TpTJET<60, 20<ApTJET, kTy
	if (sj[0][0].Pt() > 40 && sj[0][0].Pt() < 60)
	{
		auto tj = sj[kM2][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kM2][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(TVector2::Phi_0_2pi(tj.DeltaPhi(jet)) - pi) > pi / 3.)
				continue;
			if (sj[kM2][1].Pt() < jet.Pt())
				sj[kM2][1] = jet;
		}
		if (sj[kM2][1].Pt() > 20)
			disel[kM2] = true;
	}
	else
		sj[kM2] = zsj;

	//DiJetSelection 6 : 60<TpTJET<80, 20<ApTJET, kTy
	if (sj[0][0].Pt() > 60 && sj[0][0].Pt() < 80)
	{
		auto tj = sj[kM3][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kM3][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(TVector2::Phi_0_2pi(tj.DeltaPhi(jet)) - pi) > pi / 3.)
				continue;
			if (sj[kM3][1].Pt() < jet.Pt())
				sj[kM3][1] = jet;
		}
		if (sj[kM3][1].Pt() > 20)
			disel[kM3] = true;
	}
	else
		sj[kM3] = zsj;

	//DiJetSelection 7 : 80<TpTJET<110, 20<ApTJET, kTy
	if (sj[0][0].Pt() > 80 && sj[0][0].Pt() < 110)
	{
		auto tj = sj[kM4][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kM4][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(TVector2::Phi_0_2pi(tj.DeltaPhi(jet)) - pi) > pi / 3.)
				continue;
			if (sj[kM4][1].Pt() < jet.Pt())
				sj[kM4][1] = jet;
		}
		if (sj[kM4][1].Pt() > 20)
			disel[kM4] = true;
	}
	else
		sj[kM4] = zsj;

	//DiJetSelection 8 : 70<TpTJET<110, 20<ApTJET<30, kTy
	if (sj[0][0].Pt() > 70 && sj[0][0].Pt() < 110)
	{
		auto tj = sj[kM5][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kM5][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(TVector2::Phi_0_2pi(tj.DeltaPhi(jet)) - pi) > pi / 3.)
				continue;
			if (sj[kM5][1].Pt() < jet.Pt())
				sj[kM5][1] = jet;
		}
		if (sj[kM5][1].Pt() > 20 && sj[kM5][1].Pt() < 30)
			disel[kM5] = true;
	}
	else
		sj[kM5] = zsj;

	//DiJetSelection 9 : 70<TpTJET<110, 30<ApTJET<40, kTy

	if (sj[0][0].Pt() > 70 && sj[0][0].Pt() < 110)
	{
		auto tj = sj[kM6][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kM6][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(TVector2::Phi_0_2pi(tj.DeltaPhi(jet)) - pi) > pi / 3.)
				continue;
			if (sj[kM6][1].Pt() < jet.Pt())
				sj[kM6][1] = jet;
		}
		if (sj[kM6][1].Pt() > 30 && sj[kM6][1].Pt() < 40)
			disel[kM6] = true;
	}
	else
		sj[kM6] = zsj;

	//DiJetSelection 10 : 70<TpTJET<110, 40<ApTJET, kTy
	if (sj[0][0].Pt() > 70 && sj[0][0].Pt() < 110)
	{
		auto tj = sj[kM7][0] = sj[0][0]; // Leading is same as 1:LS
		sj[kM7][1] = TLorentzVector(0, 0, 0, 0);
		for (auto jet : Jets)
		{
			if (fabs(TVector2::Phi_0_2pi(tj.DeltaPhi(jet)) - pi) > pi / 3.)
				continue;
			if (sj[kM7][1].Pt() < jet.Pt())
				sj[kM7][1] = jet;
		}
		if (sj[kM7][1].Pt() > 40)
			disel[kM7] = true;
	}
	else
		sj[kM7] = zsj;
}

THnSparse *AliBKJetAnalysis::CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t *opt)
{
	const TAxis *axises[bins.size()];
	for (UInt_t i = 0; i < bins.size(); i++)
		axises[i] = &bins[i];
	THnSparse *h = fHistos->CreateTHnSparse(name, title, ndim, axises, opt);
	return h;
}

THnSparse *AliBKJetAnalysis::CreateTHnSparse(TString name, TString title, TString templ, Option_t *opt)
{
	auto o = fHistos->FindObject(templ);
	if (!o)
	{
		cout << "ERROR: no " << templ << endl;
		gSystem->Exit(1);
	}
	auto ht = dynamic_cast<THnSparse *>(o);
	const TAxis *axises[ht->GetNdimensions()];
	for (int i = 0; i < ht->GetNdimensions(); i++)
		axises[i] = ht->GetAxis(i);
	auto h = fHistos->CreateTHnSparse(name, title, ht->GetNdimensions(), axises, opt);
	return h;
}

Long64_t AliBKJetAnalysis::FillTHnSparse(TString name, std::vector<Double_t> x, Double_t w)
{
	auto hsparse = dynamic_cast<THnSparse *>(fHistos->FindObject(name));
	if (!hsparse)
	{
		cout << "ERROR : no " << name << endl;
		exit(1);
	}
	return FillTHnSparse(hsparse, x, w);
}

/*
   Long64_t AliBKJetAnalysis::Fill( TString name, std::vector<Double_t> x, Double_t w ){
   return FillTHnSparse( name, x, w );
   }
   */

Long64_t AliBKJetAnalysis::FillTHnSparse(THnSparse *h, std::vector<Double_t> x, Double_t w)
{
	if (int(x.size()) != h->GetNdimensions())
	{
		cout << "ERROR : wrong sized of array while Fill " << h->GetName() << endl;
		exit(1);
	}
	return h->Fill(&x.front(), w);
}

TAxis AliBKJetAnalysis::AxisFix(TString name, int nbin, Double_t xmin, Double_t xmax)
{
	TAxis axis(nbin, xmin, xmax);
	axis.SetName(name);
	return axis;
}

TAxis AliBKJetAnalysis::AxisStr(TString name, std::vector<TString> bin)
{
	TAxis ax = AxisFix(name, bin.size(), 0.5, bin.size() + 0.5);
	UInt_t i = 1;
	for (auto blabel : bin)
		ax.SetBinLabel(i++, blabel);
	return ax;
}

TAxis AliBKJetAnalysis::AxisVar(TString name, std::vector<Double_t> bin)
{
	TAxis axis(bin.size() - 1, &bin.front());
	axis.SetName(name);
	return axis;
}

TAxis AliBKJetAnalysis::AxisLog(TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0)
{
	int binoffset = (xmin0 < 0 || (xmin - xmin0) < 1e-9) ? 0 : 1;
	std::vector<Double_t> bin(nbin + 1 + binoffset, 0);
	double logBW3 = (log(xmax) - log(xmin)) / nbin;
	for (int ij = 0; ij <= nbin; ij++)
		bin[ij + binoffset] = xmin * exp(ij * logBW3);
	TAxis axis(nbin, &bin.front());
	axis.SetName(name);
	return axis;
}
Double_t AliBKJetAnalysis::getDiffR(double phi1, double phi2, double eta1, double eta2){
  Double_t diffPhi = TMath::Abs(phi1-phi2);
  if(diffPhi > TMath::Pi()){
    diffPhi = 2*TMath::Pi() - diffPhi;
  }
  return TMath::Sqrt(TMath::Power(diffPhi,2)+TMath::Power(eta1-eta2,2));
}