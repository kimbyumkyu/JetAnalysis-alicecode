//------------------------------------------------------------------
// When there is no time to wait for legotrain results,
// please use this macro to run jobs on grid.
// Recommended for pp or pA runs, not for AA runs.
//------------------------------------------------------------------
// Author: Beomkyu Kim
// email:  kimb@cern.ch
//
//#include "runlist.h"
//#include "AliBKJetAnalysis.h"

R__ADD_INCLUDE_PATH($ALICE_PHYSICS/OADB/macros)
#include "AddTaskPhysicsSelection.C"
#include "AddTaskCentrality.C"
R__ADD_INCLUDE_PATH($ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/)
#include "AddTaskMultSelection.C"
R__ADD_INCLUDE_PATH($ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/)
#include "AddTaskEmcalJet.C"
R__ADD_INCLUDE_PATH($ALICE_PHYSICS/PWG/EMCAL/macros/)
#include "AddTaskEmcalCorrectionTask.C"
//#include "AliBKJetAnalysis.h"
//R__LOAD_LIBRARY(AliBKJetAnalysis.cxx)
//#include "AddTaskBSDiJet.C"

class AliAnalysisManager;
class AliEmcalJetTask;
class AliBKJetAnalysis;
AliAnalysisManager *mgr = 0x0;
//AliBKJetAnalysis * AddTaskBSDiJet(TString taskname, TString option );
void run(
	const char *taskname = "Dijet", const char *option = "LHC17pAOD" // when scanning AOD, add "AOD"
	,
	const char *gridmode = "full" // or "terminate" to merge
	,
	const char *location = "local")
{
	TString foption = option;
	// add aliroot indlude path
	gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ROOTSYS")));
	gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_ROOT")));
	gROOT->ProcessLine(Form(".include %s/include", gSystem->ExpandPathName("$ALICE_PHYSICS")));

	// analysis manager
	mgr = new AliAnalysisManager(Form("%s%s", taskname, option));
	// create the alien handler and attach it to the manager
	AliAnalysisAlien *plugin = new AliAnalysisAlien();
	plugin->SetRunMode(gridmode);
	plugin->SetAPIVersion("V1.1x");
	plugin->SetAliPhysicsVersion("vAN-20200425_ROOT6-1");
	plugin->SetDropToShell(0);
	if (!foption.Contains("MC"))
		plugin->SetRunPrefix("000");
	plugin->SetNrunsPerMaster(1);
	plugin->SetOutputToRunNo();
	plugin->SetMergeViaJDL(1);

	bool isAA = false;

	//const int LHC17pRuns[] = {282343, 282342, 282341, 282340, 282314, 282313, 282312, 282309, 282307, 282306, 282305, 282304, 282303, 282302, 282247, 282230, 282229, 282227, 282224, 282206, 282189, 282147, 282146, 282127, 282126, 282125,
	//						  282123, 282122, 282120, 282119, 282118, 282099, 282098, 282078, 282051, 282050, 282031, 282030, 282025, 282021, 282016, 282008};
	const int LHC17pRuns[] = {282127};
	const int LHC15oRuns[] = {246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757,
							  246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053,
							  246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245692, 245683};

	const int LHC15nRuns[]={244628, 244627, 244626, 244619, 244618, 244617, 244542, 244540, 244531, 244484, 244483, 244482, 244481, 244480, 244456, 244453, 244377};
	//const int LHC15nRuns[]={244531};
	const int LHC13cRuns[]={195677, 195675, 195673, 195644, 195635, 195633, 195596, 195593, 195592, 195568, 195567, 195566, 195531, 195529};

	//const int LHC13dRuns[] = {195872, 195867, 195831, 195829, 195827, 195826, 195787, 195783, 195767, 195760, 195724};
	const int LHC13dRuns[] = {195872};

	if (foption.Contains("LHC15n")){
		if (!foption.Contains("MC")) {
			plugin->SetGridDataDir("/alice/data/2015/LHC15n/");
			for (int i=0; i<sizeof(LHC15nRuns)/sizeof(LHC15nRuns[0]); i++)
				plugin->AddRunNumber(LHC15nRuns[i]);
			plugin->SetDataPattern("/pass4/AOD208/*/AliAOD.root");
			isAA = false;
			plugin->SetSplitMaxInputFileNumber(300);
		}
		else {
			//plugin->AddDataFile("/alice/cern.ch/user/k/kimb/16h3.xml");
			plugin->AddDataFile("/alice/cern.ch/user/k/kimb/244531.xml");
			isAA = false;
			plugin->SetSplitMaxInputFileNumber(25);
		}
	}

	if (foption.Contains("LHC13c"))
	{
		plugin->SetGridDataDir("/alice/data/2013/LHC13c/");
		for (int i = 0; i < sizeof(LHC13cRuns) / sizeof(LHC13cRuns[0]); i++)
			//for (int i=0; i<1; i++)
			plugin->AddRunNumber(LHC13cRuns[i]);
		plugin->SetDataPattern("/pass4/AOD210/*/AliAOD.root");
		isAA = true;
	}
	if (foption.Contains("LHC13d"))
	{
		plugin->SetGridDataDir("/alice/data/2013/LHC13d/");
		for (int i = 0; i < sizeof(LHC13dRuns) / sizeof(LHC13dRuns[0]); i++)
			//for (int i=0; i<1; i++)
			plugin->AddRunNumber(LHC13dRuns[i]);
		plugin->SetDataPattern("/pass2/AOD154/*/AliAOD.root");
		isAA = true;
	}

	if (foption.Contains("LHC15o"))
	{
		if (!foption.Contains("MC"))
		{
			plugin->SetGridDataDir("/alice/data/2015/LHC15o/");
			for (int i = 10; i < 20; i++)
			{
				plugin->AddRunNumber(LHC15oRuns[i]);
			}
			plugin->SetDataPattern("/pass1/AOD194/*/AliAOD.root");
			isAA = true;
			plugin->SetSplitMaxInputFileNumber(300);
		}
		else {
			//alien_find -x LHC15oMC16j5.xml  /alice/sim/2016/LHC16j5/  /AOD200/*/AliAOD.root  > LHC15oMC16j5.xml
			plugin->AddDataFile("/alice/cern.ch/user/k/kimb/LHC15oMC16j5.xml");
			isAA = false;
			plugin->SetSplitMaxInputFileNumber(25);
		}
	}

	if (foption.Contains("LHC17p"))
	{
		if (!foption.Contains("MC"))
		{
			plugin->SetGridDataDir("/alice/data/2017/LHC17p/");
			for (int i = 0; i < sizeof(LHC17pRuns) / sizeof(LHC17pRuns[0]); i++)
			{
				plugin->AddRunNumber(LHC17pRuns[i]);
			}
			plugin->SetDataPattern("/pass1_FAST/AOD208/*/AliAOD.root");
			isAA = false;
			plugin->SetSplitMaxInputFileNumber(20);
		}
		else
		{
			if (0)
			{
				plugin->SetGridDataDir("/alice/sim/2018/LHC18j2_fast/");
				for (int i = 0; i < sizeof(LHC17pRuns) / sizeof(LHC17pRuns[0]); i++)
				{
					plugin->AddRunNumber(LHC17pRuns[i]);
				}
				//plugin->AddRunNumber(LHC17pMCRuns);
				plugin->SetDataPattern("/AOD209/*/AliAOD.root");
				isAA = false;
				plugin->SetSplitMaxInputFileNumber(20);
			}

			// alien_find -x 28234all.xml /alice/sim/2018/LHC18b8_fast/ 28234*/AOD/*/AliAOD.root  > 28234all.xml
			//plugin->AddDataFile("/alice/cern.ch/user/k/kimb/28234all.xml");

			if (1)
			{
				plugin->AddDataFile("/alice/cern.ch/user/k/kimb/28234all.xml");
				isAA = false;
				plugin->SetSplitMaxInputFileNumber(25);
			}
		}
	}

	plugin->SetGridWorkingDir(Form("%s%s", taskname, option));
	plugin->SetGridOutputDir("out");
	plugin->AddIncludePath("-I$ALICE_ROOT/include  -I$ALICE_ROOT/lib -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/lib -I$ALICE_PHYSICS/OADB/macros");
	plugin->SetAnalysisSource("AliBKJetAnalysis.cxx");
	plugin->SetAdditionalLibs("libqpythia.so libAliPythia6.so libTree.so libGeom.so libVMC.so libPhysics.so libMinuit.so libGui.so libXMLParser.so libMinuit2.so libProof.so libSTEERBase.so libESD.so libAOD.so libOADB.so libANALYSIS.so libANALYSISalice.so libCDB.so libRAWDatabase.so libSTEER.so libCORRFW.so libEMCALUtils.so libPHOSUtils.so libJETAN.so libEMCALraw.so libEMCALbase.so libEMCALrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libPWGLFforward2.so libPWGTools.so libPWGEMCALtasks.so libCGAL.so libfastjet.so libsiscone.so libsiscone_spherical.so libfastjetplugins.so libfastjettools.so libfastjetcontribfragile.so  AliBKJetAnalysis.cxx AliBKJetAnalysis.h");
	plugin->SetDefaultOutputs(kFALSE);
	plugin->SetOutputArchive();
	//plugin->SetOutputFiles("AnalysisResults.root RecTree.root");
	plugin->SetOutputFiles("AnalysisResults.root stderr stdout");
	plugin->SetMasterResubmitThreshold(90);
	plugin->SetFileForTestMode("file.text");

	// Optionally set time to live (default 30000 sec)
	plugin->SetTTL(50000);
	// Optionally set input format (default xml-single)
	plugin->SetInputFormat("xml-single");
	// Optionally modify the name of the generated JDL (default analysis.jdl)
	plugin->SetJDLName(Form("%s%s.jdl", taskname, option));
	// Optionally modify the executable name (default analysis.sh)
	plugin->SetExecutable(Form("%s%s.sh", taskname, option));
	// Optionally modify job price (default 1)
	plugin->SetPrice(1);
	// Optionally modify split mode (default 'se')
	plugin->SetSplitMode("se");

	mgr->SetGridHandler(plugin);
	AliInputEventHandler *handler;
	if (foption.Contains("AOD"))
		handler = new AliAODInputHandler();
	else
		handler = new AliESDInputHandler();

	handler->SetNeedField(1);
	mgr->SetInputEventHandler(handler);

	//__________________________________________________
	//embedding

	Int_t pthardbin = 1;
	if (foption.Contains("Emb"))
	{
		AliAnalysisTaskEmcalEmbeddingHelper *embeddingHelper = AliAnalysisTaskEmcalEmbeddingHelper::AddTaskEmcalEmbeddingHelper();
		//if (foption.Contains("LHC15o")) gSystem->Exec(Form("alien_find /alice/sim/2016/LHC16j5/ AliAOD.root | grep AOD200 | perl -nle'print \"alien://\".$_' > embfile.txt"));
		embeddingHelper->SetFileListFilename("./embfile.txt");
		embeddingHelper->SetAOD();
		embeddingHelper->SetRandomFileAccess(kTRUE);
		embeddingHelper->Initialize();
		AliAnalysisTaskBSEmbedding *preprocess = AliAnalysisTaskBSEmbedding::AddTaskBSEmbedding();
	}

	AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(foption.Contains("MC") ? true : false);
	AliMultSelectionTask *multtask = AddTaskMultSelection(false); 


	gInterpreter->LoadMacro("AliBKJetAnalysis.cxx+g");

	if (foption.Contains("Emb"))
	{
		reinterpret_cast<AliBKJetAnalysis *>(gInterpreter->ExecuteMacro(
			Form("AddTaskBKJetAnalysis.C(\"%s\",\"%s\",%d)",taskname,option,isAA)));
		reinterpret_cast<AliBKJetAnalysis *>(gInterpreter->ExecuteMacro(
			Form("AddTaskBKJetAnalysis.C(\"%s\",\"LHC15oAOD\",%d)",taskname,isAA)));
	}
	else
	{
		reinterpret_cast<AliBKJetAnalysis *>(gInterpreter->ExecuteMacro(
			Form("AddTaskBKJetAnalysis.C(\"%s\",\"%s\",%d)", taskname, option,isAA)));
	}


	mgr->SetDebugLevel(2);
	if (!mgr->InitAnalysis())
		return;
	mgr->PrintStatus();

	// start analysis
	Printf("Starting Analysis....");
	mgr->StartAnalysis(location, 1234567890, 0);
}
