AliBKJetAnalysis *AddTaskBKJetAnalysis(const char *taskname, const char *option, Int_t isaa)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    return NULL;
  }

  AliEmcalCorrectionTask *correctionTask = AddTaskEmcalCorrectionTask();
  correctionTask->SelectCollisionCandidates(AliVEvent::kAny);
  //correctionTask->SetUserConfigurationFilename("AliEmcalCorrectionConfiguration.yaml");
  correctionTask->SetUserConfigurationFilename("$ALICE_PHYSICS/PWG/EMCAL/config/AliEmcalCorrectionConfiguration.yaml");
  //correctionTask->SetUserConfigurationFilename("PWGJESampleConfig.yaml");
  correctionTask->Initialize();

  AliEmcalJetTask *jetFinderTask;
  AliEmcalJetTask *jetFinderTaskFullJet;
  AliEmcalJetTask *jetFinderTaskkt;
  AliEmcalJetTask *jetFinderTaskkine;
  AliEmcalJetTask *jetFinderTaskkineFullJet;
  AliEmcalJetTask *jetFinderTaskktFullJet;
  //TString trackcontname = option.Contains("Emb") ? "PicoTracksMer" : "tracks";
  //TString trackcontname  = "tracks";
  //TString mccontname="mcparticles";
  TString Option(option);
  jetFinderTask = AddTaskEmcalJet("usedefault", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.300, 0.005, AliJetContainer::pt_scheme, Form("Jets%s", Option.Data()), 0.);
  jetFinderTaskkt = AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, 0.2, AliJetContainer::kChargedJet, 0.15, 0.300, 0.005, AliJetContainer::pt_scheme, Form("JetsKt%s", Option.Data()), 0.);
  jetFinderTaskFullJet = AddTaskEmcalJet("usedefault", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kFullJet, 0.15, 0.300, 0.005, AliJetContainer::pt_scheme, Form("FullJets%s", Option.Data()), 0.);
  jetFinderTaskktFullJet = AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, 0.2, AliJetContainer::kFullJet, 0.15, 0.300, 0.005, AliJetContainer::pt_scheme, Form("JetsKt%s", Option.Data()), 0.);

  if (Option.Contains("Emb") || Option.Contains("MC"))
  {
    jetFinderTaskkine = AddTaskEmcalJet("mcparticles", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kChargedJet, 0.15, 0.300, 0.005, AliJetContainer::pt_scheme, Form("Jetsmc%s", Option.Data()), 0.);
    jetFinderTaskkineFullJet = AddTaskEmcalJet("mcparticles", "", AliJetContainer::antikt_algorithm, 0.4, AliJetContainer::kFullJet, 0.15, 0.300, 0.005, AliJetContainer::pt_scheme, Form("Jetsmc%s", Option.Data()), 0.);
  }

  AliAnalysisTaskRhoSparse *rhosparse = AddTaskRhoSparse("usedefault", "usedefault", "Rho", 0.2, AliEmcalJet::kTPCfid, AliJetContainer::kChargedJet, AliJetContainer::pt_scheme, kTRUE, "", "TPC", 0.15, 0.01, 0, "");
  if (Option.Contains("LHC13d") || Option.Contains("LHC13e"))
    rhosparse->SelectCollisionCandidates(AliVEvent::kEMCEJE);
  else
    rhosparse->SelectCollisionCandidates(AliVEvent::kINT7);

  // Own class
  AliBKJetAnalysis *task = new AliBKJetAnalysis(taskname, option);
  task->SetIsAA(isaa);
  task->SetLeadingParticlePtMin(5);
  //task->SelectCollisionCandidates( Option.Contains("MC") ? AliVEvent::kAny : AliVEvent::kINT7 );

  AliJetContainer *jetCont = task->AddJetContainer(jetFinderTask->GetName(), "TPCfid");
  AliJetContainer *jetContKt = task->AddJetContainer(jetFinderTaskkt->GetName(), "TPCfid");
  AliJetContainer *jetContFullJet = task->AddJetContainer(jetFinderTask->GetName(), "EMCALfid");
  AliJetContainer *jetContKtFullJet = task->AddJetContainer(jetFinderTaskktFullJet->GetName(), "EMCALfid");

  AliTrackContainer *trackCont = task->AddTrackContainer("usedefault");
  trackCont->SetFilterHybridTracks(true);
  AliClusterContainer *clusterCont = task->AddClusterContainer("usedefault");

  if (Option.Contains("Emb"))
  {
    trackCont->SetIsEmbedding(true);
    jetFinderTask->AdoptParticleContainer(trackCont);
    jetFinderTaskkt->AdoptParticleContainer(trackCont);
  }
  jetCont->ConnectParticleContainer(trackCont);
  jetCont->SetRhoName("Rho");
  jetContKt->ConnectParticleContainer(trackCont);
  jetContKtFullJet->ConnectParticleContainer(trackCont);
  jetContFullJet->ConnectParticleContainer(trackCont);
  jetContFullJet->ConnectClusterContainer(clusterCont);

  //jetCont->SetLeadingHadronType( 0 );

  AliMCParticleContainer *mcCont = 0x0;
  if (Option.Contains("Emb") || Option.Contains("MC"))
  {
    mcCont = task->AddMCParticleContainer("mcparticles");
    AliJetContainer *jetContMC = task->AddJetContainer(jetFinderTaskkine->GetName(), "TPC");
    AliJetContainer *jetContMCFullJet = task->AddJetContainer(jetFinderTaskkineFullJet->GetName(), "EMCALfid");
    jetContMC->ConnectParticleContainer(mcCont);
    jetContMCFullJet->ConnectParticleContainer(mcCont);
  }

  //jetCont->SetZLeadingCut( 0.98, 0.98 ); // FIXME: Comments me and others
  cout << "jetFinderTask->GetRadius() : " << jetFinderTask->GetRadius() << endl;
  if (jetFinderTask->GetRadius() >= 0.4)
    jetCont->SetPercAreaCut(0.6);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("%s%s", taskname, option), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s", AliAnalysisManager::GetCommonFileName()));
  if (Option.Contains("Emb") || Option.Contains("MC"))
    task->SetIsMC(true);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  return task;
}
