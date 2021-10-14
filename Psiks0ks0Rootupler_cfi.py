import FWCore.ParameterSet.Config as cms

rootuple = cms.EDAnalyzer('JPsiKs0Ks0',
                          dimuons = cms.InputTag("slimmedMuons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          secondaryVerticesPtr = cms.InputTag("slimmedKshortVertices"),
                          bslabel = cms.InputTag("offlineBeamSpot"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          prescales     = cms.InputTag("patTrigger"),
                          objects       = cms.InputTag("slimmedPatTrigger"),
                          TriggerNames  = cms.vstring(
                                                      "HLT_Dimuon25_Jpsi_v",
					            #  "HLT_Dimuon20_Jpsi_Barrel_Seagulls_v",
                                                      "HLT_DoubleMu4_3_Jpsi_Displaced_v",
                                                      "HLT_DoubleMu4_JpsiTrk_Displaced_v",
							"HLT_DoubleMu4_JpsiTrkTrk_Displaced_v"),
                          PuInfoTag      = cms.InputTag("slimmedAddPileupInfo"),
                          GenParticles = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
			  # GenParticles = cms.InputTag("genParticles"),
                          #packedGenParticles = cms.InputTag("genParticles"),

                          OnlyBest = cms.bool(False),
                          isMC = cms.bool(True),
                          OnlyGen = cms.bool(False),
                          )
