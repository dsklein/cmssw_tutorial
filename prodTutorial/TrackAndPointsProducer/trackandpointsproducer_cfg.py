import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/nfs-7/userdata/edm/53X/RelValTTbar_GEN-SIM-RECO_CMSSW_5_3_2_patch1-START53_V7A.root'
    )
)

#process.myProducerLabel = cms.EDProducer('TrackAndPointsProducer')

process.MuonTrackPoints = cms.EDProducer('TrackAndPointsProducer', src = cms.InputTag('globalMuons') )

process.TrackTrackPoints = cms.EDProducer('TrackAndPointsProducer', src = cms.InputTag('generalTracks') )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('myOutputFile.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      "keep *_generalTracks_*_*",
                                                                      "keep *_globalMuons_*_*",
                                                                      "keep *_MuonTrackPoints_*_*",
                                                                      "keep *_TrackTrackPoints_*_*")
)

  
#process.p = cms.Path(process.myProducerLabel)
process.p = cms.Path(process.MuonTrackPoints*process.TrackTrackPoints)

process.e = cms.EndPath(process.out)
