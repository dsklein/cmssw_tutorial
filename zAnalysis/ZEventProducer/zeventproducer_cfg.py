import FWCore.ParameterSet.Config as cms

process = cms.Process("ZANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/nfs-7/userdata/edm/53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_AODSIM_PU_S10_START53_V7A-v1.root'
    )
)

process.zProducer = cms.EDProducer('ZEventProducer'
)

process.tauFilter = cms.EDFilter('TauFilter'
)

process.out = cms.OutputModule("PoolOutputModule",
							   fileName = cms.untracked.string('zbabies.root'),
							   SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      "keep *_zProducer_*_*"
                                                                      )
)


  
process.p = cms.Path(process.tauFilter*process.zProducer)

process.e = cms.EndPath(process.out)
