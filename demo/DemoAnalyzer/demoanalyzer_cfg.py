import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("FWCore.MessageService.MessageLogger_cfi")
# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
	    limit = cms.untracked.int32(-1)
		)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/nfs-7/userdata/edm/53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_AODSIM_PU_S10_START53_V7A-v1.root'
        )
                            )

process.demo = cms.EDAnalyzer('DemoAnalyzer',
                              minTracks = cms.untracked.uint32(0)
                              )

#process.load("demo.DemoAnalyzer.demoanalyzer_cfi")
#process.demo.minTracks=1000

#process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histodemo.root')
                                   )

process.p = cms.Path(process.demo)
#process.p = cms.Path(process.demo*process.dump)
