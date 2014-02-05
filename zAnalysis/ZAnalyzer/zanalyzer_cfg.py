import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYZE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/home/users/dsklein/ryan_lessons/MYDEMOANALYZER/CMSSW_5_3_2_patch4/src/zAnalysis/ZEventProducer/zbabies.root'
    )
)

process.zanalyzer = cms.EDAnalyzer('ZAnalyzer'
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('zAnalysis_plots.root')
                                   )

process.p = cms.Path(process.zanalyzer)
