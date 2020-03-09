## Goal: Generate datasets {Ed10, Ed10n, Ed50, Ed50n, Ed100, Ed100n}
## along with the corresponding time-varying networks.
## How to run:
## > source('data-gen-script.R')
##
#=====================================================================================================================
##  Ed10, Ed10n gen
#=====================================================================================================================
## Remove all R objects
rm(list=ls())

source('gen-edi.R')

## Loads a template of obj 'edi.net'
load('asset/edi.net.template.RData')

## Loads obj 'true.net.adj.matrix'
load('asset/DREAM3GoldStandard_InSilicoSize10_Yeast1_TrueNet.RData')


## Generate synthetic time-varying networks.
GenLargerEdiNet(edi.net, true.net.adj.matrix, 2, 'asset/edi.net.10.RData', 'asset/edi.net.10.adj.mx.RData')

## Loads obj 'edi.net'
load('asset/edi.net.10.RData')

## Generate synthetic time-series gene expression data.
## Dataset ID: 'Ed10'.
## V = 10, T = 21, S = 4.
## No noise.
## 4 time series since Ds10 has 4 time series.
##
GenEdiData(edi.net, 0, 4, 'asset/edi-data-10.tsv')

## Generate synthetic time-series gene expression data.
## Dataset ID: 'Ed10n'.
## V = 10, T = 21, S = 4.
## Noise N(0, 0.05) like in the DREAM3 datasets.
## 4 time series since Ds10n has 4 time series.
##
GenEdiData(edi.net, 0.05, 4, 'asset/edi-data-10n.tsv')

## Remove all R objects
rm(list=ls())

#=====================================================================================================================
#=====================================================================================================================
##  Ed50, Ed50n gen
#=====================================================================================================================
source('gen-edi.R')

## Loads a template of obj 'edi.net'
load('asset/edi.net.template.RData')

## Loads obj 'true.net.adj.matrix'
load('asset/DREAM3GoldStandard_InSilicoSize50_Yeast1_TrueNet.RData')

## Generate synthetic time-varying networks.
GenLargerEdiNet(edi.net, true.net.adj.matrix, 2, 'asset/edi.net.50.RData', 'asset/edi.net.50.adj.mx.RData')

## Loads obj 'edi.net'
load('asset/edi.net.50.RData')

## Generate synthetic time-series gene expression data.
## Dataset ID: 'Ed50'.
## V = 50, T = 21, S = 23.
## No noise.
## 23 time series since Ds50 has 23 time series.
##
GenEdiData(edi.net, 0, 23, 'asset/edi-data-50.tsv')

## Generate synthetic time-series gene expression data.
## Dataset ID: 'Ed50n'.
## V = 50, T = 21, S = 23.
## Noise N(0, 0.05) like in the DREAM3 datasets.
## 23 time series since Ds50n has 23 time series.
##
GenEdiData(edi.net, 0.05, 23, 'asset/edi-data-50n.tsv')

## Remove all R objects
rm(list=ls())

#=====================================================================================================================
#=====================================================================================================================
##  Ed100, Ed100n gen
#=====================================================================================================================
  source('gen-edi.R')

## Loads a template of obj 'edi.net'
load('asset/edi.net.template.RData')

## Loads obj 'true.net.adj.matrix'
load('asset/DREAM3GoldStandard_InSilicoSize100_Yeast1_TrueNet.RData')

## Generate synthetic time-varying networks.
GenLargerEdiNet(edi.net, true.net.adj.matrix, 2, 'asset/edi.net.100.RData', 'asset/edi.net.100.adj.mx.RData')

## Loads obj 'edi.net'
load('asset/edi.net.100.RData')

## Generate synthetic time-series gene expression data.
## Dataset ID: 'Ed100'.
## V = 100, T = 21, S = 46.
## No noise.
## 46 time series since Ds100 has 46 time series.
##
GenEdiData(edi.net, 0, 46, 'asset/edi-data-100.tsv')

## Generate synthetic time-series gene expression data.
## Dataset ID: 'Ed100n'.
## V = 100, T = 21, S = 46.
## Noise N(0, 0.05) like in the DREAM3 datasets.
## 46 time series since Ds100n has 46 time series.
##
GenEdiData(edi.net, 0.05, 46, 'asset/edi-data-100n.tsv')

## Remove all R objects
rm(list=ls())

#=====================================================================================================================