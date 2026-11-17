# Script generated on Tue Dec 12 21:34:26 2025
# Edited by Wenyuan Su

library(patRoon)

# -------------------------
# initialization
# -------------------------
workPath <- "E:/patRoon/PA"
setwd(workPath)

# Set paths for external software
options(patRoon.path.pwiz = "C:/ProteoWizard")                  # ProteoWizard installation folder
options(patRoon.path.SIRIUS = "D:/Program Files/sirius")       # Optional; only required for SIRIUS
patRoon::verifyDependencies()

# Generate analysis information for raw MS files
anaInfo <- generateAnalysisInfo(paths = "E:/patRoon/PA/analyses/raw")

# Convert Thermo raw files to centroided mzML files; set FALSE to skip this step
doDataPretreatment <- TRUE
if (doDataPretreatment) {
  convertMSFiles(anaInfo = anaInfo, outPath = "analyses/mzml", from = "thermo",
                 to = "mzML", algorithm = "pwiz", centroid = "vendor")
}

# Generate analysis information for the demo mzML dataset
anaInfo <- generateAnalysisInfo(paths = "E:/patRoon/PA/analyses/mzml",
                                groups = c("Pooled_sample", "BK-STD"),
                                blanks = c("BK-STD", "BK-STD"),
                                norm_concs = c(NA, NA))

# -------------------------
# feature detection
# -------------------------
# Detect chromatographic features using OpenMS
fList <- findFeatures(anaInfo, "openms", noiseThrInt = 10000)

# Group and align features across analyses
fGroups <- groupFeatures(fList, "openms", rtalign = TRUE)

# Rule-based filtering of low-intensity, poorly reproducible and blank-associated features
fGroupsfil <- filter(fGroups, preAbsMinIntensity = 100, absMinIntensity = 10000,
                     relMinReplicateAbundance = 1, maxReplicateIntRSD = 0.75,
                     blankThreshold = 3, removeBlanks = TRUE, retentionRange = NULL,
                     mzRange = c(100, 800))

# -------------------------
# componentization
# -------------------------
# Group isotopes and adducts belonging to the same compound
components <- generateComponents(fGroupsfil, "openms", ionization = "positive")

# Select preferred molecular ions and monoisotopic peaks
# prefAdduct can be changed according to ionization mode, e.g., [M]+, [M+H]+ or [M-H]-
fGroupsSel <- selectIons(fGroupsfil, components, prefAdduct = "[M+H]+", onlyMonoIso = TRUE)

# -------------------------
# suspect screening
# -------------------------
# This section can be skipped for a fully nontarget workflow
suspList <- read.csv("C:/Users/SWY/Desktop/suspectlist.csv", stringsAsFactors = FALSE)
fGroupsSusp <- screenSuspects(fGroupsSel, suspList, onlyHits = TRUE)

# Retain suspect matches within 2 ppm
si <- screenInfo(fGroupsSusp)

if (nrow(si) > 0) {
  reference_mass <- ifelse(is.na(si$mz), si$neutralMass, si$mz)
  
  keep <- is.finite(reference_mass) &
    is.finite(si$d_mz) &
    abs(si$d_mz) <= reference_mass * 2e-6
  
  fGroupsSusp@screenInfo <- si[keep]
}

fGroupsSusp <- filter(fGroupsSusp, onlyHits = TRUE)

# -------------------------
# transformation products
# -------------------------
# Examples of different approaches for TP generation; select according to the study purpose
TPs_CTS <- generateTPsCTS(suspList, "combined_photolysis_abiotic_hydrolysis", generations = 2)
TPs_library <- generateTPs("library", parents = suspList)
TPs_biotransformer <- generateTPs("biotransformer", type = "env", generations = 2, parents = suspList)

# Example: export transformation products generated from the internal library
tp_edges <- as.data.table(TPs_library)
data.table::fwrite(tp_edges, file = "library.tsv", sep = "\t")

# -------------------------
# MS/MS peak lists
# -------------------------
# Generate averaged MS peak lists; precursorMzWindow = NULL can be used for DIA data
avgMSListParams <- getDefAvgPListParams(clusterMzWindow = 0.005)
mslists <- generateMSPeakLists(fGroupsSusp, "mzr", maxMSRtWindow = 12,
                               precursorMzWindow = NULL, avgFeatParams = avgMSListParams,
                               avgFGroupParams = avgMSListParams)

# Remove weak MS/MS peaks and retain the 30 most intense fragments
mslists <- filter(mslists, absMSIntThr = NULL, absMSMSIntThr = NULL,
                  relMSIntThr = NULL, withMSMS = TRUE, relMSMSIntThr = 0.02,
                  topMSPeaks = NULL, topMSMSPeaks = 30)

# -------------------------
# diagnostic MS/MS filtering
# -------------------------
mz_tolerance_ppm <- 5

characteristic_fragments <- c(98.9842, 251.0468, 265.0624, 327.0781)
neutral_losses <- c(18.0106, 28.0313, 42.0470, 56.0626)

# Match within 5 ppm of the theoretical target mass
equalMZ <- function(target_mz, observed_mz) {
  tolerance_da <- target_mz * mz_tolerance_ppm * 1e-6
  any(abs(observed_mz - target_mz) <= tolerance_da, na.rm = TRUE)
}

# Match neutral losses within 5 ppm of each theoretical loss mass
findNeutralLosses <- function(observed_mz) {
  if (length(observed_mz) < 2) {
    return(rep(FALSE, length(neutral_losses)))
  }
  
  mass_diffs <- abs(combn(
    observed_mz, 2,
    FUN = function(x) x[1] - x[2]
  ))
  
  sapply(neutral_losses, function(loss) equalMZ(loss, mass_diffs))
}

# Retain spectra containing at least two diagnostic features
min_diagnostic_features <- 2

mslistsF <- delete(mslists, j = function(pl, grp, ana, type) {
  if (type != "MSMS") return(integer(0))
  
  num_characteristic_fragments <- sum(sapply(
    characteristic_fragments, equalMZ, observed_mz = pl$mz
  ))
  num_neutral_losses <- sum(findNeutralLosses(pl$mz))
  num_diagnostic_features <- num_characteristic_fragments +
    num_neutral_losses
  
  if (num_diagnostic_features < min_diagnostic_features) {
    return(seq_len(nrow(pl)))
  }
  
  return(integer(0))
})

# -------------------------
# formula and structural annotation
# -------------------------
# Generate molecular formula candidates using GenForm
formulas <- generateFormulas(fGroupsSusp, mslistsF, "genform", relMzDev = 2,
                             MSMode = "msms", batchSize = 64, maxCandidates = 1000,
                             elements = "C[100]H[200]O[10]N[10]P[3]S[3]F[10]Cl[10]Br[10]",
                             oc = TRUE, calculateFeatures = FALSE, timeout = 60)

# Generate structural candidates using MetFrag and PubChemLite
compounds <- generateCompounds(fGroupsSusp, mslistsF, "metfrag", method = "CL",
                               dbRelMzDev = 2, fragRelMzDev = 5, 
                               database = "pubchemlite", maxCandidatesToStop = 100)

# Combine suspect matches, formula candidates and structural annotations
fGroupsa <- annotateSuspects(fGroupsSusp, formulas = formulas, compounds = compounds,
                             MSPeakLists = mslistsF)

# Incorporate molecular formula scores into compound candidate ranking
compounds <- addFormulaScoring(compounds, formulas, updateScore = TRUE)

# -------------------------
# reporting
# -------------------------
# Advanced report settings can be modified in report.yml
report(fGroupsa, MSPeakLists = mslistsF, formulas = formulas, compounds = compounds,
       components = NULL, TPs = NULL, parallel = TRUE, openReport = FALSE)