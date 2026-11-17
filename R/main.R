#' Encode DNA 5-mer Sequences for Prediction
#'
#' @description
#' Converts DNA sequence strings into nucleotide-position categorical factors.
#' Each nucleotide position is represented as a factor with levels A, T, C, and G,
#' matching the encoding used to train the m6A random forest model.
#'
#' @param dna_strings A character vector containing DNA sequences of equal length.
#'
#' @return A data frame in which each column represents one nucleotide position.
#' Each column is a factor with levels c("A", "T", "C", "G").
#'
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "GAACA"))
#'
#' @keywords internal
#'
#' @import randomForest
#'
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Sites for Multiple Samples
#'
#' @description
#' Predicts m6A-positive probabilities and binary m6A statuses for multiple
#' samples using a fitted random forest model. The input data frame must contain
#' the seven features used by the model, including a DNA 5-mer sequence.
#'
#' @param ml_fit A fitted random forest classification model.
#' @param feature_df A data frame containing the columns gc_content, RNA_type,
#' RNA_region, exon_length, distance_to_junction, evolutionary_conservation,
#' and DNA_5mer.
#' @param positive_threshold A numeric probability threshold used to classify a
#' sample as Positive. The default is 0.5.
#'
#' @return The input data frame with two additional columns:
#' predicted_m6A_prob and predicted_m6A_status.
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction"))
#' example_data <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                      package = "m6APrediction"))
#' head(prediction_multiple(ml_fit, example_data, positive_threshold = 0.6))
#'
#' @export
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df)
}

#' Predict an m6A Site for a Single Sample
#'
#' @description
#' Predicts the m6A-positive probability and binary m6A status for one sample
#' from seven input features using a fitted random forest model.
#'
#' @param ml_fit A fitted random forest classification model.
#' @param gc_content A numeric value for GC content.
#' @param RNA_type A character value specifying RNA type: mRNA, lincRNA,
#' lncRNA, or pseudogene.
#' @param RNA_region A character value specifying RNA region: CDS, intron,
#' 3'UTR, or 5'UTR.
#' @param exon_length A numeric value for exon length.
#' @param distance_to_junction A numeric value for distance to junction.
#' @param evolutionary_conservation A numeric value for evolutionary conservation.
#' @param DNA_5mer A five-base DNA sequence represented by A, T, C, and G.
#' @param positive_threshold A numeric probability threshold used to classify a
#' sample as Positive. The default is 0.5.
#'
#' @return A named character vector containing predicted_m6A_prob and
#' predicted_m6A_status.
#'
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction"))
#' prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#'
#' @export
#'
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector)
}
