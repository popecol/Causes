# Set default species
message("Setting species variables...")
species <- "PP"
inter_species <- "PJ"

species_df <- read.csv2(file = "data/species.csv")
species_dict <- setNames(species_df$name, species_df$ID)
inter_species_name <- species_dict[inter_species]


#' Print selected species
#'
#' Helper function used to print selected species
#'
#' @param species_focal character; focal species
#' @param species_inter character; interacting species
#' @param data_dir character; path to a folder containing folders with species data
#' @param print_inter_species  logical; if `TRUE` interactive species if set
#'
#' @return
#' @export
#'
#' @examples
#'
#' set_species()
#' set_species(print_inter_species = TRUE)
species_info <- function(species_focal = species, species_inter = inter_species,
                        data_dir = "data/", print_inter_species = FALSE,
                        dict_species = species_dict) {

  available_species <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

  # print info
  cat(paste("Available species:", paste(available_species, collapse = ", "),
            sep = " "), "\n")
  cat(paste0("Selected species: ", species_focal, " (", dict_species[species_focal] ,")"), "\n")

  if(print_inter_species)
    cat(paste0("Selected interacting species: ", species_inter, " (", dict_species[species_inter] ,")"), "\n")

  message("\nNote: You can change selected species in set_species.R and use source(\"R/set_species.R\") to apply changes")

  }
