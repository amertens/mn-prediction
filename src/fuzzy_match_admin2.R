# Load required libraries
library(dplyr)
library(stringdist)
library(RecordLinkage)

# Function to clean district names
clean_district_names <- function(names) {
  names %>%
    # Remove \r\n line breaks
    gsub("\\r\\n", " ", .) %>%
    # Remove common suffixes
    gsub("\\s+(Municipal|Metropolis|Metropolitan)$", "", .) %>%
    # Remove extra whitespace
    trimws() %>%
    # Convert to title case for consistency
    stringr::str_to_title() %>%
    # Replace hyphens with spaces for better matching
    gsub("-", " ", .)
}


# Clean both datasets
data1_clean <- data.frame(
  original_district = unique(data1$Admin2),
  clean_district = clean_district_names(unique(data1$Admin2)),
  stringsAsFactors = FALSE
)

# Clean both datasets
data2_clean <- data.frame(
  original_district = unique(data2$District),
  clean_district = clean_district_names(unique(data2$District)),
  stringsAsFactors = FALSE
)

df_clean <- data.frame(
  original_admin2 = unique(df$Admin2),
  clean_admin2 = clean_district_names(unique(df$Admin2)),
  stringsAsFactors = FALSE
)

# Function to perform fuzzy matching
fuzzy_match_districts <- function(data2_clean, df_clean, threshold = 0.8) {

  # Create distance matrix using Jaro-Winkler similarity
  dist_matrix <- stringdistmatrix(
    data2_clean$clean_district,
    df_clean$clean_admin2,
    method = "jw"  # Jaro-Winkler
  )

  # Convert to similarity (1 - distance)
  sim_matrix <- 1 - dist_matrix

  # Find best matches
  matches <- data.frame(
    data2_district = character(),
    df_admin2 = character(),
    similarity = numeric(),
    match_type = character(),
    stringsAsFactors = FALSE
  )

  for(i in 1:nrow(data2_clean)) {
    # Find the best match for each district in data2
    best_match_idx <- which.max(sim_matrix[i, ])
    best_similarity <- sim_matrix[i, best_match_idx]

    match_type <- case_when(
      best_similarity >= 0.95 ~ "Exact/Very High",
      best_similarity >= threshold ~ "Good Match",
      TRUE ~ "Poor Match"
    )

    matches <- rbind(matches, data.frame(
      data2_district = data2_clean$original_district[i],
      df_admin2 = df_clean$original_admin2[best_match_idx],
      data2_clean = data2_clean$clean_district[i],
      df_clean = df_clean$clean_admin2[best_match_idx],
      similarity = best_similarity,
      match_type = match_type,
      stringsAsFactors = FALSE
    ))
  }

  return(matches)
}

# Perform the matching
matches <- fuzzy_match_districts(data2_clean, df_clean, threshold = 0.8)

matches1 <- fuzzy_match_districts(data1_clean, df_clean, threshold = 0.8)

# Sort by similarity (best matches first)
matches <- matches[order(matches$similarity, decreasing = TRUE), ]

# View results
print("=== MATCHING RESULTS ===")
print(paste("Total districts in data2:", nrow(data2_clean)))
print(paste("Total admin2 regions in df:", nrow(df_clean)))
print(paste("Good matches (>=0.8):", sum(matches$similarity >= 0.8)))
print(paste("Poor matches (<0.8):", sum(matches$similarity < 0.8)))

# Show summary by match type
print("\n=== MATCH QUALITY SUMMARY ===")
print(table(matches$match_type))

# Display first 20 matches
print("\n=== FIRST 20 MATCHES ===")
print(matches[1:20, c("data2_district", "df_admin2", "similarity", "match_type")])

# Show poor matches that need manual review
poor_matches <- matches[matches$similarity < 0.8, ]
if(nrow(poor_matches) > 0) {
  print("\n=== POOR MATCHES NEEDING MANUAL REVIEW ===")
  print(poor_matches[, c("data2_district", "df_admin2", "similarity")])
}

# Create a lookup table for joining datasets
lookup_table <- matches %>%
  select(data2_district, df_admin2, similarity, match_type) %>%
  rename(
    District = data2_district,
    Admin2 = df_admin2
  )

# Example of how to use the lookup table to join your datasets
# Assuming your original datasets are 'data2' and 'df'
join_example <- function() {
  # Add the matched Admin2 to data2
  data2_matched <- data2 %>%
    left_join(lookup_table, by = "District")

  # Now you can join with df using the matched Admin2
  combined_data <- data2_matched %>%
    left_join(df, by = "Admin2")

  return(combined_data)
}

# # Alternative: Create a manual correction list for problematic matches
# create_manual_corrections <- function(poor_matches) {
#   manual_corrections <- data.frame(
#     District = poor_matches$data2_district,
#     Suggested_Admin2 = poor_matches$df_admin2,
#     Action_Needed = "REVIEW",
#     stringsAsFactors = FALSE
#   )
#   return(manual_corrections)
# }
#
# create_manual_corrections(poor_matches)
#
# # Export results for manual review if needed
# if(nrow(poor_matches) > 0) {
#   manual_review <- create_manual_corrections(poor_matches)
#   print("\n=== EXPORT MANUAL REVIEW TABLE ===")
#   print("Use this to manually correct poor matches:")
#   print(manual_review)
# }
#
# # Function to apply manual corrections
# apply_manual_corrections <- function(lookup_table, corrections) {
#   # corrections should be a data.frame with columns: District, Corrected_Admin2
#   corrected_lookup <- lookup_table
#   for(i in 1:nrow(corrections)) {
#     idx <- which(corrected_lookup$District == corrections$District[i])
#     if(length(idx) > 0) {
#       corrected_lookup$Admin2[idx] <- corrections$Corrected_Admin2[i]
#       corrected_lookup$match_type[idx] <- "Manual Correction"
#     }
#   }
#   return(corrected_lookup)
# }
