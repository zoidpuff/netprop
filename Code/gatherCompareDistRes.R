# Define a recursive function to insert an element into the nested list
insertIntoNestedList <- function(lst, keys, value) {
  if(length(keys) == 1) {
    lst[[keys]] <- value
  } else {
    if(!is.list(lst[[keys[1]]])) {
      lst[[keys[1]]] <- list()
    }
    lst[[keys[1]]] <- insertIntoNestedList(lst[[keys[1]]], keys[-1], value)
  }
  return(lst)
}

# Folder path containing the .rdata files
folderPath <- "/home/gummi/netprop/compareDist"

# Initialize the master list
masterRes <- list()

# List all .rdata files in the folder
files <- list.files(path = folderPath, pattern = "\\.rdata$", full.names = TRUE)

# Iterate over the files and load the 'res' object into the master list
for (file in files) {
  # Load the .rdata file
  load(file)
  
  # Extract the filename without extension and split by '_'
  filename <- tools::file_path_sans_ext(basename(file))
  parts <- unlist(strsplit(filename, "_"))
  
  # Remove the first component (common prefix) and any other not needed parts
  # In this example, we keep everything as is, assuming 'res' is what we want and is present in the .rdata file
  keys <- parts[-1] # Adjust if necessary
  
  # Insert the 'res' object into the nested list
  masterRes <- insertIntoNestedList(masterRes, keys, res)
}

save(masterRes, file = "masterRes.rdata")

# At this point, masterRes is structured as per your filename components
# To access an element, use something like:
# masterRes[["AllTraitsOverallScores"]][["noNorm"]][["0"]][["TRUE"]][["FALSE"]][["cosineDist"]]


