######################################
### Section 2 - DOWNLOAD Tmax DATA ###
######################################

# input  data: "data/stations.rda"
# output data: "data/Y365.rda"

### DATA ACCESSED FOR THIS WORK ON 
### 23-SEPTEMBER-2022

### Note: ECA makes periodic updates,
###       so the data used in this work
###       and those accessed at another
###       time may differ

# Clear workspace
rm(list = ls())

# Read list of stations of interest. Based on previous work
load("data/stations.rda")

# Download Tmax dataset from ECA
# Increase timeout time
getOption("timeout")
options(timeout = 2400)

# Around 600 MB are required
tmp_dir <- tempdir()
if (!file.exists(tmp_dir)) dir.create(path = tmp_dir, recursive = T)
tmp_file <- file.path(tmp_dir, "ECA_blend_tx.zip")
download.file(url = "https://knmi-ecad-assets-prd.s3.amazonaws.com/download/ECA_blend_tx.zip",
              destfile = tmp_file)

# Get list of zipped files
zip_ls <- unzip(tmp_file, list = TRUE)$Name
head(zip_ls)
# Get number ID
num_ls <- grep("TX_STAID", zip_ls)
num_id <- sapply(num_ls, function(x) unlist(strsplit(zip_ls[x], "TX_STAID")[[1]][2]))
num_id <- as.integer(sapply(1:length(num_id), function(x) unlist(strsplit(num_id[x], ".txt")[[1]][1])))

# Find the ones already in the stations data.frame and extract them
raw_dir <- file.path(tmp_dir, "Raw")
if (!dir.exists(raw_dir)) dir.create(path = raw_dir, recursive = TRUE)
# Extract files of interest
# Info file
unzip(zipfile = tmp_file,
      files = zip_ls[num_ls[match(stations$STAID, num_id)]],
      exdir = raw_dir)

#Extract all the files in the folder Iberia_TX from that correspond to the data from each station
tx_ls <- list.files(path = raw_dir, full.names = TRUE)

# Create empty data.frame
# Daily date range
target_date <- seq(from = as.Date("1960-01-01"),
                   to = as.Date("2021-12-31"),
                   by = "day")
tail(target_date)
target_date <- target_date[grep("-02-29", target_date, invert = TRUE)]
# Output file
tx_mat <- as.data.frame(matrix(data = as.numeric(NA),
                               nrow = length(target_date),
                               ncol = nrow(stations)))
rownames(tx_mat) <- target_date
colnames(tx_mat) <- stations$STAID

# Get data
for (ii in 1:length(tx_ls)){
  cat(paste0("..", ii))
  # Read file
  aux <- read.csv(file = tx_ls[ii], header = TRUE, skip = 20)
  # Remove scores without good label data (0)
  nozero <- which(aux$Q_TX != 0)
  if (length(nozero) > 0) aux$TX[nozero] <- NA
  # Take the days of interest
  aux$DATE <- as.Date(strptime(aux$DATE, "%Y%m%d"))
  tx_mat[, ii] <- aux$TX[match(as.Date(rownames(tx_mat)), aux$DATE)]
}

# Reshape data
LL <- 365
TT <- nrow(tx_mat) / LL
SS <- nrow(stations)
# Transpose 't' and 'l'
Y365 <- apply(array(data = as.matrix(tx_mat), dim = c(LL, TT, SS)), c(1, 3), t) / 10

# Write object
save(Y365, file = "data/Y365.rda", compress = "xz")

# Delete temporary directory
unlink(tmp_dir, recursive = TRUE)
