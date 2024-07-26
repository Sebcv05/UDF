#!/bin/bash
#Bash script to get latest version of UDF 

#Step 0: Remove old UDF files and initialize UDF

rm -r build
rm CMakeL*
rm -r src
rm -r build
rm -r include
rm libconverge_udf.so

source $HOME/3.0.25.sh
cvg_udf_init_3025


# Variables
HOME_DIR=$HOME
STR1="/Desktop/Remote/db/UDF/"    # Path to your local repository
REPO_DIR=$HOME$STR1
SRC="/src"
INC="/include"
WORKING_DIR=$(pwd)                       # Path to your working directory

# Step 1: Navigate to the repository directory
cd $REPO_DIR

# Step 2: Pull the latest changes from the repository
echo "Pulling the latest changes from the repository..."
git pull origin main


# Step 3: Copy the latest code to the working directory
echo "Copying the latest code to the working directory..."

rsync -av --exclude='.git' "$REPO_DIR$SRC" "$WORKING_DIR" || { echo "Failed to copy files"; exit 1; }
rsync -av --exclude='.git' "$REPO_DIR$INC" "$WORKING_DIR" || { echo "Failed to copy files"; exit 1; }

# Step 4: Navigate to the working directory
cd "$WORKING_DIR" || { echo "Failed to navigate to working directory"; exit 1; }



# Step 5: Compile the UDF 
echo "Compiling the code..."
cd build
cmake ..
make
cp libconverge_udf.so ..
cd ..

# Step 6: Indicate that the script has finished
echo "Update and compilation complete."