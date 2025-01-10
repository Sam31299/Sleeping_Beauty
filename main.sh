#!/bin/bash

# Create results directories if they don't exist
mkdir -p results/circos
mkdir -p results/lollipop
mkdir -p results/table

# Echo start of analysis
echo "Starting CIS analysis pipeline..."
echo "================================"

# Function to run R script and check status
run_r_script() {
    local script=$1
    echo "Running $script..."
    Rscript "lib/$script"
    if [ $? -eq 0 ]; then
        echo "✓ Successfully completed $script"
    else
        echo "✗ Error running $script"
        exit 1
    fi
    echo "--------------------------------"
}

# Run each analysis script in sequence
echo "Running table analysis..."
run_r_script "table.R"

echo "Running co-occurrence analysis..."
run_r_script "SB-co-occurrence.R"

echo "Running circos analysis..."
run_r_script "circos.R"

echo "Running lollipop analysis..."
run_r_script "lollipop.R"

echo "================================"
echo "Analysis pipeline complete!"