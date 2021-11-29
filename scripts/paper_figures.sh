#!/bin/bash

#
# Generate all the paper figures, one after the other.
#

# Tissue visualisation
./scripts/fig_expression.sh

# Chemoaffinity model, random choice of gradient expression. WT.
./scripts/fig_G0_wildtype.sh

# Chemoaffinity model, hand-tuned gradient expression. WT.
./scripts/fig_G_wildtype.sh

# Chemoaffinity model, hand-tuned gradient expression. Surgical manipulations
./scripts/fig_G_surgical.sh

# GC model, Surgical manipulations
./scripts/fig_GC_surgical.sh

# GI model, Surgical manipulations
#./scripts/fig_GI_surgical.sh

# GJ model, Surgical manipulations
#./scripts/fig_GJ_surgical.sh
