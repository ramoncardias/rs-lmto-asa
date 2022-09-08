#!/bin/bash

# (1) will generate `output`
../../new_source/rs.x 

# (2) will generate `output-bulk`
../../new_source/rs.x input.nml output=output-bulk atoms=atoms-bulk.nml

# (3) will generate `output-self`
../../new_source/rs.x input.nml output=output-self self=self-different.nml

# (4) will generate `output-bulk-self`
../../new_source/rs.x input.nml output=output-bulk-self atoms=atoms-bulk.nml self=self-different.nml