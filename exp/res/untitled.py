#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os.path
import os
import sys
import pandas as pd

dataset_name = sys.argv[1]
timeout = int(sys.argv[2])

df = pd.read_csv(dataset_name)



print "Mean: " + str(df.mean())
